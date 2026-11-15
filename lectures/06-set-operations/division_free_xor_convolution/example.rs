// XOR convolution without division: computing C_U = sum_{S xor T = U} A_S B_T
// for length-2^N arrays using only ring operations (add, subtract,
// multiply) -- no multiplication by an inverse of 2 or any other element,
// unlike the standard Walsh-Hadamard-transform solution. This matters
// whenever the coefficient ring has no inverse for 2 (e.g. Z/2^k Z, used
// below specifically because 2 is a zero divisor there and so provably has
// no inverse), or more generally whenever "no division" is a requirement.
//
// The approach: view A, B as multilinear polynomials in commuting
// variables x_0,...,x_{N-1} (A(x) = sum_S A_S * prod_{i in S} x_i). XOR
// convolution is exactly polynomial multiplication mod (x_i^2 - 1, all i).
// Substituting x_i = y_i - 1 turns each relation into y_i^2 - 2y_i, and
// multiplying mod (y_i^2 - 2y_i, all i) is a "ranked OR-convolution": the
// same machinery used for subset convolution, tracking |S|+|T| as an
// extra polynomial variable t so that y_i^2 = 2y_i can be applied via
// integer powers of 2 (an ordinary ring element, needing no division) at
// the end rather than during the convolution itself. Converting between
// the x- and y-bases is a per-variable change of basis, computable with
// N passes of a SOS-transform-shaped loop.
//
// All arithmetic below is mod 2^20 -- a power of two specifically because
// 2 is then a zero divisor (2 * 2^19 = 0) and therefore has no
// multiplicative inverse, so a Walsh-Hadamard-transform-based solution
// (which divides by 2 up to N times) could not be used in this ring at
// all, while the algorithm here works unmodified.

// PDFBEGIN
const MOD_BITS: u32 = 20;

const MOD: i64 = 1 << MOD_BITS;

const MASK: i64 = MOD - 1;

fn add(a: i64, b: i64) -> i64 {
    (a + b) & MASK
}

fn sub(a: i64, b: i64) -> i64 {
    (a - b) & MASK
}

fn mul(a: i64, b: i64) -> i64 {
    (a & MASK).wrapping_mul(b & MASK) & MASK
}

/// Subset zeta transform (SOS): B_T = sum_{S subseteq T} A_S. In place, O(n*2^n).
fn subset_zeta(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 1 {
                a[t] = add(a[t], a[t ^ (1 << i)]);
            }
        }
    }
    a
}

/// Subset Mobius transform, the inverse of the above. In place, O(n*2^n).
fn subset_mobius(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 1 {
                a[t] = sub(a[t], a[t ^ (1 << i)]);
            }
        }
    }
    a
}

/// Converts the coefficients of a multilinear polynomial from the x-basis
/// (A_S = coefficient of prod_{i in S} x_i) to the y-basis under x_i = y_i - 1,
/// i.e. A(x) = A1(y). Per variable: a0 + a1*x_i = a0 + a1*(y_i-1)
/// = (a0-a1) + a1*y_i, so (a0,a1) -> (a0-a1, a1). O(n*2^n).
fn transform_x_to_y(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 0 {
                let hi = t | (1 << i);
                a[t] = sub(a[t], a[hi]);
            }
        }
    }
    a
}

/// The inverse of `transform_x_to_y`: converts from the y-basis back to the
/// x-basis under y_i = x_i + 1, i.e. C1(y) = C(x). Per variable:
/// b0 + b1*y_i = b0 + b1*(x_i+1) = (b0+b1) + b1*x_i, so (b0,b1) -> (b0+b1, b1).
fn transform_y_to_x(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 0 {
                let hi = t | (1 << i);
                a[t] = add(a[t], a[hi]);
            }
        }
    }
    a
}

/// Problem 2: C_U = sum_{S union T = U} A_S * B_T * 2^(|S|+|T|-|U|), i.e.
/// polynomial multiplication mod (y_i^2 - 2y_i, all i). Solved exactly like
/// subset convolution, but keeping every rank r=|S|+|T| (not just r=|U|) so
/// that the final 2^(|S|+|T|-|U|) factor can be folded in afterward using
/// only integer powers of 2 (repeated multiplication -- no division). O(n^2*2^n).
fn ranked_or_convolution_pow2(n: usize, a1: &[i64], b1: &[i64]) -> Vec<i64> {
    let size = 1usize << n;

    // Split each array by popcount into n+1 "rank slices," then take the
    // subset zeta transform of each slice independently.
    let mut za: Vec<Vec<i64>> = vec![vec![0i64; size]; n + 1];
    let mut zb: Vec<Vec<i64>> = vec![vec![0i64; size]; n + 1];
    for s in 0..size {
        let k = (s as u32).count_ones() as usize;
        za[k][s] = a1[s];
        zb[k][s] = b1[s];
    }
    for k in 0..=n {
        za[k] = subset_zeta(n, std::mem::take(&mut za[k]));
        zb[k] = subset_zeta(n, std::mem::take(&mut zb[k]));
    }

    // At each U, ordinary-polynomial-multiply the two rank-indexed
    // sequences (za[.][U], zb[.][U]) to get conv[U][r] = sum_{k+l=r}
    // za[k][U]*zb[l][U] = sum_{S,T subseteq U} A1_S*B1_T*[|S|+|T|=r].
    let mut conv = vec![vec![0i64; 2 * n + 1]; size];
    for u in 0..size {
        for k in 0..=n {
            let zak = za[k][u];
            if zak == 0 {
                continue;
            }
            for l in 0..=n {
                let zbl = zb[l][u];
                if zbl == 0 {
                    continue;
                }
                conv[u][k + l] = add(conv[u][k + l], mul(zak, zbl));
            }
        }
    }

    // For each rank r, subset-Mobius-transform the column conv[.][r] to
    // recover chat[r][U] = sum_{S union T = U, |S|+|T|=r} A1_S*B1_T.
    let mut chat = vec![vec![0i64; size]; 2 * n + 1];
    for r in 0..=2 * n {
        let col: Vec<i64> = (0..size).map(|u| conv[u][r]).collect();
        chat[r] = subset_mobius(n, col);
    }

    // Fold in 2^(r - |U|) (always a nonnegative power, since |S|+|T| >= |S union T|).
    let mut c1 = vec![0i64; size];
    for u in 0..size {
        let pu = (u as u32).count_ones() as usize;
        let mut acc = 0i64;
        let mut pow2 = 1i64;
        for r in pu..=2 * n {
            acc = add(acc, mul(chat[r][u], pow2));
            pow2 = mul(pow2, 2);
        }
        c1[u] = acc;
    }
    c1
}

/// The full division-free XOR convolution algorithm.
fn xor_convolution_division_free(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let a1 = transform_x_to_y(n, a.to_vec());
    let b1 = transform_x_to_y(n, b.to_vec());
    let c1 = ranked_or_convolution_pow2(n, &a1, &b1);
    transform_y_to_x(n, c1)
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

/// Brute-force O(4^n) reference for XOR convolution.
fn xor_convolution_naive(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut c = vec![0i64; size];
    for s in 0..size {
        for t in 0..size {
            c[s ^ t] = add(c[s ^ t], mul(a[s], b[t]));
        }
    }
    c
}

/// Brute-force O(4^n) reference for Problem 2 (used to test
/// `ranked_or_convolution_pow2` in isolation).
fn problem2_naive(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut c = vec![0i64; size];
    for s in 0..size {
        for t in 0..size {
            let u = s | t;
            let w = (s as u32).count_ones() + (t as u32).count_ones() - (u as u32).count_ones();
            let mut pow2 = 1i64;
            for _ in 0..w {
                pow2 = mul(pow2, 2);
            }
            c[u] = add(c[u], mul(mul(a[s], b[t]), pow2));
        }
    }
    c
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
}

fn main() {
    let mut rng = Rng(2026091901);

    // --- Sanity check: transform_y_to_x undoes transform_x_to_y exactly ---
    let mut roundtrip_ok = true;
    for _ in 0..100 {
        let n = rng.range(0, 8) as usize;
        let size = 1usize << n;
        let a: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let y = transform_x_to_y(n, a.clone());
        let back = transform_y_to_x(n, y);
        if back != a {
            roundtrip_ok = false;
            println!("BASIS ROUND-TRIP MISMATCH n={n}");
        }
    }
    println!("Basis change x<->y round-trips exactly, over 100 random trials: {}", if roundtrip_ok { "PASS" } else { "FAIL" });

    // --- Problem 2 in isolation: ranked OR-convolution vs. brute force ---
    let mut problem2_ok = true;
    for _ in 0..100 {
        let n = rng.range(0, 6) as usize;
        let size = 1usize << n;
        let a: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let b: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let via_ranked = ranked_or_convolution_pow2(n, &a, &b);
        let via_naive = problem2_naive(n, &a, &b);
        if via_ranked != via_naive {
            problem2_ok = false;
            println!("PROBLEM2 MISMATCH n={n}");
        }
    }
    println!("Problem 2 (ranked OR-convolution with 2^(|S|+|T|-|U|) weights) matches brute force, over 100 random trials: {}", if problem2_ok { "PASS" } else { "FAIL" });

    // --- Full pipeline: division-free XOR convolution vs. brute force ---
    let mut full_ok = true;
    for _ in 0..100 {
        let n = rng.range(0, 6) as usize;
        let size = 1usize << n;
        let a: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let b: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let via_division_free = xor_convolution_division_free(n, &a, &b);
        let via_naive = xor_convolution_naive(n, &a, &b);
        if via_division_free != via_naive {
            full_ok = false;
            println!("XOR CONVOLUTION MISMATCH n={n}");
        }
    }
    println!(
        "Division-free XOR convolution (mod 2^{MOD_BITS}, where 2 is a zero divisor with no inverse) matches brute force, over 100 random trials: {}",
        if full_ok { "PASS" } else { "FAIL" }
    );

    let all_ok = roundtrip_ok && problem2_ok && full_ok;
    println!("\nAll checks: {}", if all_ok { "PASS" } else { "FAIL" });
}
