// Kth term of a linearly recurrent sequence: two O(M(d) log K)-time
// algorithms, both reducing the problem to fast exponentiation of a
// degree-d object. Fiduccia works with x^K modulo the recurrence's
// characteristic polynomial; Bostan-Mori works with the generating
// function's rational form P(x)/Q(x) directly, halving the target index
// each round via an even/odd split. Polynomial arithmetic here is the
// naive O(d^2) kind, so the overall cost is O(d^2 log K), matching what
// the lecture describes as the straightforward (non-FFT) implementation.

// PDFBEGIN
const MOD: i64 = 998244353;

fn add(a: i64, b: i64) -> i64 {
    (a + b) % MOD
}

fn mul(a: i64, b: i64) -> i64 {
    (a % MOD) * (b % MOD) % MOD
}

fn sub(a: i64, b: i64) -> i64 {
    ((a - b) % MOD + MOD) % MOD
}

// ---------------------------------------------------------------------
// Fiduccia's algorithm
// ---------------------------------------------------------------------

/// Multiply a degree-<d polynomial by x, reducing modulo the
/// characteristic polynomial Gamma(x) = x^d - c[0]x^{d-1} - ... - c[d-1],
/// using x^d = sum_j c[j] x^{d-1-j}. This single step is O(d).
fn mul_x_mod_gamma(poly: &[i64], c: &[i64]) -> Vec<i64> {
    let d = poly.len();
    let carry = poly[d - 1];
    let mut out = vec![0i64; d];
    for i in (1..d).rev() {
        out[i] = poly[i - 1];
    }
    for j in 0..d {
        out[d - 1 - j] = add(out[d - 1 - j], mul(carry, c[j]));
    }
    out
}

/// Multiply two degree-<d polynomials and reduce modulo Gamma(x), by
/// convolving then repeatedly substituting x^e (e >= d) with the
/// recurrence's own right-hand side, from the top degree down.
fn mulmod_gamma(a: &[i64], b: &[i64], c: &[i64]) -> Vec<i64> {
    let d = a.len();
    let mut raw = vec![0i64; 2 * d - 1];
    for i in 0..d {
        if a[i] == 0 {
            continue;
        }
        for j in 0..d {
            raw[i + j] = add(raw[i + j], mul(a[i], b[j]));
        }
    }
    for e in (d..raw.len()).rev() {
        let coeff = raw[e];
        if coeff == 0 {
            continue;
        }
        raw[e] = 0;
        for j in 0..d {
            raw[e - 1 - j] = add(raw[e - 1 - j], mul(coeff, c[j]));
        }
    }
    raw.truncate(d);
    raw
}

/// x^k mod Gamma(x), via the doubling recursion in which the "multiply by
/// x" step (for odd k) uses the cheap O(d) helper above.
fn xk_mod_gamma(k: u64, d: usize, c: &[i64]) -> Vec<i64> {
    if (k as usize) < d {
        let mut r = vec![0i64; d];
        r[k as usize] = 1;
        return r;
    }
    let half = xk_mod_gamma(k / 2, d, c);
    let mut r = mulmod_gamma(&half, &half, c);
    if k % 2 == 1 {
        r = mul_x_mod_gamma(&r, c);
    }
    r
}

/// A_K for a sequence satisfying A_i = c[0]A_{i-1} + ... + c[d-1]A_{i-d},
/// given initial values A_0..A_{d-1}.
fn fiduccia(init: &[i64], c: &[i64], k: u64) -> i64 {
    let d = init.len();
    if (k as usize) < d {
        return init[k as usize];
    }
    let r = xk_mod_gamma(k, d, c);
    (0..d).fold(0i64, |acc, i| add(acc, mul(r[i], init[i])))
}

// ---------------------------------------------------------------------
// Bostan-Mori's algorithm
// ---------------------------------------------------------------------

fn poly_mul(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0i64; a.len() + b.len() - 1];
    for i in 0..a.len() {
        if a[i] == 0 {
            continue;
        }
        for j in 0..b.len() {
            out[i + j] = add(out[i + j], mul(a[i], b[j]));
        }
    }
    out
}

fn negate_odd(a: &[i64]) -> Vec<i64> {
    a.iter().enumerate().map(|(i, &v)| if i % 2 == 1 { sub(0, v) } else { v }).collect()
}

fn even_part(a: &[i64]) -> Vec<i64> {
    a.iter().step_by(2).cloned().collect()
}

fn odd_part(a: &[i64]) -> Vec<i64> {
    if a.len() < 2 {
        return vec![0];
    }
    a[1..].iter().step_by(2).cloned().collect()
}

/// [x^k] P(x)/Q(x), given Q(0) = 1.
fn bostan_mori(mut p: Vec<i64>, mut q: Vec<i64>, mut k: u64) -> i64 {
    while k > 0 {
        let qneg = negate_odd(&q);
        let u = poly_mul(&p, &qneg);
        let w = poly_mul(&q, &qneg);
        p = if k % 2 == 0 { even_part(&u) } else { odd_part(&u) };
        q = even_part(&w);
        k /= 2;
    }
    p.first().cloned().unwrap_or(0)
}

/// Same recurrence as `fiduccia`, solved via the generating function
/// A(x) = P(x)/Q(x) with Q(x) = 1 - c[0]x - ... - c[d-1]x^d.
fn kth_term_bostan_mori(init: &[i64], c: &[i64], k: u64) -> i64 {
    let d = init.len();
    let mut q = vec![0i64; d + 1];
    q[0] = 1;
    for j in 0..d {
        q[j + 1] = sub(0, c[j]);
    }
    let full = poly_mul(init, &q);
    let p: Vec<i64> = full[..d].to_vec();
    bostan_mori(p, q, k)
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

// ---------------------------------------------------------------------
// Direct O(K) unrolling, for cross-checking on small K.
// ---------------------------------------------------------------------

fn direct(init: &[i64], c: &[i64], k: u64) -> i64 {
    let d = init.len();
    if (k as usize) < d {
        return init[k as usize];
    }
    let mut a = init.to_vec();
    for i in d..=(k as usize) {
        let mut next = 0i64;
        for j in 0..d {
            next = add(next, mul(c[j], a[i - 1 - j]));
        }
        a.push(next);
    }
    a[k as usize]
}

fn main() {
    // --- The lecture's own worked example: A_i = A_{i-1}+A_{i-2}+A_{i-3},
    // A_0=3, A_1=2, A_2=1. Expected A_5 = 16.
    let init = vec![3i64, 2, 1];
    let c = vec![1i64, 1, 1];
    let a5_direct = direct(&init, &c, 5);
    let a5_fid = fiduccia(&init, &c, 5);
    let a5_bm = kth_term_bostan_mori(&init, &c, 5);
    println!("worked example: A_5 = {a5_direct} (direct), {a5_fid} (Fiduccia), {a5_bm} (Bostan-Mori)");
    assert_eq!(a5_direct, 16);
    assert_eq!(a5_fid, 16);
    assert_eq!(a5_bm, 16);

    // --- Fibonacci, cross-checked against direct computation -----------
    let fib_init = vec![0i64, 1];
    let fib_c = vec![1i64, 1];
    let mut fib_ok = true;
    for k in 0..60u64 {
        let d = direct(&fib_init, &fib_c, k);
        let f = fiduccia(&fib_init, &fib_c, k);
        let b = kth_term_bostan_mori(&fib_init, &fib_c, k);
        if d != f || d != b {
            fib_ok = false;
            println!("FIBONACCI MISMATCH k={k}: direct={d} fiduccia={f} bostan_mori={b}");
        }
    }
    println!("Fibonacci cross-check (direct vs. both fast algorithms): {}", if fib_ok { "PASS" } else { "FAIL" });

    // --- Random small-order recurrences, many random K -------------------
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
    }
    let mut rng = Rng(42);
    let mut random_ok = true;
    for _ in 0..300 {
        let d = 1 + (rng.next() % 6) as usize;
        let init: Vec<i64> = (0..d).map(|_| (rng.next() % 1000) as i64).collect();
        let c: Vec<i64> = (0..d).map(|_| (rng.next() % 1000) as i64).collect();
        let k = rng.next() % 500;
        let dv = direct(&init, &c, k);
        let fv = fiduccia(&init, &c, k);
        let bv = kth_term_bostan_mori(&init, &c, k);
        if dv != fv || dv != bv {
            random_ok = false;
            println!("RANDOM MISMATCH d={d} k={k}: direct={dv} fiduccia={fv} bostan_mori={bv}");
        }
    }
    println!("random small-order recurrences (direct vs. both fast algorithms): {}", if random_ok { "PASS" } else { "FAIL" });

    // --- Huge K: the two fast algorithms must still agree with each other
    let huge_k = 1_000_000_000_000u64;
    let f = fiduccia(&fib_init, &fib_c, huge_k);
    let b = kth_term_bostan_mori(&fib_init, &fib_c, huge_k);
    println!("\nFibonacci(10^12) mod {MOD}: Fiduccia = {f}, Bostan-Mori = {b}, agree = {}", f == b);
}
