// Zeta and Mobius transforms with respect to the subset/superset relation
// on bitmasks. The subset zeta transform B_T = sum_{S subseteq T} A_S
// (a.k.a. "sum over subsets", SOS) and its inverse, the subset Mobius
// transform, are computed here three ways -- an O(4^n) definition-level
// double loop, an O(3^n) submask-enumeration loop, and the O(n*2^n) fast
// zeta/Mobius transform that flips one bit's role (fixed value -> wildcard)
// at a time -- and cross-checked against each other. The superset versions
// are analogous with the inclusion direction reversed. Finally, OR- and
// AND-convolution are built on top: C_U = sum_{S|T=U} A_S*B_T reduces to a
// pointwise product after a subset zeta transform (since S|T=U implies, and
// nearly is implied by, S,T subseteq U) followed by a subset Mobius
// transform back; AND-convolution is the superset-relation mirror image.
// All of these are cross-checked against brute-force O(4^n) definitions.

/// O(4^n): definition-level double loop.
fn subset_zeta_naive(n: usize, a: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut b = vec![0i64; size];
    for t in 0..size {
        for s in 0..size {
            if s & t == s {
                b[t] += a[s];
            }
        }
    }
    b
}

/// O(3^n): enumerate submasks of T directly via the standard
/// `for (s = t; ; s = (s-1)&t)` trick.
fn subset_zeta_submask_enum(n: usize, a: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut b = vec![0i64; size];
    for t in 0..size {
        let mut s = t;
        loop {
            b[t] += a[s];
            if s == 0 {
                break;
            }
            s = (s.wrapping_sub(1)) & t;
        }
    }
    b
}

/// O(n*2^n): fast zeta transform (FZT). In place.
fn subset_zeta_fast(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 1 {
                a[t] += a[t ^ (1 << i)];
            }
        }
    }
    a
}

/// O(n*2^n): fast Mobius transform (FMT), the exact inverse of the fast
/// zeta transform above.
fn subset_mobius_fast(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 1 {
                a[t] -= a[t ^ (1 << i)];
            }
        }
    }
    a
}

/// O(4^n): definition-level double loop for the superset zeta transform,
/// B_T = sum_{S superseteq T} A_S.
fn superset_zeta_naive(n: usize, a: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut b = vec![0i64; size];
    for t in 0..size {
        for s in 0..size {
            if s & t == t {
                b[t] += a[s];
            }
        }
    }
    b
}

/// O(n*2^n): fast superset zeta transform (mirror image of the subset one:
/// flip the roles of 0 and 1 in which bits get "widened" to wildcards).
fn superset_zeta_fast(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 0 {
                a[t] += a[t ^ (1 << i)];
            }
        }
    }
    a
}

/// O(n*2^n): fast superset Mobius transform, the inverse of the above.
fn superset_mobius_fast(n: usize, mut a: Vec<i64>) -> Vec<i64> {
    let size = 1usize << n;
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 0 {
                a[t] -= a[t ^ (1 << i)];
            }
        }
    }
    a
}

/// OR-convolution: C_U = sum_{S|T=U} A_S * B_T, via subset zeta -> pointwise
/// product -> subset Mobius.
fn or_convolution(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let za = subset_zeta_fast(n, a.to_vec());
    let zb = subset_zeta_fast(n, b.to_vec());
    let d: Vec<i64> = za.iter().zip(zb.iter()).map(|(&x, &y)| x * y).collect();
    subset_mobius_fast(n, d)
}

/// AND-convolution: C_U = sum_{S&T=U} A_S * B_T, via superset zeta ->
/// pointwise product -> superset Mobius.
fn and_convolution(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let za = superset_zeta_fast(n, a.to_vec());
    let zb = superset_zeta_fast(n, b.to_vec());
    let d: Vec<i64> = za.iter().zip(zb.iter()).map(|(&x, &y)| x * y).collect();
    superset_mobius_fast(n, d)
}

/// Brute-force O(4^n) reference for OR-convolution.
fn or_convolution_naive(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut c = vec![0i64; size];
    for s in 0..size {
        for t in 0..size {
            c[s | t] += a[s] * b[t];
        }
    }
    c
}

/// Brute-force O(4^n) reference for AND-convolution.
fn and_convolution_naive(n: usize, a: &[i64], b: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    let mut c = vec![0i64; size];
    for s in 0..size {
        for t in 0..size {
            c[s & t] += a[s] * b[t];
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
    let mut rng = Rng(2026091801);

    // --- Worked example from the lecture: A = (3,1,4,1,5,9,2,6), n=3 ---
    let a_example = vec![3i64, 1, 4, 1, 5, 9, 2, 6];
    let expected_b = vec![3i64, 4, 7, 9, 8, 18, 14, 31];
    let b_via_naive = subset_zeta_naive(3, &a_example);
    let b_via_submask = subset_zeta_submask_enum(3, &a_example);
    let b_via_fast = subset_zeta_fast(3, a_example.clone());
    let example_ok = b_via_naive == expected_b && b_via_submask == expected_b && b_via_fast == expected_b;
    println!("Worked example (n=3): all three zeta transform implementations match the lecture's own B: {}", if example_ok { "PASS" } else { "FAIL" });

    // --- Round-trip: subset Mobius(subset zeta(A)) == A ---
    let mut roundtrip_ok = true;
    let mut cross_impl_ok = true;
    for _ in 0..100 {
        let n = rng.range(0, 8) as usize;
        let size = 1usize << n;
        let a: Vec<i64> = (0..size).map(|_| rng.range(-50, 50)).collect();

        let b_naive = subset_zeta_naive(n, &a);
        let b_submask = subset_zeta_submask_enum(n, &a);
        let b_fast = subset_zeta_fast(n, a.clone());
        if b_naive != b_submask || b_naive != b_fast {
            cross_impl_ok = false;
            println!("SUBSET ZETA CROSS-IMPL MISMATCH n={n}");
        }

        let a_back = subset_mobius_fast(n, b_fast.clone());
        if a_back != a {
            roundtrip_ok = false;
            println!("SUBSET ROUND-TRIP MISMATCH n={n}");
        }

        // Superset transforms: cross-check fast against naive, and round-trip.
        let sup_naive = superset_zeta_naive(n, &a);
        let sup_fast = superset_zeta_fast(n, a.clone());
        if sup_naive != sup_fast {
            cross_impl_ok = false;
            println!("SUPERSET ZETA CROSS-IMPL MISMATCH n={n}");
        }
        let sup_back = superset_mobius_fast(n, sup_fast.clone());
        if sup_back != a {
            roundtrip_ok = false;
            println!("SUPERSET ROUND-TRIP MISMATCH n={n}");
        }
    }
    println!("Subset/superset zeta transform: naive (O(4^n)), submask-enum (O(3^n)), and fast (O(n*2^n)) all agree, over 100 random trials: {}", if cross_impl_ok { "PASS" } else { "FAIL" });
    println!("Mobius transform exactly inverts the zeta transform (subset and superset), over 100 random trials: {}", if roundtrip_ok { "PASS" } else { "FAIL" });

    // --- OR- and AND-convolution vs. brute force ---
    let mut or_ok = true;
    let mut and_ok = true;
    for _ in 0..100 {
        let n = rng.range(0, 7) as usize;
        let size = 1usize << n;
        let a: Vec<i64> = (0..size).map(|_| rng.range(-10, 10)).collect();
        let b: Vec<i64> = (0..size).map(|_| rng.range(-10, 10)).collect();

        let or_fast = or_convolution(n, &a, &b);
        let or_naive = or_convolution_naive(n, &a, &b);
        if or_fast != or_naive {
            or_ok = false;
            println!("OR-CONVOLUTION MISMATCH n={n}");
        }

        let and_fast = and_convolution(n, &a, &b);
        let and_naive = and_convolution_naive(n, &a, &b);
        if and_fast != and_naive {
            and_ok = false;
            println!("AND-CONVOLUTION MISMATCH n={n}");
        }
    }
    println!("OR-convolution (zeta -> pointwise product -> Mobius) matches brute force, over 100 random trials: {}", if or_ok { "PASS" } else { "FAIL" });
    println!("AND-convolution (superset zeta -> pointwise product -> superset Mobius) matches brute force, over 100 random trials: {}", if and_ok { "PASS" } else { "FAIL" });

    let all_ok = example_ok && cross_impl_ok && roundtrip_ok && or_ok && and_ok;
    println!("\nAll checks: {}", if all_ok { "PASS" } else { "FAIL" });
}
