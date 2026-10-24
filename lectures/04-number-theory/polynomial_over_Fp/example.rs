// Polynomials over F_p: the division/remainder/factor theorems that hold
// because F_p is a field (every nonzero element is invertible), the
// bound "a degree-n polynomial has at most n roots," the resulting
// factorization x^p - x = prod_{a in F_p} (x - a) (hence Wilson's
// theorem), and the Schwartz-Zippel lemma bounding how often a nonzero
// multivariate polynomial vanishes at a random point.
//
// This file verifies all of these computationally: polynomial long
// division against a brute-force check, the root-count bound against
// naive root-finding by evaluation, x^p - x's factorization by directly
// expanding prod (x-a) and comparing coefficients, Wilson's theorem by
// direct computation, and Schwartz-Zippel by measuring vanishing
// frequency of random nonzero polynomials against the d/p bound.

// PDFBEGIN
const P: i64 = 17; // small prime, so exhaustive checks over all of F_p are cheap

fn add(a: i64, b: i64) -> i64 {
    (a + b).rem_euclid(P)
}

fn sub(a: i64, b: i64) -> i64 {
    (a - b).rem_euclid(P)
}

fn mul(a: i64, b: i64) -> i64 {
    (a.rem_euclid(P) * b.rem_euclid(P)).rem_euclid(P)
}

fn mod_pow(mut base: i64, mut exp: i64, m: i64) -> i64 {
    let mut result = 1i64;
    base = base.rem_euclid(m);
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base).rem_euclid(m);
        }
        base = (base * base).rem_euclid(m);
        exp >>= 1;
    }
    result
}

fn mod_inv(a: i64) -> i64 {
    mod_pow(a, P - 2, P)
}

/// Coefficients stored low-degree first: poly[i] is the coefficient of x^i.
type Poly = Vec<i64>;

fn trim(mut a: Poly) -> Poly {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
    a
}

fn degree(a: &Poly) -> i64 {
    let a = trim(a.clone());
    if a.len() == 1 && a[0] == 0 {
        -1 // the zero polynomial: degree -infinity, represented as -1 here
    } else {
        (a.len() - 1) as i64
    }
}

fn eval_poly(a: &Poly, x: i64) -> i64 {
    // Horner's method.
    let mut acc = 0i64;
    for &c in a.iter().rev() {
        acc = add(mul(acc, x), c);
    }
    acc
}

/// Polynomial long division: returns (quotient, remainder) with
/// deg(remainder) < deg(divisor). Requires divisor != 0.
fn poly_divmod(f: &Poly, g: &Poly) -> (Poly, Poly) {
    assert!(degree(g) >= 0, "division by the zero polynomial");
    let mut r = trim(f.clone());
    let dg = degree(g);
    let lead_inv = mod_inv(*g.last().unwrap());
    let mut q = vec![0i64; 1];
    while degree(&r) >= dg {
        let dr = degree(&r) as usize;
        let shift = dr - dg as usize;
        let coeff = mul(*r.last().unwrap(), lead_inv);
        if q.len() < shift + 1 {
            q.resize(shift + 1, 0);
        }
        q[shift] = coeff;
        for (i, &gi) in g.iter().enumerate() {
            r[shift + i] = sub(r[shift + i], mul(coeff, gi));
        }
        r = trim(r);
        if degree(&r) < 0 {
            break;
        }
    }
    (trim(q), trim(r))
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if degree(a) < 0 || degree(b) < 0 {
        return vec![0];
    }
    let mut out = vec![0i64; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] = add(out[i + j], mul(ai, bj));
        }
    }
    trim(out)
}

/// All roots of f in F_p, found by brute-force evaluation (P is small).
fn roots_by_evaluation(f: &Poly) -> Vec<i64> {
    (0..P).filter(|&a| eval_poly(f, a) == 0).collect()
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn main() {
    // --- Remainder theorem: f(x) mod (x - a) == f(a), for many random f, a ---
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
    let mut rng = Rng(20260320);
    let mut remainder_ok = true;
    for _ in 0..300 {
        let deg = rng.range(0, 6) as usize;
        let f: Poly = (0..=deg).map(|_| rng.range(0, P - 1)).collect();
        let f = trim(f);
        let a = rng.range(0, P - 1);
        let divisor = vec![sub(0, a), 1]; // x - a
        let (_, r) = poly_divmod(&f, &divisor);
        let r_const = if degree(&r) < 0 { 0 } else { r[0] };
        if r_const != eval_poly(&f, a) {
            remainder_ok = false;
            println!("REMAINDER MISMATCH f={:?} a={a}: r={r_const} f(a)={}", f, eval_poly(&f, a));
        }
    }
    println!(
        "Remainder theorem (f(x) mod (x-a) == f(a)) over {} random trials: {}",
        300,
        if remainder_ok { "PASS" } else { "FAIL" }
    );

    // --- Factor theorem + root count bound: brute-force roots vs. the
    //     factorization f(x) = c * prod (x - root_i) * cofactor, checked
    //     by re-expanding and comparing to f, and root count <= degree. ---
    let mut factor_ok = true;
    let mut root_bound_ok = true;
    for _ in 0..300 {
        let deg = rng.range(1, 6) as usize;
        let f: Poly = {
            let mut v: Poly = (0..=deg).map(|_| rng.range(0, P - 1)).collect();
            if v[deg] == 0 {
                v[deg] = 1; // ensure the intended degree is exact
            }
            v
        };
        let roots = roots_by_evaluation(&f);
        if roots.len() > degree(&f) as usize {
            root_bound_ok = false;
            println!("ROOT COUNT BOUND VIOLATED f={:?} roots={:?}", f, roots);
        }
        // Reconstruct prod (x - root) for the found roots and check it divides f exactly.
        let mut prod: Poly = vec![1];
        for &r in &roots {
            prod = poly_mul(&prod, &vec![sub(0, r), 1]);
        }
        if degree(&prod) >= 0 {
            let (q, rem) = poly_divmod(&f, &prod);
            if degree(&rem) >= 0 {
                factor_ok = false;
                println!("FACTOR THEOREM: prod(x-root) does not divide f; f={:?} roots={:?}", f, roots);
            }
            // Re-expanding q * prod should reproduce f exactly.
            let reexpanded = poly_mul(&q, &prod);
            if reexpanded != trim(f.clone()) {
                factor_ok = false;
                println!("FACTOR THEOREM: re-expansion mismatch f={:?}", f);
            }
        }
    }
    println!("Root count bound (#roots <= degree) over 300 random trials: {}", if root_bound_ok { "PASS" } else { "FAIL" });
    println!("Factor theorem (prod(x-root_i) divides f, re-expands to f) over 300 random trials: {}", if factor_ok { "PASS" } else { "FAIL" });

    // --- x^p - x = prod_{a in F_p} (x - a), checked by direct expansion ---
    let mut rhs: Poly = vec![1];
    for a in 0..P {
        rhs = poly_mul(&rhs, &vec![sub(0, a), 1]);
    }
    let mut lhs = vec![0i64; (P as usize) + 1];
    lhs[P as usize] = 1;
    lhs[1] = sub(lhs[1], 1); // x^p - x
    let xp_x_ok = trim(lhs) == trim(rhs.clone());
    println!("x^p - x == prod_{{a in F_p}}(x-a), checked by direct expansion (p={P}): {}", if xp_x_ok { "PASS" } else { "FAIL" });

    // --- Wilson's theorem: compare (p-1)! mod p to the coefficient of x^1
    //     in prod_{a in F_p}(x-a), which should match by the corollary
    //     above (both should equal -1 mod p, i.e. p-1). ---
    let mut factorial = 1i64;
    for k in 1..P {
        factorial = mul(factorial, k);
    }
    let coeff_x1 = if rhs.len() > 1 { rhs[1] } else { 0 };
    let wilson_ok = factorial == P - 1 && coeff_x1 == P - 1;
    println!("Wilson's theorem ((p-1)! == -1 mod p, matching the x^1 coefficient above): {}", if wilson_ok { "PASS" } else { "FAIL" });

    // --- Schwartz-Zippel lemma: for random nonzero low-degree multivariate
    //     polynomials (here bivariate, f(x,y) built from monomials up to
    //     total degree d), estimate P[f(a,b)=0] over random (a,b) and check
    //     it stays <= d/p (with slack for statistical noise). ---
    let mut sz_ok = true;
    for _ in 0..30 {
        let d = rng.range(1, 5) as usize; // total degree bound
        // Build a random nonzero polynomial as a sum of monomials x^i * y^j with i+j <= d.
        let mut terms: Vec<(usize, usize, i64)> = Vec::new();
        for i in 0..=d {
            for j in 0..=(d - i) {
                let c = rng.range(0, P - 1);
                if c != 0 {
                    terms.push((i, j, c));
                }
            }
        }
        if terms.is_empty() {
            continue;
        }
        let eval_f = |x: i64, y: i64| -> i64 {
            let mut acc = 0i64;
            for &(i, j, c) in &terms {
                acc = add(acc, mul(c, mul(mod_pow(x, i as i64, P), mod_pow(y, j as i64, P))));
            }
            acc
        };
        let mut zero_count = 0u64;
        let total = (P * P) as u64; // exhaustive over all of F_p x F_p, since p is small
        for x in 0..P {
            for y in 0..P {
                if eval_f(x, y) == 0 {
                    zero_count += 1;
                }
            }
        }
        let observed = zero_count as f64 / total as f64;
        let bound = d as f64 / P as f64;
        if observed > bound + 1e-9 {
            sz_ok = false;
            println!("SCHWARTZ-ZIPPEL VIOLATED d={d}: observed={observed:.4} bound={bound:.4}");
        }
    }
    println!("Schwartz-Zippel lemma (P[f=0] <= d/p, checked exhaustively over F_p x F_p): {}", if sz_ok { "PASS" } else { "FAIL" });
}
