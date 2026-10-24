// DP acceleration via eigenfunctions: some polynomial-valued DP
// recurrences A_{i+1} = f(A_i) are linear maps f on the space of
// polynomials, and if f happens to have a known eigenbasis -- functions
// v with f(v) = lambda*v -- then M iterations collapse to "expand A_0 in
// that basis, multiply each coefficient by lambda^M, recombine," turning
// an O(M) loop into an O(log M) computation (plus the cost of the
// basis-change itself).
//
// This file fully implements and tests Problem 2 from the lecture, whose
// eigenbasis is the elementary (1+x)^k (found via Taylor shift, no
// composition machinery needed). For Problems 3 and 4, whose eigenbases
// involve e^{kx} and x^k e^{-x}, it numerically verifies the eigenvalue
// identities themselves (via a from-scratch power-series exp), including
// the truncation subtlety the lecture calls out: e^{kx} mod x^{N+1} is
// only an APPROXIMATE eigenfunction once you truncate too early, which is
// exactly why those algorithms pad the working degree to N+M before
// truncating back down.

// PDFBEGIN
const MOD: i64 = 998244353;

fn add(a: i64, b: i64) -> i64 {
    (a + b).rem_euclid(MOD)
}

fn mul(a: i64, b: i64) -> i64 {
    (a.rem_euclid(MOD) * b.rem_euclid(MOD)).rem_euclid(MOD)
}

fn mod_pow(mut base: i64, mut exp: i64) -> i64 {
    let mut result = 1i64;
    base = base.rem_euclid(MOD);
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul(result, base);
        }
        base = mul(base, base);
        exp >>= 1;
    }
    result
}

fn mod_inv(a: i64) -> i64 {
    mod_pow(a, MOD - 2)
}

// ---------------------------------------------------------------------
// Problem 2: A_{i+1} = (1+x) A_i'.  Eigenfunctions (1+x)^k, eigenvalue k.
// ---------------------------------------------------------------------

/// Pascal's triangle mod p, up to row `n`.
fn binomial_table(n: usize) -> Vec<Vec<i64>> {
    let mut c = vec![vec![0i64; n + 1]; n + 1];
    for i in 0..=n {
        c[i][0] = 1;
        for j in 1..=i {
            c[i][j] = add(c[i - 1][j - 1], if j <= i - 1 { c[i - 1][j] } else { 0 });
        }
    }
    c
}

/// Taylor shift: given P of degree <= n, compute Q(x) = P(x + t), O(n^2).
fn taylor_shift(a: &[i64], t: i64) -> Vec<i64> {
    let n = a.len() - 1;
    let c = binomial_table(n);
    let mut t_pow = vec![1i64; n + 1];
    for i in 1..=n {
        t_pow[i] = mul(t_pow[i - 1], t);
    }
    let mut q = vec![0i64; n + 1];
    for i in 0..=n {
        let mut s = 0i64;
        for j in i..=n {
            s = add(s, mul(a[j], mul(c[j][i], t_pow[j - i])));
        }
        q[i] = s;
    }
    q
}

/// Solve Problem 2 via the eigenfunction method: O(N^2) here (Taylor
/// shift implemented naively); O(N log N) with a fast Taylor shift.
fn solve_problem2_eigenfunction(a0: &[i64], m: u64) -> Vec<i64> {
    // Expand A_0(x) = sum_k c_k (1+x)^k: substitute x -> x-1, i.e. shift by -1.
    let mut c = taylor_shift(a0, MOD - 1);
    for (k, ck) in c.iter_mut().enumerate() {
        *ck = mul(*ck, mod_pow(k as i64, m as i64));
    }
    // Recombine: A_M(x) = sum_k c_k (1+x)^k = D(x+1), shift by +1.
    taylor_shift(&c, 1)
}

/// The definition, directly: simulate the DP table. O(N*M). Used only to
/// cross-check the eigenfunction-based solution above.
fn solve_problem2_direct(a0: &[i64], m: u64) -> Vec<i64> {
    let n = a0.len() - 1;
    let mut dp = a0.to_vec();
    for _ in 0..m {
        let mut next = vec![0i64; n + 1];
        for j in 0..=n {
            let dp_j1 = if j + 1 <= n { dp[j + 1] } else { 0 };
            next[j] = add(mul(j as i64, dp[j]), mul((j + 1) as i64, dp_j1));
        }
        dp = next;
    }
    dp
}

// ---------------------------------------------------------------------
// Power series exp, from scratch (O(n^2)), used to numerically verify
// the eigenfunction claims in Problems 3 and 4 without building a full
// composition/reversion engine.
// ---------------------------------------------------------------------

/// exp(f) mod x^n, given f[0] = 0. Standard E' = f' E recurrence:
/// k*E_k = sum_{i=1}^{k} i*f_i*E_{k-i}.
fn power_series_exp(f: &[i64], n: usize) -> Vec<i64> {
    assert_eq!(f.get(0).copied().unwrap_or(0), 0);
    let mut e = vec![0i64; n];
    e[0] = 1;
    for k in 1..n {
        let mut s = 0i64;
        for i in 1..=k.min(f.len() - 1) {
            s = add(s, mul(mul(i as i64, f[i]), e[k - i]));
        }
        e[k] = mul(s, mod_inv(k as i64));
    }
    e
}

fn derivative(a: &[i64]) -> Vec<i64> {
    (1..a.len()).map(|i| mul(i as i64, a[i])).collect()
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    (0..n).map(|i| add(a.get(i).copied().unwrap_or(0), b.get(i).copied().unwrap_or(0))).collect()
}

fn poly_scale(a: &[i64], s: i64) -> Vec<i64> {
    a.iter().map(|&x| mul(x, s)).collect()
}

/// x^k * e^{-x}, as a power series truncated to `len` terms.
fn x_pow_k_times_exp_neg_x(k: usize, len: usize) -> Vec<i64> {
    let mut f = vec![0i64; len.max(1)];
    if len > 1 {
        f[1] = MOD - 1; // f(x) = -x
    }
    let e = power_series_exp(&f, len.saturating_sub(k).max(1));
    let mut out = vec![0i64; len];
    for (i, &v) in e.iter().enumerate() {
        if k + i < len {
            out[k + i] = v;
        }
    }
    out
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn main() {
    // --- Problem 2: eigenfunction method vs. direct DP simulation --------
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
    let mut rng = Rng(31415);
    let mut p2_ok = true;
    for _ in 0..200 {
        let n = 1 + (rng.next() % 8) as usize;
        let m = rng.range(0, 20) as u64;
        let a0: Vec<i64> = (0..=n).map(|_| rng.range(0, MOD - 1)).collect();
        let via_eigen = solve_problem2_eigenfunction(&a0, m);
        let via_direct = solve_problem2_direct(&a0, m);
        if via_eigen != via_direct {
            p2_ok = false;
            println!("PROBLEM2 MISMATCH a0={:?} m={m}: eigen={:?} direct={:?}", a0, via_eigen, via_direct);
        }
    }
    println!("Problem 2 (A_i+1=(1+x)A_i'): eigenfunction method matches direct DP simulation: {}", if p2_ok { "PASS" } else { "FAIL" });

    // --- Problem 3's eigenfunction claim: e^{kx} has eigenvalue (k+1) ----
    // under f(A) = A + A', PROVIDED enough terms are kept -- verify both
    // that it holds with enough padding, and that it visibly breaks at
    // the top degree without it (the lecture's own stated caveat).
    let mut p3_exact_ok = true;
    let mut p3_caveat_demonstrated = false;
    for trial in 0..50 {
        let k = (trial % 7) as i64;
        let n = 6usize; // "requested" precision: coefficients 0..=n
        let pad = 3; // extra terms kept beyond n, to be safe
        let f: Vec<i64> = { let mut v = vec![0i64; n + pad + 2]; if v.len() > 1 { v[1] = k; } v };
        let a_full = power_series_exp(&f, n + pad + 1); // e^{kx} mod x^{n+pad+1}
        let a_deriv = derivative(&a_full); // exact for indices 0..=n+pad-1
        let lhs = poly_add(&a_full[..n + pad], &a_deriv[..n + pad]);
        let rhs = poly_scale(&a_full[..n + pad], k + 1);
        if lhs != rhs {
            p3_exact_ok = false;
            println!("PROBLEM3 EIGEN MISMATCH k={k}: lhs={:?} rhs={:?}", lhs, rhs);
        }

        // Now the under-padded version: truncate A to exactly n+1 terms,
        // compute A' from THAT truncation (losing the true top term),
        // and check the identity fails at degree n (but still holds for
        // degrees 0..n-1).
        let a_short = &a_full[..=n];
        let a_short_deriv_true_len = derivative(a_short); // length n, indices 0..n-1
        let lhs_short = poly_add(a_short, &a_short_deriv_true_len);
        let rhs_short = poly_scale(a_short, k + 1);
        let matches_low = lhs_short[..n] == rhs_short[..n];
        let matches_top = lhs_short[n] == rhs_short[n];
        if matches_low && !matches_top {
            p3_caveat_demonstrated = true;
        } else if !matches_low {
            println!("PROBLEM3 CAVEAT: unexpected low-degree mismatch k={k}");
        }
    }
    println!("Problem 3 eigenfunction identity A+A'=(k+1)A holds exactly with enough padding: {}", if p3_exact_ok { "PASS" } else { "FAIL" });
    println!("...and (as the lecture warns) breaks specifically at the top degree without padding: {}", if p3_caveat_demonstrated { "CONFIRMED" } else { "NOT OBSERVED" });

    // --- Problem 4's eigenfunction claim: x^k e^{-x} has eigenvalue k ----
    // under f(A) = x(A + A').
    let mut p4_ok = true;
    for trial in 0..50 {
        let k = (trial % 5) as usize;
        let n = 6usize;
        let pad = 3;
        let len = k + n + pad + 2;
        let a = x_pow_k_times_exp_neg_x(k, len);
        let a_deriv = derivative(&a); // length len-1, exact throughout
        // f(A) = x*(A + A'): compute (A+A') then shift up by one degree (multiply by x).
        let sum = poly_add(&a[..len - 1], &a_deriv);
        let mut f_a = vec![0i64; len];
        for (i, &v) in sum.iter().enumerate() {
            if i + 1 < len {
                f_a[i + 1] = v;
            }
        }
        let want = poly_scale(&a, k as i64);
        // Compare on a safely-interior range (well clear of both ends).
        let check_len = k + n;
        if f_a[..check_len] != want[..check_len] {
            p4_ok = false;
            println!("PROBLEM4 EIGEN MISMATCH k={k}: f(A)={:?} k*A={:?}", &f_a[..check_len], &want[..check_len]);
        }
    }
    println!("Problem 4 eigenfunction identity x(A+A')=kA holds (for x^k e^-x) with enough padding: {}", if p4_ok { "PASS" } else { "FAIL" });
}
