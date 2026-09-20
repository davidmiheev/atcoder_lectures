// Bostan-Mori, advanced: extracting a whole INTERVAL of coefficients
// [x^[L,R]] P(x)/Q(x) in O(M(d) log R), not just a single one. This is the
// "MSB-first" variant of Bostan-Mori: instead of halving the target index
// from the bottom up (LSB-first, the ordinary Kth-term algorithm), it
// halves the interval's right end from the top down. Applications:
// shifting a linearly recurrent sequence's generating function, and
// computing x^K mod Gamma(x) (the object at the heart of Fiduccia's
// algorithm) via this machinery instead of binary exponentiation.

const MOD: i64 = 998244353;

// rem_euclid (not plain %) throughout: inputs here are sometimes raw,
// unreduced values (e.g. straight from test-data generation), and plain
// `%` keeps the dividend's sign in Rust, silently leaving negative
// "coefficients" that are only correct as a bit pattern, not as a value
// in [0, MOD).
fn add(a: i64, b: i64) -> i64 {
    (a + b).rem_euclid(MOD)
}
fn sub(a: i64, b: i64) -> i64 {
    (a - b).rem_euclid(MOD)
}
fn mul(a: i64, b: i64) -> i64 {
    (a.rem_euclid(MOD) * b.rem_euclid(MOD)).rem_euclid(MOD)
}

fn poly_mul(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        if av == 0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            out[i + j] = add(out[i + j], mul(av, bv));
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

/// ceil(a / b) for b > 0, correct for negative a.
fn ceil_div(a: i64, b: i64) -> i64 {
    -((-a).div_euclid(b))
}

/// Every interval-extraction result is returned as a PLAIN polynomial,
/// reindexed to start at x^0: index `j` holds the true coefficient at
/// degree `l + j` for the `[l, r]` interval that produced it. Once
/// returned, the array is treated purely by its own 0-based index in any
/// further multiplication -- exactly as the source's `[x^[L,R]] f(x)`
/// notation intends (it produces a fresh polynomial, not a "shifted
/// view" that remembers where it came from). Mixing up these two
/// readings is the classic off-by-shift bug in this algorithm.

/// [x^[L,R]] 1/Q(x), for Q of degree <= d with Q(0) = 1, R >= 0 (L may be
/// negative -- degrees below 0 are just 0, since 1/Q(x) is an ordinary
/// power series). `q` must have exactly d+1 entries (zero-padded).
///
/// This is the MSB-first Bostan-Mori recursion: at each step, halve R
/// (from the top) using the same Q(x)Q(-x) = V(x^2) identity as the
/// ordinary (LSB-first) algorithm, but keep track of a whole interval
/// instead of a single target index.
fn range_coefficients_of_inverse(q: &[i64], d: usize, l: i64, r: i64) -> Vec<i64> {
    if r == 0 {
        // [x^[L,0]] 1/Q(x): only degree 0 can be nonzero (it's 1, since
        // Q(0)=1), and everything from L up to -1 is a negative-degree
        // coefficient, hence 0. So the answer is the single monomial
        // x^(-L) in the OUTPUT index, i.e. a 1 at output position -L.
        let len = (0 - l + 1) as usize;
        let mut out = vec![0i64; len];
        out[(-l) as usize] = 1;
        return out;
    }

    let l_prime = ceil_div(l - d as i64, 2);
    let r_prime = r.div_euclid(2);

    let qneg = negate_odd(q);
    let v = even_part(&poly_mul(q, &qneg)); // degree <= d, since q has d+1 entries

    // W(x) = [x^[L',R']] 1/V(x), a plain polynomial: W[t] = [x^(L'+t)] 1/V(x).
    let w = range_coefficients_of_inverse(&v, d, l_prime, r_prime);

    // S(x) = [x^[L-d,R]] 1/V(x^2), built from W: since 1/V(x^2) is even,
    // [x^(L-d+j)] 1/V(x^2) is nonzero only when L-d+j = 2*(L'+t) for some
    // t, and then equals W[t]. shift = 2L'-L+d is that j at t=0.
    let s_len = (r - (l - d as i64) + 1) as usize;
    let mut s = vec![0i64; s_len];
    let shift = 2 * l_prime - l + d as i64;
    for (t, &wt) in w.iter().enumerate() {
        let j = shift + 2 * t as i64;
        if j >= 0 && (j as usize) < s_len {
            s[j as usize] = wt;
        }
    }

    // Q(-x) * S(x), as plain polynomials; then take LITERAL indices
    // [d, R-L+d] of that product (not a degree-shifted read).
    let prod = poly_mul(&qneg, &s);
    let mut out = vec![0i64; (r - l + 1) as usize];
    for (p, out_p) in out.iter_mut().enumerate() {
        let idx = d as i64 + p as i64;
        if idx >= 0 && (idx as usize) < prod.len() {
            *out_p = prod[idx as usize];
        }
    }
    out
}

/// [x^[L,R]] P(x)/Q(x), for P of degree <= d-1, Q of degree <= d,
/// Q(0) = 1. (Problem 1, reduced to Problem 2 above.)
fn range_coefficients_of_rational(p: &[i64], q: &[i64], d: usize, l: i64, r: i64) -> Vec<i64> {
    // U(x) = [x^[L-d+1,R]] 1/Q(x): U[j] = [x^(L-d+1+j)] 1/Q(x).
    let u = range_coefficients_of_inverse(q, d, l - d as i64 + 1, r);
    // P(x)*U(x): index m <-> true degree (L-d+1+m). Want true degree
    // L+p, i.e. m = p+d-1.
    let prod = poly_mul(p, &u);
    (0..=(r - l) as usize)
        .map(|p_idx| {
            let idx = p_idx as i64 + d as i64 - 1;
            if idx >= 0 && (idx as usize) < prod.len() {
                prod[idx as usize]
            } else {
                0
            }
        })
        .collect()
}

// --- Reference implementations, for cross-checking ----------------------

/// O(d^2 * N): the first `up_to`+1 coefficients of 1/Q(x) via the
/// standard recurrence, given Q(0) = 1.
fn brute_force_inverse(q: &[i64], up_to: usize) -> Vec<i64> {
    let mut inv = vec![0i64; up_to + 1];
    inv[0] = 1;
    for k in 1..=up_to {
        let mut s = 0i64;
        for i in 1..=k.min(q.len() - 1) {
            s = add(s, mul(q[i], inv[k - i]));
        }
        inv[k] = sub(0, s);
    }
    inv
}

fn brute_force_range_of_inverse(q: &[i64], l: i64, r: i64) -> Vec<i64> {
    let inv = brute_force_inverse(q, r.max(0) as usize);
    (l..=r).map(|deg| if deg < 0 { 0 } else { inv[deg as usize] }).collect()
}

/// Direct recurrence unrolling of A_i = c_1 A_{i-1} + ... + c_d A_{i-d},
/// from A(x) = P(x)/Q(x) with Q(x) = 1 - c_1 x - ... - c_d x^d.
fn direct_sequence(p: &[i64], q: &[i64], up_to: usize) -> Vec<i64> {
    let d = q.len() - 1;
    let c: Vec<i64> = (1..=d).map(|j| sub(0, q[j])).collect();
    let mut a = vec![0i64; up_to + 1];
    for i in 0..=up_to {
        let p_i = if i < p.len() { p[i] } else { 0 };
        if i < d {
            // A(x)Q(x) = P(x), Q(0)=1: A_i = P_i - sum_{j<i} A_j Q_{i-j}.
            let mut s = p_i.rem_euclid(MOD);
            for j in 0..i {
                s = sub(s, mul(a[j], q[i - j]));
            }
            a[i] = s;
        } else {
            // The recurrence proper: A_i = c_1 A_{i-1} + ... + c_d A_{i-d}.
            let mut v = 0i64;
            for j in 0..d {
                v = add(v, mul(c[j], a[i - 1 - j]));
            }
            a[i] = v;
        }
    }
    a
}

/// x^K mod Gamma(x) via ordinary binary exponentiation (Fiduccia), used
/// only to cross-check the Bostan-Mori-based method below.
fn xk_mod_gamma_fiduccia(k: u64, d: usize, c: &[i64]) -> Vec<i64> {
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
    fn mulmod_gamma(a: &[i64], b: &[i64], c: &[i64]) -> Vec<i64> {
        let d = a.len();
        let mut raw = vec![0i64; 2 * d - 1];
        for i in 0..d {
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
    fn rec(k: u64, d: usize, c: &[i64]) -> Vec<i64> {
        if (k as usize) < d {
            let mut r = vec![0i64; d];
            r[k as usize] = 1;
            return r;
        }
        let half = rec(k / 2, d, c);
        let mut r = mulmod_gamma(&half, &half, c);
        if k % 2 == 1 {
            r = mul_x_mod_gamma(&r, c);
        }
        r
    }
    rec(k, d, c)
}

/// x^K mod Gamma(x) via interval coefficient extraction (Problem 5): U(x)
/// = [x^[K-d+1,K]] 1/Q(x), then reverse (Q(x)U(x) mod x^d).
fn xk_mod_gamma_via_bostan_mori(q: &[i64], d: usize, k: i64) -> Vec<i64> {
    let u = range_coefficients_of_inverse(q, d, k - d as i64 + 1, k);
    let qu = poly_mul(q, &u);
    let r_tilde: Vec<i64> = qu[..d].to_vec();
    r_tilde.into_iter().rev().collect()
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
    let mut rng = Rng(2025);

    // --- Core: [x^[L,R]] 1/Q(x), including negative L, vs. brute force --
    let mut inv_ok = true;
    for _ in 0..300 {
        let d = 1 + (rng.next() % 5) as usize;
        let mut q = vec![0i64; d + 1];
        q[0] = 1;
        for j in 1..=d {
            q[j] = rng.range(-5, 5);
        }
        let r = rng.range(0, 20);
        let l = rng.range(-10, r);
        let got = range_coefficients_of_inverse(&q, d, l, r);
        let want = brute_force_range_of_inverse(&q, l, r);
        if got != want {
            inv_ok = false;
            println!("INVERSE MISMATCH q={:?} l={l} r={r}: got={:?} want={:?}", q, got, want);
        }
    }
    println!("[x^[L,R]] 1/Q(x) matches brute-force power series inversion (L possibly negative): {}", if inv_ok { "PASS" } else { "FAIL" });

    // --- P(x)/Q(x) interval extraction vs. direct recurrence unrolling --
    let mut rational_ok = true;
    for _ in 0..300 {
        let d = 1 + (rng.next() % 5) as usize;
        let mut q = vec![0i64; d + 1];
        q[0] = 1;
        for j in 1..=d {
            q[j] = rng.range(-5, 5);
        }
        let p: Vec<i64> = (0..d).map(|_| rng.range(-5, 5)).collect();
        let r = rng.range(0, 25);
        let l = rng.range(0, r); // keep L >= 0 to compare against the sequence directly
        let got = range_coefficients_of_rational(&p, &q, d, l, r);
        let seq = direct_sequence(&p, &q, r as usize);
        let want: Vec<i64> = (l..=r).map(|i| seq[i as usize]).collect();
        if got != want {
            rational_ok = false;
            println!("RATIONAL MISMATCH p={:?} q={:?} l={l} r={r}: got={:?} want={:?}", p, q, got, want);
        }
    }
    println!("[x^[L,R]] P(x)/Q(x) matches direct recurrence unrolling: {}", if rational_ok { "PASS" } else { "FAIL" });

    // --- Problem 5: x^K mod Gamma(x), Bostan-Mori vs. Fiduccia -----------
    let mut modgamma_ok = true;
    for _ in 0..200 {
        let d = 1 + (rng.next() % 5) as usize;
        let mut q = vec![0i64; d + 1];
        q[0] = 1;
        for j in 1..=d {
            q[j] = rng.range(-5, 5);
        }
        let c: Vec<i64> = (1..=d).map(|j| sub(0, q[j])).collect();
        let k = rng.range(0, 2000) as u64;
        let via_fiduccia = xk_mod_gamma_fiduccia(k, d, &c);
        let via_bostan_mori = xk_mod_gamma_via_bostan_mori(&q, d, k as i64);
        if via_fiduccia != via_bostan_mori {
            modgamma_ok = false;
            println!("MODGAMMA MISMATCH q={:?} k={k}: fiduccia={:?} bostan_mori={:?}", q, via_fiduccia, via_bostan_mori);
        }
    }
    println!("x^K mod Gamma(x): interval-extraction method matches Fiduccia's binary exponentiation: {}", if modgamma_ok { "PASS" } else { "FAIL" });

    // --- Problem 4: shifting a linearly recurrent sequence ---------------
    let mut shift_ok = true;
    for _ in 0..150 {
        let d = 1 + (rng.next() % 4) as usize;
        let mut q = vec![0i64; d + 1];
        q[0] = 1;
        for j in 1..=d {
            q[j] = rng.range(-4, 4);
        }
        let p: Vec<i64> = (0..d).map(|_| rng.range(-4, 4)).collect();
        let k = rng.range(0, 100) as usize;

        // T = [x^[K,K+d-1]] P(x)/Q(x); S = (T * Q) truncated to degree d-1.
        let t = range_coefficients_of_rational(&p, &q, d, k as i64, (k + d - 1) as i64);
        let mut s = poly_mul(&t, &q);
        s.truncate(d);

        // S(x)/Q(x) should generate (A_K, A_{K+1}, ...).
        let shifted_seq = direct_sequence(&s, &q, 2 * d);
        let original_seq = direct_sequence(&p, &q, k + 2 * d);
        let want: Vec<i64> = (0..=2 * d).map(|i| original_seq[k + i]).collect();
        if shifted_seq != want {
            shift_ok = false;
            println!("SHIFT MISMATCH p={:?} q={:?} k={k}: shifted={:?} want={:?}", p, q, shifted_seq, want);
        }
    }
    println!("shifted generating function S(x)/Q(x) reproduces (A_K, A_K+1, ...): {}", if shift_ok { "PASS" } else { "FAIL" });
}
