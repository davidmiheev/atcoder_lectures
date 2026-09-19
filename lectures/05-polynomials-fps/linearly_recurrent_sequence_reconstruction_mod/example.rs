// Reeds-Sloane algorithm: recover a minimal-order linear recurrence for a
// sequence over Z/mZ, for arbitrary (not necessarily prime) modulus m.
//
// Strategy: factor m = prod p_j^{e_j}, solve the recovery problem modulo
// each prime power separately (the hard part -- ordinary Berlekamp-Massey
// does not work since Z/p^eZ is not a field), then recombine coefficient-
// wise with the Chinese Remainder Theorem.

fn mod_pow(mut base: i64, mut exp: u32, m: i64) -> i64 {
    let mut result = 1i64 % m;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % m;
        }
        base = base * base % m;
        exp >>= 1;
    }
    result
}

fn ext_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = ext_gcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}

/// Inverse of a unit `a` modulo `m` (requires gcd(a, m) = 1).
fn mod_inv(a: i64, m: i64) -> i64 {
    let (g, x, _) = ext_gcd(((a % m) + m) % m, m);
    assert_eq!(g, 1, "{} is not invertible mod {}", a, m);
    ((x % m) + m) % m
}

fn factorize(mut m: i64) -> Vec<(i64, u32)> {
    let mut factors = Vec::new();
    let mut p = 2i64;
    while p * p <= m {
        if m % p == 0 {
            let mut e = 0;
            while m % p == 0 {
                m /= p;
                e += 1;
            }
            factors.push((p, e));
        }
        p += 1;
    }
    if m > 1 {
        factors.push((m, 1));
    }
    factors
}

/// Decompose a nonzero residue x in [0, p^e) as x = t * p^u with t a unit
/// mod p^e. Returns the sentinel (e, 1) for x = 0 ("infinite valuation").
fn valuation_unit(x: i64, p: i64, e: u32, pe: i64) -> (u32, i64) {
    let mut x = ((x % pe) + pe) % pe;
    if x == 0 {
        return (e, 1);
    }
    let mut u = 0u32;
    while x % p == 0 {
        x /= p;
        u += 1;
    }
    (u, x)
}

/// Reeds-Sloane recovery of the minimal (L, Q(x)) with Q(0) = 1, over the
/// ring Z/p^eZ, for the sequence `s`. Returns (L, Q) where Q has L + 1
/// coefficients, Q[0] = 1.
fn reeds_sloane_prime_power(s: &[i64], p: i64, e: u32) -> (usize, Vec<i64>) {
    let pe = p.pow(e);
    let s: Vec<i64> = s.iter().map(|&x| ((x % pe) + pe) % pe).collect();
    let n_terms = s.len();
    let e_us = e as usize;

    // Slot j tracks (Q_j, L_j) with Q_j(0) = p^j, the minimal-order
    // recurrence at "valuation level" j seen so far.
    let mut q: Vec<Vec<i64>> = (0..e_us).map(|j| vec![mod_pow(p, j as u32, pe)]).collect();
    let mut l: Vec<usize> = vec![0; e_us];
    // Auxiliary "last improvement" data per slot: (B, L_B, n_B, t_B).
    let mut aux: Vec<Option<(Vec<i64>, usize, usize, i64)>> = vec![None; e_us];

    for n in 0..n_terms {
        // Discrepancy of each slot's current recurrence at position n.
        let mut delta = vec![0i64; e_us];
        for j in 0..e_us {
            let mut d = 0i64;
            for (i, &qi) in q[j].iter().enumerate() {
                if i > n {
                    break;
                }
                d = (d + qi * s[n - i]) % pe;
            }
            delta[j] = d;
        }
        let mut t = vec![0i64; e_us];
        let mut u = vec![0u32; e_us];
        for j in 0..e_us {
            let (uj, tj) = valuation_unit(delta[j], p, e, pe);
            u[j] = uj;
            t[j] = tj;
        }

        let mut q_new = q.clone();
        let mut l_new = l.clone();

        for j in 0..e_us {
            if u[j] as usize == e_us {
                continue; // discrepancy is exactly 0: nothing to fix
            }
            let k = e_us - 1 - u[j] as usize;
            match &aux[k] {
                None => {
                    // No donor available yet: record the failure, keep Q_j.
                    l_new[j] = n + 1;
                }
                Some((bk, lbk, nbk, tbk)) => {
                    let c = t[j] * mod_inv(*tbk, pe) % pe;
                    let shift = n - nbk;
                    l_new[j] = l[j].max(lbk + shift);
                    let mut poly = vec![0i64; l_new[j] + 1];
                    for (i, &qi) in q[j].iter().enumerate() {
                        poly[i] = qi;
                    }
                    for (i, &bi) in bk.iter().enumerate() {
                        let idx = shift + i;
                        poly[idx] = ((poly[idx] - c * bi) % pe + pe) % pe;
                    }
                    q_new[j] = poly;
                }
            }
            if q_new[j].len() != l_new[j] + 1 {
                q_new[j].resize(l_new[j] + 1, 0);
            }
        }

        // Save auxiliary data using the OLD (pre-update) donor slot k.
        let mut aux_new = aux.clone();
        for j in 0..e_us {
            if l[j] < l_new[j] {
                let k = e_us - 1 - u[j] as usize;
                aux_new[j] = Some((q[k].clone(), l[k], n, t[k]));
            }
        }

        q = q_new;
        l = l_new;
        aux = aux_new;
    }

    (l[0], q[0].clone())
}

fn crt_pair(r1: i64, m1: i64, r2: i64, m2: i64) -> (i64, i64) {
    let (g, x, _) = ext_gcd(m1, m2);
    debug_assert_eq!(g, 1, "moduli must be coprime");
    let m = m1 * m2;
    let diff = ((r2 - r1) % m2 + m2) % m2;
    let t = (diff * x) % m2;
    let r = ((r1 + m1 * t) % m + m) % m;
    (r, m)
}

/// Full Reeds-Sloane recovery for an arbitrary modulus m: solve modulo each
/// prime power factor, then recombine every coefficient via CRT.
fn reeds_sloane(s: &[i64], m: i64) -> (usize, Vec<i64>) {
    let factors = factorize(m);
    let per_factor: Vec<(i64, usize, Vec<i64>)> = factors
        .iter()
        .map(|&(p, e)| {
            let (l, q) = reeds_sloane_prime_power(s, p, e);
            (p.pow(e), l, q)
        })
        .collect();

    let l = per_factor.iter().map(|&(_, l, _)| l).max().unwrap();

    let mut q = vec![0i64; l + 1];
    let mut modulus = 1i64;
    q[0] = 0;
    for i in 0..=l {
        let (mut r, mut cur_mod) = (0i64, 1i64);
        for (pe, lj, qj) in &per_factor {
            let coeff = if i <= *lj { qj[i] } else { 0 };
            let (nr, nmod) = crt_pair(r, cur_mod, coeff, *pe);
            r = nr;
            cur_mod = nmod;
        }
        q[i] = r;
        modulus = cur_mod;
    }
    debug_assert_eq!(modulus, m);
    (l, q)
}

/// Extend a sequence using S_i = -sum_{j=1}^{L} Q[j] * S[i-j] (mod m), given
/// Q[0] = 1. Returns `count` further terms after the known prefix.
fn extrapolate(known: &[i64], q: &[i64], m: i64, count: usize) -> Vec<i64> {
    let mut s = known.to_vec();
    let l = q.len() - 1;
    for _ in 0..count {
        let n = s.len();
        let mut acc = 0i64;
        for j in 1..=l {
            if j <= n {
                acc = (acc + q[j] * s[n - j]) % m;
            }
        }
        s.push(((-acc) % m + m) % m);
    }
    s[known.len()..].to_vec()
}

fn main() {
    // --- Example 1: composite, non-prime-power modulus ---------------
    // Fibonacci-like sequence mod 12 = 4 * 3, satisfying S_i = S_{i-1} + S_{i-2}.
    let m = 12i64;
    let mut fib = vec![1i64, 1];
    for i in 2..20 {
        let next = (fib[i - 1] + fib[i - 2]) % m;
        fib.push(next);
    }
    let known = &fib[..10];
    let (l, q) = reeds_sloane(known, m);
    println!("mod 12 Fibonacci-like: recovered L = {}, Q = {:?}", l, q);
    let predicted = extrapolate(known, &q, m, fib.len() - known.len());
    println!("predicted tail: {:?}", predicted);
    println!("actual tail:    {:?}", &fib[known.len()..]);
    println!("match: {}\n", predicted == fib[known.len()..]);

    // --- Example 2: the classic non-uniqueness witness ----------------
    // S = (2,2,2,...) mod 4 satisfies both S_i = S_{i-1} and S_i = 3 S_{i-1}.
    let s2 = vec![2i64; 8];
    let (l2, q2) = reeds_sloane(&s2, 4);
    println!("mod 4 constant sequence: recovered L = {}, Q = {:?}", l2, q2);
    let more = extrapolate(&s2, &q2, 4, 4);
    println!("predicted continuation: {:?} (all should be 2)\n", more);

    // --- Example 3: prime modulus, sanity check against a field --------
    // Over Z/5Z this must reduce to ordinary Berlekamp-Massey behaviour.
    let m3 = 5i64;
    let mut trib = vec![1i64, 1, 2];
    for i in 3..12 {
        let next = ((trib[i - 1] + trib[i - 2] + trib[i - 3]) % m3 + m3) % m3;
        trib.push(next);
    }
    let known3 = &trib[..8];
    let (l3, q3) = reeds_sloane(known3, m3);
    println!("mod 5 tribonacci-like: recovered L = {}, Q = {:?}", l3, q3);
    let pred3 = extrapolate(known3, &q3, m3, trib.len() - known3.len());
    println!("match: {}", pred3 == trib[known3.len()..]);
}
