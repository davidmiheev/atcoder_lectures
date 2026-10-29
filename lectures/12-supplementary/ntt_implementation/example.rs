// The Number Theoretic Transform (NTT): a Discrete Fourier Transform
// carried out in Z/pZ instead of C, for a prime p with p-1 divisible by a
// large power of two. It gives the same O(n log n) fast-convolution trick
// as the classical FFT, but with exact integer arithmetic -- no rounding
// error, which is why it's the default for polynomial multiplication
// modulo a prime in competitive programming.
//
// This lecture is a supplement, not part of AtCoder's Algorithm Lectures
// index: several lectures in this repository (the Kth-term and advanced
// Bostan-Mori lectures in particular) refer to fast, NTT-based polynomial
// multiplication without building the transform itself. This fills that
// gap with a from-scratch, from-first-principles implementation.

// PDFBEGIN
const MOD: i64 = 998244353; // = 119 * 2^23 + 1, a common NTT-friendly prime

const PRIMITIVE_ROOT: i64 = 3; // a primitive root mod 998244353

fn mod_pow(mut base: i64, mut exp: u64, m: i64) -> i64 {
    let mut result = 1i64;
    base = base.rem_euclid(m);
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % m;
        }
        base = base * base % m;
        exp >>= 1;
    }
    result
}

fn mod_inverse(a: i64, m: i64) -> i64 {
    // Fermat's little theorem: a^(m-2) = a^-1 (mod m), for prime m.
    mod_pow(a, (m - 2) as u64, m)
}

/// In-place iterative NTT (Cooley-Tukey), length `a.len()` a power of 2.
/// `invert = false` is the forward transform; `invert = true` is the
/// inverse (including the final 1/n scaling).
fn ntt(a: &mut [i64], invert: bool) {
    let n = a.len();
    debug_assert!(n.is_power_of_two());

    // Bit-reversal permutation: puts the array in the order the iterative
    // butterfly network expects.
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }

    let mut len = 2usize;
    while len <= n {
        // A primitive `len`-th root of unity mod p (or its inverse).
        let mut w_len = mod_pow(PRIMITIVE_ROOT, ((MOD - 1) / len as i64) as u64, MOD);
        if invert {
            w_len = mod_inverse(w_len, MOD);
        }
        let mut i = 0;
        while i < n {
            let mut w = 1i64;
            for k in 0..len / 2 {
                let u = a[i + k];
                let v = a[i + k + len / 2] * w % MOD;
                a[i + k] = (u + v) % MOD;
                a[i + k + len / 2] = (u - v).rem_euclid(MOD);
                w = w * w_len % MOD;
            }
            i += len;
        }
        len <<= 1;
    }

    if invert {
        let n_inv = mod_inverse(n as i64, MOD);
        for x in a.iter_mut() {
            *x = *x * n_inv % MOD;
        }
    }
}

/// Polynomial multiplication mod MOD via NTT: O(n log n) instead of the
/// naive O(n^2).
fn multiply_via_ntt(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let result_len = a.len() + b.len() - 1;
    let mut size = 1usize;
    while size < result_len {
        size <<= 1;
    }

    let mut fa: Vec<i64> = a.iter().map(|&x| x.rem_euclid(MOD)).collect();
    let mut fb: Vec<i64> = b.iter().map(|&x| x.rem_euclid(MOD)).collect();
    fa.resize(size, 0);
    fb.resize(size, 0);

    ntt(&mut fa, false);
    ntt(&mut fb, false);
    for i in 0..size {
        fa[i] = fa[i] * fb[i] % MOD;
    }
    ntt(&mut fa, true);

    fa.truncate(result_len);
    fa
}

// ---- For comparison: the same transform over the complex numbers -------

/// In-place iterative complex FFT (Cooley-Tukey). Same bit-reversal and
/// butterfly structure as `ntt`; only the field differs. Twiddle factors
/// are evaluated directly with cos/sin (accuracy matters here, unlike in
/// the NTT where every root is exact).
fn fft(a: &mut [(f64, f64)], invert: bool) {
    let n = a.len();
    debug_assert!(n.is_power_of_two());

    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }

    let mut len = 2usize;
    while len <= n {
        let sign = if invert { -1.0 } else { 1.0 };
        let ang = sign * 2.0 * std::f64::consts::PI / len as f64;
        let mut i = 0;
        while i < n {
            for k in 0..len / 2 {
                let (wr, wi) = ((ang * k as f64).cos(), (ang * k as f64).sin());
                let (ur, ui) = a[i + k];
                let (xr, xi) = a[i + k + len / 2];
                let (vr, vi) = (xr * wr - xi * wi, xr * wi + xi * wr);
                a[i + k] = (ur + vr, ui + vi);
                a[i + k + len / 2] = (ur - vr, ui - vi);
            }
            i += len;
        }
        len <<= 1;
    }

    if invert {
        for x in a.iter_mut() {
            x.0 /= n as f64;
            x.1 /= n as f64;
        }
    }
}

/// Real polynomial product via the complex FFT. The result is only an
/// approximation of the true product; see the lecture for when rounding
/// it to integers is safe.
fn multiply_via_fft(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let result_len = a.len() + b.len() - 1;
    let mut size = 1usize;
    while size < result_len {
        size <<= 1;
    }
    let mut fa: Vec<(f64, f64)> = a.iter().map(|&x| (x, 0.0)).collect();
    let mut fb: Vec<(f64, f64)> = b.iter().map(|&x| (x, 0.0)).collect();
    fa.resize(size, (0.0, 0.0));
    fb.resize(size, (0.0, 0.0));
    fft(&mut fa, false);
    fft(&mut fb, false);
    for i in 0..size {
        let (ar, ai) = fa[i];
        let (br, bi) = fb[i];
        fa[i] = (ar * br - ai * bi, ar * bi + ai * br);
    }
    fft(&mut fa, true);
    fa.iter().take(result_len).map(|c| c.0).collect()
}

/// Product mod MOD using only floating-point FFTs, by splitting every
/// coefficient into 15-bit halves so each partial product stays far below
/// 2^53. Needs 4 forward and 3 inverse transforms (7 in total), against
/// 3 for `multiply_via_ntt`.
fn multiply_mod_via_fft_split(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let result_len = a.len() + b.len() - 1;
    let mut size = 1usize;
    while size < result_len {
        size <<= 1;
    }
    let split = |v: &[i64]| -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
        let mut lo = vec![(0.0, 0.0); size];
        let mut hi = vec![(0.0, 0.0); size];
        for (i, &x) in v.iter().enumerate() {
            let x = x.rem_euclid(MOD);
            lo[i] = ((x & 0x7FFF) as f64, 0.0);
            hi[i] = ((x >> 15) as f64, 0.0);
        }
        fft(&mut lo, false);
        fft(&mut hi, false);
        (lo, hi)
    };
    let (a0, a1) = split(a);
    let (b0, b1) = split(b);
    let mul = |x: &[(f64, f64)], y: &[(f64, f64)]| -> Vec<(f64, f64)> {
        x.iter().zip(y).map(|(&(ar, ai), &(br, bi))| (ar * br - ai * bi, ar * bi + ai * br)).collect()
    };
    let add = |x: Vec<(f64, f64)>, y: Vec<(f64, f64)>| -> Vec<(f64, f64)> {
        x.iter().zip(&y).map(|(&(ar, ai), &(br, bi))| (ar + br, ai + bi)).collect()
    };
    let mut c00 = mul(&a0, &b0);
    let mut c01 = add(mul(&a0, &b1), mul(&a1, &b0));
    let mut c11 = mul(&a1, &b1);
    fft(&mut c00, true);
    fft(&mut c01, true);
    fft(&mut c11, true);
    let shift15 = 1i64 << 15;
    let shift30 = (1i64 << 30) % MOD;
    (0..result_len)
        .map(|k| {
            let r00 = (c00[k].0.round() as i64).rem_euclid(MOD);
            let r01 = (c01[k].0.round() as i64).rem_euclid(MOD);
            let r11 = (c11[k].0.round() as i64).rem_euclid(MOD);
            (r00 + r01 * shift15 % MOD + r11 * shift30 % MOD) % MOD
        })
        .collect()
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

/// The definition, directly: O(n*m). Used only to cross-check the NTT.
fn multiply_naive(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        let av = av.rem_euclid(MOD);
        if av == 0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            out[i + j] = (out[i + j] + av * bv.rem_euclid(MOD)) % MOD;
        }
    }
    out
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
    // --- Sanity check: the roots of unity really are roots of unity ------
    let n = 8usize;
    let w = mod_pow(PRIMITIVE_ROOT, ((MOD - 1) / n as i64) as u64, MOD);
    let w_pow_n = mod_pow(w, n as u64, MOD);
    let w_pow_half = mod_pow(w, (n / 2) as u64, MOD);
    println!("primitive {n}-th root w={w}: w^{n} = {w_pow_n} (want 1), w^{} = {w_pow_half} (want {} = -1 mod p)", n / 2, MOD - 1);
    assert_eq!(w_pow_n, 1);
    assert_eq!(w_pow_half, MOD - 1);

    // --- NTT-based multiplication vs. the naive definition ---------------
    let mut rng = Rng(7);
    let mut mul_ok = true;
    for _ in 0..300 {
        let na = 1 + (rng.next() % 40) as usize;
        let nb = 1 + (rng.next() % 40) as usize;
        let a: Vec<i64> = (0..na).map(|_| rng.range(0, MOD - 1)).collect();
        let b: Vec<i64> = (0..nb).map(|_| rng.range(0, MOD - 1)).collect();
        let via_ntt = multiply_via_ntt(&a, &b);
        let naive = multiply_naive(&a, &b);
        if via_ntt != naive {
            mul_ok = false;
            println!("MULTIPLY MISMATCH na={na} nb={nb}: ntt={:?} naive={:?}", via_ntt, naive);
        }
    }
    println!("NTT-based polynomial multiplication matches the naive O(nm) definition: {}", if mul_ok { "PASS" } else { "FAIL" });

    // --- Round-trip: NTT then inverse NTT recovers the original ----------
    let mut roundtrip_ok = true;
    for _ in 0..100 {
        let k = 1 + (rng.next() % 8) as u32;
        let size = 1usize << k;
        let original: Vec<i64> = (0..size).map(|_| rng.range(0, MOD - 1)).collect();
        let mut a = original.clone();
        ntt(&mut a, false);
        ntt(&mut a, true);
        if a != original {
            roundtrip_ok = false;
            println!("ROUND-TRIP MISMATCH size={size}: got={:?} want={:?}", a, original);
        }
    }
    println!("forward NTT followed by inverse NTT recovers the original array: {}", if roundtrip_ok { "PASS" } else { "FAIL" });

    // --- A concrete worked example -----------------------------------------
    let a = vec![1i64, 2, 3]; // 1 + 2x + 3x^2
    let b = vec![4i64, 5, 6]; // 4 + 5x + 6x^2
    println!("\n(1+2x+3x^2)(4+5x+6x^2) via NTT   = {:?}", multiply_via_ntt(&a, &b));
    println!("(1+2x+3x^2)(4+5x+6x^2) via naive = {:?}", multiply_naive(&a, &b));
    // Expected by hand: 4 + 13x + 28x^2 + 27x^3 + 18x^4

    // --- NTT vs FFT: roots of unity in Z/17Z mirror the complex ones -----
    {
        let p = 17i64;
        let w = mod_pow(3, ((p - 1) / 8) as u64, p); // 3 is a primitive root mod 17
        let powers: Vec<i64> = (0..8).map(|k| mod_pow(w, k, p)).collect();
        let ok = w == 9 && powers == vec![1, 9, 13, 15, 16, 8, 4, 2] && mod_pow(w, 4, p) == p - 1 && mod_pow(w, 8, p) == 1;
        println!("\nroots of unity mod 17: w=9 has powers {:?}, w^4 = -1: {}", powers, if ok { "PASS" } else { "FAIL" });
    }

    // --- NTT vs FFT: when is rounding an FFT result safe? -----------------
    // Deterministic pseudo-random coefficients (same LCG in Rust and Python).
    fn lcg_values(seed: u64, count: usize, bound: u64) -> Vec<i64> {
        let mut state = seed;
        (0..count)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((state >> 33) % bound) as i64
            })
            .collect()
    }
    // Rounding noise makes exact counts of wrong coefficients chaotic once the
    // values pass 2^53, so report them by category only.
    fn describe(wrong: usize, total: usize) -> String {
        if wrong == 0 {
            "none".to_string()
        } else if wrong * 100 >= total * 99 {
            "essentially all (>= 99%)".to_string()
        } else {
            format!("about {}%", ((wrong * 100 + total / 2) / total + 5) / 10 * 10)
        }
    }
    let n = 1024usize;
    println!("\nplain FFT product of two length-{n} polynomials with coefficients in [0, B):");
    for bits in [8u32, 14, 18, 22, 30] {
        let bound = 1u64 << bits;
        let a = lcg_values(11 + bits as u64, n, bound);
        let b = lcg_values(97 + bits as u64, n, bound);
        let af: Vec<f64> = a.iter().map(|&x| x as f64).collect();
        let bf: Vec<f64> = b.iter().map(|&x| x as f64).collect();
        let approx = multiply_via_fft(&af, &bf);
        let mut exact = vec![0u128; 2 * n - 1];
        for i in 0..n {
            for j in 0..n {
                exact[i + j] += a[i] as u128 * b[j] as u128;
            }
        }
        let mut max_dev = 0.0f64;
        let mut wrong = 0usize;
        for k in 0..approx.len() {
            max_dev = max_dev.max((approx[k] - approx[k].round()).abs());
            if approx[k].round() as i128 != exact[k] as i128 {
                wrong += 1;
            }
        }
        let dev = if max_dev == 0.0 { "0".to_string() } else { format!("~1e{}", max_dev.log10().floor() as i32) };
        println!("  B = 2^{bits:<2}: max distance to nearest integer {dev:>6}, wrong coefficients: {}", describe(wrong, approx.len()));
    }

    // --- NTT vs FFT: arbitrary values mod p -------------------------------
    {
        let a = lcg_values(5, 1024, MOD as u64);
        let b = lcg_values(6, 1024, MOD as u64);
        let via_ntt = multiply_via_ntt(&a, &b);
        let af: Vec<f64> = a.iter().map(|&x| x as f64).collect();
        let bf: Vec<f64> = b.iter().map(|&x| x as f64).collect();
        let plain: Vec<i64> = multiply_via_fft(&af, &bf).iter().map(|x| (x.round() as i64).rem_euclid(MOD)).collect();
        let plain_wrong = plain.iter().zip(&via_ntt).filter(|(x, y)| x != y).count();
        let split = multiply_mod_via_fft_split(&a, &b);
        println!("\nfull-range coefficients mod p (length 1024 each):");
        println!("  plain FFT then reduce mod p: wrong residues: {}", describe(plain_wrong, via_ntt.len()));
        println!("  FFT with 15-bit splitting equals the NTT result: {}", if split == via_ntt { "PASS" } else { "FAIL" });
    }
}
