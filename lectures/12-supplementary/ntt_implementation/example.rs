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
}
