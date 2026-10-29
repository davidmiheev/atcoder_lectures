"""The Number Theoretic Transform (NTT): a Discrete Fourier Transform
carried out in Z/pZ instead of C, for a prime p with p-1 divisible by a
large power of two. It gives the same O(n log n) fast-convolution trick
as the classical FFT, but with exact integer arithmetic -- no rounding
error, which is why it's the default for polynomial multiplication
modulo a prime in competitive programming.

This lecture is a supplement, not part of AtCoder's Algorithm Lectures
index: several lectures in this repository (the Kth-term and advanced
Bostan-Mori lectures in particular) refer to fast, NTT-based polynomial
multiplication without building the transform itself. This fills that
gap with a from-scratch, from-first-principles implementation.
"""

# PDFBEGIN
from __future__ import annotations

import math
import random

MOD = 998244353  # = 119 * 2^23 + 1, a common NTT-friendly prime

PRIMITIVE_ROOT = 3  # a primitive root mod 998244353

def mod_inverse(a: int, m: int) -> int:
    """Fermat's little theorem: a^(m-2) = a^-1 (mod m), for prime m."""
    return pow(a, m - 2, m)

def ntt(a: list[int], invert: bool) -> None:
    """In-place iterative NTT (Cooley-Tukey), length len(a) a power of 2.
    invert=False is the forward transform; invert=True is the inverse
    (including the final 1/n scaling).
    """
    n = len(a)
    assert n & (n - 1) == 0, "length must be a power of 2"

    # Bit-reversal permutation: puts the array in the order the iterative
    # butterfly network expects.
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]

    length = 2
    while length <= n:
        w_len = pow(PRIMITIVE_ROOT, (MOD - 1) // length, MOD)
        if invert:
            w_len = mod_inverse(w_len, MOD)
        i = 0
        while i < n:
            w = 1
            half = length // 2
            for k in range(half):
                u = a[i + k]
                v = a[i + k + half] * w % MOD
                a[i + k] = (u + v) % MOD
                a[i + k + half] = (u - v) % MOD
                w = w * w_len % MOD
            i += length
        length <<= 1

    if invert:
        n_inv = mod_inverse(n, MOD)
        for i in range(n):
            a[i] = a[i] * n_inv % MOD

def multiply_via_ntt(a: list[int], b: list[int]) -> list[int]:
    """Polynomial multiplication mod MOD via NTT: O(n log n) instead of
    the naive O(n^2)."""
    if not a or not b:
        return []
    result_len = len(a) + len(b) - 1
    size = 1
    while size < result_len:
        size <<= 1

    fa = [x % MOD for x in a] + [0] * (size - len(a))
    fb = [x % MOD for x in b] + [0] * (size - len(b))

    ntt(fa, False)
    ntt(fb, False)
    for i in range(size):
        fa[i] = fa[i] * fb[i] % MOD
    ntt(fa, True)

    return fa[:result_len]

# ---- For comparison: the same transform over the complex numbers -------

def fft(a: list[complex], invert: bool) -> None:
    """In-place iterative complex FFT (Cooley-Tukey). Same bit-reversal and
    butterfly structure as `ntt`; only the field differs. Twiddle factors
    are evaluated directly with cos/sin (accuracy matters here, unlike in
    the NTT where every root is exact)."""
    n = len(a)
    assert n & (n - 1) == 0, "length must be a power of 2"

    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]

    length = 2
    while length <= n:
        sign = -1.0 if invert else 1.0
        ang = sign * 2.0 * math.pi / length
        half = length // 2
        for i in range(0, n, length):
            for k in range(half):
                w = complex(math.cos(ang * k), math.sin(ang * k))
                u = a[i + k]
                v = a[i + k + half] * w
                a[i + k] = u + v
                a[i + k + half] = u - v
        length <<= 1

    if invert:
        for i in range(n):
            a[i] = complex(a[i].real / n, a[i].imag / n)

def multiply_via_fft(a: list[float], b: list[float]) -> list[float]:
    """Real polynomial product via the complex FFT. The result is only an
    approximation of the true product; see the lecture for when rounding
    it to integers is safe."""
    if not a or not b:
        return []
    result_len = len(a) + len(b) - 1
    size = 1
    while size < result_len:
        size <<= 1
    fa = [complex(x, 0.0) for x in a] + [0j] * (size - len(a))
    fb = [complex(x, 0.0) for x in b] + [0j] * (size - len(b))
    fft(fa, False)
    fft(fb, False)
    for i in range(size):
        fa[i] = fa[i] * fb[i]
    fft(fa, True)
    return [fa[i].real for i in range(result_len)]

def multiply_mod_via_fft_split(a: list[int], b: list[int]) -> list[int]:
    """Product mod MOD using only floating-point FFTs, by splitting every
    coefficient into 15-bit halves so each partial product stays far below
    2^53. Needs 4 forward and 3 inverse transforms (7 in total), against 3
    for `multiply_via_ntt`."""
    if not a or not b:
        return []
    result_len = len(a) + len(b) - 1
    size = 1
    while size < result_len:
        size <<= 1

    def split(v: list[int]) -> tuple[list[complex], list[complex]]:
        lo = [0j] * size
        hi = [0j] * size
        for i, x in enumerate(v):
            x %= MOD
            lo[i] = complex(x & 0x7FFF, 0.0)
            hi[i] = complex(x >> 15, 0.0)
        fft(lo, False)
        fft(hi, False)
        return lo, hi

    a0, a1 = split(a)
    b0, b1 = split(b)
    c00 = [a0[i] * b0[i] for i in range(size)]
    c01 = [a0[i] * b1[i] + a1[i] * b0[i] for i in range(size)]
    c11 = [a1[i] * b1[i] for i in range(size)]
    fft(c00, True)
    fft(c01, True)
    fft(c11, True)
    shift15 = 1 << 15
    shift30 = (1 << 30) % MOD
    out = []
    for k in range(result_len):
        r00 = round(c00[k].real) % MOD
        r01 = round(c01[k].real) % MOD
        r11 = round(c11[k].real) % MOD
        out.append((r00 + r01 * shift15 % MOD + r11 * shift30 % MOD) % MOD)
    return out

# PDFEND -- everything below is the test harness and is not printed in the PDF

def multiply_naive(a: list[int], b: list[int]) -> list[int]:
    """The definition, directly: O(n*m). Used only to cross-check the NTT."""
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        av %= MOD
        if av == 0:
            continue
        for j, bv in enumerate(b):
            out[i + j] = (out[i + j] + av * (bv % MOD)) % MOD
    return out

def main() -> None:
    # --- Sanity check: the roots of unity really are roots of unity ------
    n = 8
    w = pow(PRIMITIVE_ROOT, (MOD - 1) // n, MOD)
    w_pow_n = pow(w, n, MOD)
    w_pow_half = pow(w, n // 2, MOD)
    print(f"primitive {n}-th root w={w}: w^{n} = {w_pow_n} (want 1), w^{n // 2} = {w_pow_half} (want {MOD - 1} = -1 mod p)")
    assert w_pow_n == 1
    assert w_pow_half == MOD - 1

    # --- NTT-based multiplication vs. the naive definition ---------------
    rng = random.Random(7)
    mul_ok = True
    for _ in range(300):
        na = rng.randint(1, 40)
        nb = rng.randint(1, 40)
        a = [rng.randint(0, MOD - 1) for _ in range(na)]
        b = [rng.randint(0, MOD - 1) for _ in range(nb)]
        via_ntt = multiply_via_ntt(a, b)
        naive = multiply_naive(a, b)
        if via_ntt != naive:
            mul_ok = False
            print(f"MULTIPLY MISMATCH na={na} nb={nb}: ntt={via_ntt} naive={naive}")
    print(f"NTT-based polynomial multiplication matches the naive O(nm) definition: {'PASS' if mul_ok else 'FAIL'}")

    # --- Round-trip: NTT then inverse NTT recovers the original ----------
    roundtrip_ok = True
    for _ in range(100):
        k = rng.randint(1, 8)
        size = 1 << k
        original = [rng.randint(0, MOD - 1) for _ in range(size)]
        a = original[:]
        ntt(a, False)
        ntt(a, True)
        if a != original:
            roundtrip_ok = False
            print(f"ROUND-TRIP MISMATCH size={size}: got={a} want={original}")
    print(f"forward NTT followed by inverse NTT recovers the original array: {'PASS' if roundtrip_ok else 'FAIL'}")

    # --- A concrete worked example -----------------------------------------
    a = [1, 2, 3]  # 1 + 2x + 3x^2
    b = [4, 5, 6]  # 4 + 5x + 6x^2
    print(f"\n(1+2x+3x^2)(4+5x+6x^2) via NTT   = {multiply_via_ntt(a, b)}")
    print(f"(1+2x+3x^2)(4+5x+6x^2) via naive = {multiply_naive(a, b)}")
    # Expected by hand: 4 + 13x + 28x^2 + 27x^3 + 18x^4

    # --- NTT vs FFT: roots of unity in Z/17Z mirror the complex ones -----
    p = 17
    w = pow(3, (p - 1) // 8, p)  # 3 is a primitive root mod 17
    powers = [pow(w, k, p) for k in range(8)]
    ok = w == 9 and powers == [1, 9, 13, 15, 16, 8, 4, 2] and pow(w, 4, p) == p - 1 and pow(w, 8, p) == 1
    print(f"\nroots of unity mod 17: w=9 has powers {powers}, w^4 = -1: {'PASS' if ok else 'FAIL'}")

    # --- NTT vs FFT: when is rounding an FFT result safe? -----------------
    # Deterministic pseudo-random coefficients (same LCG in Rust and Python).
    def lcg_values(seed: int, count: int, bound: int) -> list[int]:
        state = seed
        out = []
        for _ in range(count):
            state = (state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
            out.append((state >> 33) % bound)
        return out

    # Rounding noise makes exact counts of wrong coefficients chaotic once the
    # values pass 2^53, so report them by category only.
    def describe(wrong: int, total: int) -> str:
        if wrong == 0:
            return "none"
        if wrong * 100 >= total * 99:
            return "essentially all (>= 99%)"
        return f"about {((wrong * 100 + total // 2) // total + 5) // 10 * 10}%"

    n = 1024
    print(f"\nplain FFT product of two length-{n} polynomials with coefficients in [0, B):")
    for bits in (8, 14, 18, 22, 30):
        bound = 1 << bits
        a = lcg_values(11 + bits, n, bound)
        b = lcg_values(97 + bits, n, bound)
        approx = multiply_via_fft([float(x) for x in a], [float(x) for x in b])
        exact = [0] * (2 * n - 1)
        for i in range(n):
            ai = a[i]
            for j in range(n):
                exact[i + j] += ai * b[j]
        max_dev = 0.0
        wrong = 0
        for k, x in enumerate(approx):
            max_dev = max(max_dev, abs(x - round(x)))
            if round(x) != exact[k]:
                wrong += 1
        dev = "0" if max_dev == 0.0 else f"~1e{math.floor(math.log10(max_dev))}"
        print(f"  B = 2^{bits:<2}: max distance to nearest integer {dev:>6}, wrong coefficients: {describe(wrong, len(approx))}")

    # --- NTT vs FFT: arbitrary values mod p -------------------------------
    a = lcg_values(5, 1024, MOD)
    b = lcg_values(6, 1024, MOD)
    via_ntt = multiply_via_ntt(a, b)
    plain = [round(x) % MOD for x in multiply_via_fft([float(x) for x in a], [float(x) for x in b])]
    plain_wrong = sum(1 for x, y in zip(plain, via_ntt) if x != y)
    split = multiply_mod_via_fft_split(a, b)
    print("\nfull-range coefficients mod p (length 1024 each):")
    print(f"  plain FFT then reduce mod p: wrong residues: {describe(plain_wrong, len(via_ntt))}")
    print(f"  FFT with 15-bit splitting equals the NTT result: {'PASS' if split == via_ntt else 'FAIL'}")

if __name__ == "__main__":
    main()
