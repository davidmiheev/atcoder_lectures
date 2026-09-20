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
from __future__ import annotations

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


if __name__ == "__main__":
    main()
