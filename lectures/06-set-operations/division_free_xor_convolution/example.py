"""XOR convolution without division: computing C_U = sum_{S xor T = U} A_S B_T
for length-2^N arrays using only ring operations (add, subtract, multiply)
-- no multiplication by an inverse of 2 or any other element, unlike the
standard Walsh-Hadamard-transform solution. This matters whenever the
coefficient ring has no inverse for 2 (e.g. Z/2^k Z, used below
specifically because 2 is a zero divisor there and so provably has no
inverse), or more generally whenever "no division" is a requirement.

The approach: view A, B as multilinear polynomials in commuting variables
x_0,...,x_{N-1} (A(x) = sum_S A_S * prod_{i in S} x_i). XOR convolution is
exactly polynomial multiplication mod (x_i^2 - 1, all i). Substituting
x_i = y_i - 1 turns each relation into y_i^2 - 2y_i, and multiplying mod
(y_i^2 - 2y_i, all i) is a "ranked OR-convolution": the same machinery
used for subset convolution, tracking |S|+|T| as an extra polynomial
variable t so that y_i^2 = 2y_i can be applied via integer powers of 2 (an
ordinary ring element, needing no division) at the end rather than during
the convolution itself. Converting between the x- and y-bases is a
per-variable change of basis, computable with N passes of a
SOS-transform-shaped loop.

All arithmetic below is mod 2^20 -- a power of two specifically because 2
is then a zero divisor (2 * 2^19 = 0) and therefore has no multiplicative
inverse, so a Walsh-Hadamard-transform-based solution (which divides by 2
up to N times) could not be used in this ring at all, while the algorithm
here works unmodified.
"""

# PDFBEGIN
from __future__ import annotations

import random

MOD_BITS = 20

MOD = 1 << MOD_BITS

MASK = MOD - 1

def add(a: int, b: int) -> int:
    return (a + b) & MASK

def sub(a: int, b: int) -> int:
    return (a - b) & MASK

def mul(a: int, b: int) -> int:
    return ((a & MASK) * (b & MASK)) & MASK

def subset_zeta(n: int, a: list[int]) -> list[int]:
    """Subset zeta transform (SOS): B_T = sum_{S subseteq T} A_S. O(n*2^n)."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if (t >> i) & 1:
                a[t] = add(a[t], a[t ^ (1 << i)])
    return a

def subset_mobius(n: int, a: list[int]) -> list[int]:
    """Subset Mobius transform, the inverse of the above. O(n*2^n)."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if (t >> i) & 1:
                a[t] = sub(a[t], a[t ^ (1 << i)])
    return a

def transform_x_to_y(n: int, a: list[int]) -> list[int]:
    """Converts coefficients from the x-basis to the y-basis under
    x_i = y_i - 1: a0 + a1*x_i = (a0-a1) + a1*y_i, so (a0,a1) -> (a0-a1,a1)."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if not (t >> i) & 1:
                hi = t | (1 << i)
                a[t] = sub(a[t], a[hi])
    return a

def transform_y_to_x(n: int, a: list[int]) -> list[int]:
    """Inverse of transform_x_to_y, under y_i = x_i + 1:
    b0 + b1*y_i = (b0+b1) + b1*x_i, so (b0,b1) -> (b0+b1,b1)."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if not (t >> i) & 1:
                hi = t | (1 << i)
                a[t] = add(a[t], a[hi])
    return a

def ranked_or_convolution_pow2(n: int, a1: list[int], b1: list[int]) -> list[int]:
    """Problem 2: C_U = sum_{S union T = U} A_S*B_T*2^(|S|+|T|-|U|), i.e.
    polynomial multiplication mod (y_i^2 - 2y_i, all i). Solved exactly
    like subset convolution, but keeping every rank r=|S|+|T| (not just
    r=|U|) so the final 2^(|S|+|T|-|U|) factor can be folded in afterward
    using only integer powers of 2. O(n^2*2^n)."""
    size = 1 << n

    za = [[0] * size for _ in range(n + 1)]
    zb = [[0] * size for _ in range(n + 1)]
    for s in range(size):
        k = bin(s).count("1")
        za[k][s] = a1[s]
        zb[k][s] = b1[s]
    for k in range(n + 1):
        za[k] = subset_zeta(n, za[k])
        zb[k] = subset_zeta(n, zb[k])

    conv = [[0] * (2 * n + 1) for _ in range(size)]
    for u in range(size):
        for k in range(n + 1):
            zak = za[k][u]
            if zak == 0:
                continue
            for l in range(n + 1):
                zbl = zb[l][u]
                if zbl == 0:
                    continue
                conv[u][k + l] = add(conv[u][k + l], mul(zak, zbl))

    chat = [None] * (2 * n + 1)
    for r in range(2 * n + 1):
        col = [conv[u][r] for u in range(size)]
        chat[r] = subset_mobius(n, col)

    c1 = [0] * size
    for u in range(size):
        pu = bin(u).count("1")
        acc = 0
        pow2 = 1
        for r in range(pu, 2 * n + 1):
            acc = add(acc, mul(chat[r][u], pow2))
            pow2 = mul(pow2, 2)
        c1[u] = acc
    return c1

def xor_convolution_division_free(n: int, a: list[int], b: list[int]) -> list[int]:
    """The full division-free XOR convolution algorithm."""
    a1 = transform_x_to_y(n, a)
    b1 = transform_x_to_y(n, b)
    c1 = ranked_or_convolution_pow2(n, a1, b1)
    return transform_y_to_x(n, c1)

# PDFEND -- everything below is the test harness and is not printed in the PDF

def xor_convolution_naive(n: int, a: list[int], b: list[int]) -> list[int]:
    """Brute-force O(4^n) reference for XOR convolution."""
    size = 1 << n
    c = [0] * size
    for s in range(size):
        for t in range(size):
            c[s ^ t] = add(c[s ^ t], mul(a[s], b[t]))
    return c

def problem2_naive(n: int, a: list[int], b: list[int]) -> list[int]:
    """Brute-force O(4^n) reference for Problem 2."""
    size = 1 << n
    c = [0] * size
    for s in range(size):
        for t in range(size):
            u = s | t
            w = bin(s).count("1") + bin(t).count("1") - bin(u).count("1")
            pow2 = 1
            for _ in range(w):
                pow2 = mul(pow2, 2)
            c[u] = add(c[u], mul(mul(a[s], b[t]), pow2))
    return c

def main() -> None:
    rng = random.Random(2026091901)

    # --- Sanity check: transform_y_to_x undoes transform_x_to_y exactly ---
    roundtrip_ok = True
    for _ in range(100):
        n = rng.randint(0, 8)
        size = 1 << n
        a = [rng.randint(0, MOD - 1) for _ in range(size)]
        y = transform_x_to_y(n, a)
        back = transform_y_to_x(n, y)
        if back != a:
            roundtrip_ok = False
            print(f"BASIS ROUND-TRIP MISMATCH n={n}")
    print(f"Basis change x<->y round-trips exactly, over 100 random trials: {'PASS' if roundtrip_ok else 'FAIL'}")

    # --- Problem 2 in isolation: ranked OR-convolution vs. brute force ---
    problem2_ok = True
    for _ in range(100):
        n = rng.randint(0, 6)
        size = 1 << n
        a = [rng.randint(0, MOD - 1) for _ in range(size)]
        b = [rng.randint(0, MOD - 1) for _ in range(size)]
        via_ranked = ranked_or_convolution_pow2(n, a, b)
        via_naive = problem2_naive(n, a, b)
        if via_ranked != via_naive:
            problem2_ok = False
            print(f"PROBLEM2 MISMATCH n={n}")
    print(f"Problem 2 (ranked OR-convolution with 2^(|S|+|T|-|U|) weights) matches brute force, over 100 random trials: {'PASS' if problem2_ok else 'FAIL'}")

    # --- Full pipeline: division-free XOR convolution vs. brute force ---
    full_ok = True
    for _ in range(100):
        n = rng.randint(0, 6)
        size = 1 << n
        a = [rng.randint(0, MOD - 1) for _ in range(size)]
        b = [rng.randint(0, MOD - 1) for _ in range(size)]
        via_division_free = xor_convolution_division_free(n, a, b)
        via_naive = xor_convolution_naive(n, a, b)
        if via_division_free != via_naive:
            full_ok = False
            print(f"XOR CONVOLUTION MISMATCH n={n}")
    print(
        f"Division-free XOR convolution (mod 2^{MOD_BITS}, where 2 is a zero divisor with no inverse) matches brute force, over 100 random trials: {'PASS' if full_ok else 'FAIL'}"
    )

    all_ok = roundtrip_ok and problem2_ok and full_ok
    print(f"\nAll checks: {'PASS' if all_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
