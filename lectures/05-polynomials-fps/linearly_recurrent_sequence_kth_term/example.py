"""Kth term of a linearly recurrent sequence: two O(M(d) log K)-time
algorithms, both reducing the problem to fast exponentiation of a
degree-d object. Fiduccia works with x^K modulo the recurrence's
characteristic polynomial; Bostan-Mori works with the generating
function's rational form P(x)/Q(x) directly, halving the target index
each round via an even/odd split. Polynomial arithmetic here is the naive
O(d^2) kind, so the overall cost is O(d^2 log K), matching what the
lecture describes as the straightforward (non-FFT) implementation.
"""

# PDFBEGIN
from __future__ import annotations

MOD = 998244353

# ---------------------------------------------------------------------
# Fiduccia's algorithm
# ---------------------------------------------------------------------


def mul_x_mod_gamma(poly: list[int], c: list[int]) -> list[int]:
    """Multiply a degree-<d polynomial by x, reducing modulo the
    characteristic polynomial Gamma(x) = x^d - c[0]x^{d-1} - ... - c[d-1],
    using x^d = sum_j c[j] x^{d-1-j}. This single step is O(d).
    """
    d = len(poly)
    carry = poly[d - 1]
    out = [0] * d
    for i in range(d - 1, 0, -1):
        out[i] = poly[i - 1]
    for j in range(d):
        out[d - 1 - j] = (out[d - 1 - j] + carry * c[j]) % MOD
    return out

def mulmod_gamma(a: list[int], b: list[int], c: list[int]) -> list[int]:
    """Multiply two degree-<d polynomials and reduce modulo Gamma(x), by
    convolving then repeatedly substituting x^e (e >= d) with the
    recurrence's own right-hand side, from the top degree down.
    """
    d = len(a)
    raw = [0] * (2 * d - 1)
    for i in range(d):
        if a[i] == 0:
            continue
        for j in range(d):
            raw[i + j] = (raw[i + j] + a[i] * b[j]) % MOD
    for e in range(len(raw) - 1, d - 1, -1):
        coeff = raw[e]
        if coeff == 0:
            continue
        raw[e] = 0
        for j in range(d):
            raw[e - 1 - j] = (raw[e - 1 - j] + coeff * c[j]) % MOD
    return raw[:d]

def xk_mod_gamma(k: int, d: int, c: list[int]) -> list[int]:
    """x^k mod Gamma(x), via the doubling recursion in which the "multiply
    by x" step (for odd k) uses the cheap O(d) helper above.
    """
    if k < d:
        r = [0] * d
        r[k] = 1
        return r
    half = xk_mod_gamma(k // 2, d, c)
    r = mulmod_gamma(half, half, c)
    if k % 2 == 1:
        r = mul_x_mod_gamma(r, c)
    return r

def fiduccia(init: list[int], c: list[int], k: int) -> int:
    """A_K for a sequence satisfying A_i = c[0]A_{i-1} + ... +
    c[d-1]A_{i-d}, given initial values A_0..A_{d-1}.
    """
    d = len(init)
    if k < d:
        return init[k]
    r = xk_mod_gamma(k, d, c)
    return sum(r[i] * init[i] for i in range(d)) % MOD

# ---------------------------------------------------------------------
# Bostan-Mori's algorithm
# ---------------------------------------------------------------------


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        if av == 0:
            continue
        for j, bv in enumerate(b):
            out[i + j] = (out[i + j] + av * bv) % MOD
    return out

def negate_odd(a: list[int]) -> list[int]:
    return [(-v % MOD) if i % 2 == 1 else v for i, v in enumerate(a)]

def even_part(a: list[int]) -> list[int]:
    return a[0::2]

def odd_part(a: list[int]) -> list[int]:
    return a[1::2] if len(a) >= 2 else [0]

def bostan_mori(p: list[int], q: list[int], k: int) -> int:
    """[x^k] P(x)/Q(x), given Q(0) = 1."""
    while k > 0:
        qneg = negate_odd(q)
        u = poly_mul(p, qneg)
        w = poly_mul(q, qneg)
        p = even_part(u) if k % 2 == 0 else odd_part(u)
        q = even_part(w)
        k //= 2
    return p[0] if p else 0

def kth_term_bostan_mori(init: list[int], c: list[int], k: int) -> int:
    """Same recurrence as `fiduccia`, solved via the generating function
    A(x) = P(x)/Q(x) with Q(x) = 1 - c[0]x - ... - c[d-1]x^d.
    """
    d = len(init)
    q = [1] + [(-c[j]) % MOD for j in range(d)]
    full = poly_mul(init, q)
    p = full[:d]
    return bostan_mori(p, q, k)

# PDFEND -- everything below is the test harness and is not printed in the PDF

# ---------------------------------------------------------------------
# Direct O(K) unrolling, for cross-checking on small K.
# ---------------------------------------------------------------------


def direct(init: list[int], c: list[int], k: int) -> int:
    d = len(init)
    if k < d:
        return init[k]
    a = init[:]
    for i in range(d, k + 1):
        nxt = sum(c[j] * a[i - 1 - j] for j in range(d)) % MOD
        a.append(nxt)
    return a[k]

def main() -> None:
    # --- The lecture's own worked example: A_i = A_{i-1}+A_{i-2}+A_{i-3},
    # A_0=3, A_1=2, A_2=1. Expected A_5 = 16.
    init = [3, 2, 1]
    c = [1, 1, 1]
    a5_direct = direct(init, c, 5)
    a5_fid = fiduccia(init, c, 5)
    a5_bm = kth_term_bostan_mori(init, c, 5)
    print(f"worked example: A_5 = {a5_direct} (direct), {a5_fid} (Fiduccia), {a5_bm} (Bostan-Mori)")
    assert a5_direct == a5_fid == a5_bm == 16

    # --- Fibonacci, cross-checked against direct computation -----------
    fib_init = [0, 1]
    fib_c = [1, 1]
    fib_ok = True
    for k in range(60):
        d = direct(fib_init, fib_c, k)
        f = fiduccia(fib_init, fib_c, k)
        b = kth_term_bostan_mori(fib_init, fib_c, k)
        if d != f or d != b:
            fib_ok = False
            print(f"FIBONACCI MISMATCH k={k}: direct={d} fiduccia={f} bostan_mori={b}")
    print(f"Fibonacci cross-check (direct vs. both fast algorithms): {'PASS' if fib_ok else 'FAIL'}")

    # --- Random small-order recurrences, many random K -------------------
    import random

    rng = random.Random(42)
    random_ok = True
    for _ in range(300):
        d = rng.randint(1, 6)
        init_r = [rng.randint(0, 999) for _ in range(d)]
        c_r = [rng.randint(0, 999) for _ in range(d)]
        k = rng.randint(0, 499)
        dv = direct(init_r, c_r, k)
        fv = fiduccia(init_r, c_r, k)
        bv = kth_term_bostan_mori(init_r, c_r, k)
        if dv != fv or dv != bv:
            random_ok = False
            print(f"RANDOM MISMATCH d={d} k={k}: direct={dv} fiduccia={fv} bostan_mori={bv}")
    print(f"random small-order recurrences (direct vs. both fast algorithms): {'PASS' if random_ok else 'FAIL'}")

    # --- Huge K: the two fast algorithms must still agree with each other
    huge_k = 1_000_000_000_000
    f = fiduccia(fib_init, fib_c, huge_k)
    b = kth_term_bostan_mori(fib_init, fib_c, huge_k)
    print(f"\nFibonacci(10^12) mod {MOD}: Fiduccia = {f}, Bostan-Mori = {b}, agree = {f == b}")

if __name__ == "__main__":
    main()
