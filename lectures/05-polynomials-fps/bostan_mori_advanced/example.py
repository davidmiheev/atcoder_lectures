"""Bostan-Mori, advanced: extracting a whole INTERVAL of coefficients
[x^[L,R]] P(x)/Q(x) in O(M(d) log R), not just a single one. This is the
"MSB-first" variant of Bostan-Mori: instead of halving the target index
from the bottom up (LSB-first, the ordinary Kth-term algorithm), it
halves the interval's right end from the top down. Applications: shifting
a linearly recurrent sequence's generating function, and computing
x^K mod Gamma(x) (the object at the heart of Fiduccia's algorithm) via
this machinery instead of binary exponentiation.
"""

# PDFBEGIN
from __future__ import annotations

MOD = 998244353

def add(a: int, b: int) -> int:
    return (a + b) % MOD

def sub(a: int, b: int) -> int:
    return (a - b) % MOD

def mul(a: int, b: int) -> int:
    return (a % MOD) * (b % MOD) % MOD

def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        if av == 0:
            continue
        for j, bv in enumerate(b):
            out[i + j] = add(out[i + j], mul(av, bv))
    return out

def negate_odd(a: list[int]) -> list[int]:
    return [sub(0, v) if i % 2 == 1 else v for i, v in enumerate(a)]

def even_part(a: list[int]) -> list[int]:
    return a[0::2]

def ceil_div(a: int, b: int) -> int:
    """ceil(a / b) for b > 0, correct for negative a."""
    return -(-a // b)

# Every interval-extraction result is returned as a PLAIN polynomial,
# reindexed to start at x^0: index j holds the true coefficient at degree
# l + j for the [l, r] interval that produced it. Once returned, the
# array is treated purely by its own 0-based index in any further
# multiplication -- exactly as the source's "[x^[L,R]] f(x)" notation
# intends (it produces a fresh polynomial, not a "shifted view" that
# remembers where it came from). Mixing up these two readings is the
# classic off-by-shift bug in this algorithm.


def range_coefficients_of_inverse(q: list[int], d: int, l: int, r: int) -> list[int]:
    """[x^[L,R]] 1/Q(x), for Q of degree <= d with Q(0) = 1, R >= 0 (L may
    be negative -- degrees below 0 are just 0, since 1/Q(x) is an
    ordinary power series). `q` must have exactly d+1 entries.

    This is the MSB-first Bostan-Mori recursion: at each step, halve R
    (from the top) using the same Q(x)Q(-x) = V(x^2) identity as the
    ordinary (LSB-first) algorithm, but keep track of a whole interval
    instead of a single target index.
    """
    if r == 0:
        # [x^[L,0]] 1/Q(x): only degree 0 can be nonzero (it's 1, since
        # Q(0)=1), and everything from L up to -1 is a negative-degree
        # coefficient, hence 0. So the answer is the single monomial
        # x^(-L) in the OUTPUT index, i.e. a 1 at output position -L.
        out = [0] * (0 - l + 1)
        out[-l] = 1
        return out

    l_prime = ceil_div(l - d, 2)
    r_prime = r // 2

    qneg = negate_odd(q)
    v = even_part(poly_mul(q, qneg))  # degree <= d, since q has d+1 entries

    # W(x) = [x^[L',R']] 1/V(x), a plain polynomial: W[t] = [x^(L'+t)] 1/V(x).
    w = range_coefficients_of_inverse(v, d, l_prime, r_prime)

    # S(x) = [x^[L-d,R]] 1/V(x^2), built from W: since 1/V(x^2) is even,
    # [x^(L-d+j)] 1/V(x^2) is nonzero only when L-d+j = 2*(L'+t) for some
    # t, and then equals W[t]. shift = 2L'-L+d is that j at t=0.
    s_len = r - (l - d) + 1
    s = [0] * s_len
    shift = 2 * l_prime - l + d
    for t, wt in enumerate(w):
        j = shift + 2 * t
        if 0 <= j < s_len:
            s[j] = wt

    # Q(-x) * S(x), as plain polynomials; then take LITERAL indices
    # [d, R-L+d] of that product (not a degree-shifted read).
    prod = poly_mul(qneg, s)
    out = [0] * (r - l + 1)
    for p in range(len(out)):
        idx = d + p
        if 0 <= idx < len(prod):
            out[p] = prod[idx]
    return out

def range_coefficients_of_rational(p: list[int], q: list[int], d: int, l: int, r: int) -> list[int]:
    """[x^[L,R]] P(x)/Q(x), for P of degree <= d-1, Q of degree <= d,
    Q(0) = 1. (Problem 1, reduced to Problem 2 above.)
    """
    # U(x) = [x^[L-d+1,R]] 1/Q(x): U[j] = [x^(L-d+1+j)] 1/Q(x).
    u = range_coefficients_of_inverse(q, d, l - d + 1, r)
    # P(x)*U(x): index m <-> true degree (L-d+1+m). Want true degree
    # L+p, i.e. m = p+d-1.
    prod = poly_mul(p, u)
    out = []
    for p_idx in range(r - l + 1):
        idx = p_idx + d - 1
        out.append(prod[idx] if 0 <= idx < len(prod) else 0)
    return out

def xk_mod_gamma_fiduccia(k: int, d: int, c: list[int]) -> list[int]:
    """x^K mod Gamma(x) via ordinary binary exponentiation (Fiduccia),
    used only to cross-check the Bostan-Mori-based method below.
    """

    def mul_x_mod_gamma(poly: list[int]) -> list[int]:
        d = len(poly)
        carry = poly[d - 1]
        out = [0] * d
        for i in range(d - 1, 0, -1):
            out[i] = poly[i - 1]
        for j in range(d):
            out[d - 1 - j] = add(out[d - 1 - j], mul(carry, c[j]))
        return out

    def mulmod_gamma(a: list[int], b: list[int]) -> list[int]:
        d = len(a)
        raw = [0] * (2 * d - 1)
        for i in range(d):
            for j in range(d):
                raw[i + j] = add(raw[i + j], mul(a[i], b[j]))
        for e in range(len(raw) - 1, d - 1, -1):
            coeff = raw[e]
            if coeff == 0:
                continue
            raw[e] = 0
            for j in range(d):
                raw[e - 1 - j] = add(raw[e - 1 - j], mul(coeff, c[j]))
        return raw[:d]

    def rec(k: int) -> list[int]:
        if k < d:
            r = [0] * d
            r[k] = 1
            return r
        half = rec(k // 2)
        r = mulmod_gamma(half, half)
        if k % 2 == 1:
            r = mul_x_mod_gamma(r)
        return r

    return rec(k)

def xk_mod_gamma_via_bostan_mori(q: list[int], d: int, k: int) -> list[int]:
    """x^K mod Gamma(x) via interval coefficient extraction (Problem 5):
    U(x) = [x^[K-d+1,K]] 1/Q(x), then reverse (Q(x)U(x) mod x^d).
    """
    u = range_coefficients_of_inverse(q, d, k - d + 1, k)
    qu = poly_mul(q, u)
    r_tilde = qu[:d]
    return r_tilde[::-1]

# PDFEND -- everything below is the test harness and is not printed in the PDF

# --- Reference implementations, for cross-checking ----------------------


def brute_force_inverse(q: list[int], up_to: int) -> list[int]:
    """O(d^2 * N): the first up_to+1 coefficients of 1/Q(x) via the
    standard recurrence, given Q(0) = 1.
    """
    inv = [0] * (up_to + 1)
    inv[0] = 1
    for k in range(1, up_to + 1):
        s = 0
        for i in range(1, min(k, len(q) - 1) + 1):
            s = add(s, mul(q[i], inv[k - i]))
        inv[k] = sub(0, s)
    return inv

def brute_force_range_of_inverse(q: list[int], l: int, r: int) -> list[int]:
    inv = brute_force_inverse(q, max(r, 0))
    return [0 if deg < 0 else inv[deg] for deg in range(l, r + 1)]

def direct_sequence(p: list[int], q: list[int], up_to: int) -> list[int]:
    """Direct recurrence unrolling of A_i = c_1 A_{i-1} + ... + c_d A_{i-d},
    from A(x) = P(x)/Q(x) with Q(x) = 1 - c_1 x - ... - c_d x^d.
    """
    d = len(q) - 1
    c = [sub(0, q[j]) for j in range(1, d + 1)]
    a = [0] * (up_to + 1)
    for i in range(up_to + 1):
        p_i = p[i] if i < len(p) else 0
        if i < d:
            # A(x)Q(x) = P(x), Q(0)=1: A_i = P_i - sum_{j<i} A_j Q_{i-j}.
            s = p_i % MOD
            for j in range(i):
                s = sub(s, mul(a[j], q[i - j]))
            a[i] = s
        else:
            v = 0
            for j in range(d):
                v = add(v, mul(c[j], a[i - 1 - j]))
            a[i] = v
    return a

class Rng:
    def __init__(self, seed: int) -> None:
        self.state = seed

    def next(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        return z ^ (z >> 31)

    def range(self, lo: int, hi: int) -> int:
        return lo + (self.next() % (hi - lo + 1))

def main() -> None:
    rng = Rng(2025)

    # --- Core: [x^[L,R]] 1/Q(x), including negative L, vs. brute force --
    inv_ok = True
    for _ in range(300):
        d = 1 + rng.next() % 5
        q = [1] + [rng.range(-5, 5) for _ in range(d)]
        r = rng.range(0, 20)
        l = rng.range(-10, r)
        got = range_coefficients_of_inverse(q, d, l, r)
        want = brute_force_range_of_inverse(q, l, r)
        if got != want:
            inv_ok = False
            print(f"INVERSE MISMATCH q={q} l={l} r={r}: got={got} want={want}")
    print(f"[x^[L,R]] 1/Q(x) matches brute-force power series inversion (L possibly negative): {'PASS' if inv_ok else 'FAIL'}")

    # --- P(x)/Q(x) interval extraction vs. direct recurrence unrolling --
    rational_ok = True
    for _ in range(300):
        d = 1 + rng.next() % 5
        q = [1] + [rng.range(-5, 5) for _ in range(d)]
        p = [rng.range(-5, 5) for _ in range(d)]
        r = rng.range(0, 25)
        l = rng.range(0, r)  # keep L >= 0 to compare against the sequence directly
        got = range_coefficients_of_rational(p, q, d, l, r)
        seq = direct_sequence(p, q, r)
        want = [seq[i] for i in range(l, r + 1)]
        if got != want:
            rational_ok = False
            print(f"RATIONAL MISMATCH p={p} q={q} l={l} r={r}: got={got} want={want}")
    print(f"[x^[L,R]] P(x)/Q(x) matches direct recurrence unrolling: {'PASS' if rational_ok else 'FAIL'}")

    # --- Problem 5: x^K mod Gamma(x), Bostan-Mori vs. Fiduccia -----------
    modgamma_ok = True
    for _ in range(200):
        d = 1 + rng.next() % 5
        q = [1] + [rng.range(-5, 5) for _ in range(d)]
        c = [sub(0, q[j]) for j in range(1, d + 1)]
        k = rng.range(0, 2000)
        via_fiduccia = xk_mod_gamma_fiduccia(k, d, c)
        via_bostan_mori = xk_mod_gamma_via_bostan_mori(q, d, k)
        if via_fiduccia != via_bostan_mori:
            modgamma_ok = False
            print(f"MODGAMMA MISMATCH q={q} k={k}: fiduccia={via_fiduccia} bostan_mori={via_bostan_mori}")
    print(f"x^K mod Gamma(x): interval-extraction method matches Fiduccia's binary exponentiation: {'PASS' if modgamma_ok else 'FAIL'}")

    # --- Problem 4: shifting a linearly recurrent sequence ---------------
    shift_ok = True
    for _ in range(150):
        d = 1 + rng.next() % 4
        q = [1] + [rng.range(-4, 4) for _ in range(d)]
        p = [rng.range(-4, 4) for _ in range(d)]
        k = rng.range(0, 100)

        # T = [x^[K,K+d-1]] P(x)/Q(x); S = (T * Q) truncated to degree d-1.
        t = range_coefficients_of_rational(p, q, d, k, k + d - 1)
        s = poly_mul(t, q)[:d]

        # S(x)/Q(x) should generate (A_K, A_{K+1}, ...).
        shifted_seq = direct_sequence(s, q, 2 * d)
        original_seq = direct_sequence(p, q, k + 2 * d)
        want = [original_seq[k + i] for i in range(2 * d + 1)]
        if shifted_seq != want:
            shift_ok = False
            print(f"SHIFT MISMATCH p={p} q={q} k={k}: shifted={shifted_seq} want={want}")
    print(f"shifted generating function S(x)/Q(x) reproduces (A_K, A_K+1, ...): {'PASS' if shift_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
