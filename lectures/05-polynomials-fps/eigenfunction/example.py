"""DP acceleration via eigenfunctions: some polynomial-valued DP
recurrences A_{i+1} = f(A_i) are linear maps f on the space of
polynomials, and if f happens to have a known eigenbasis -- functions v
with f(v) = lambda*v -- then M iterations collapse to "expand A_0 in that
basis, multiply each coefficient by lambda^M, recombine," turning an O(M)
loop into an O(log M) computation (plus the cost of the basis-change
itself).

This file fully implements and tests Problem 2 from the lecture, whose
eigenbasis is the elementary (1+x)^k (found via Taylor shift, no
composition machinery needed). For Problems 3 and 4, whose eigenbases
involve e^{kx} and x^k e^{-x}, it numerically verifies the eigenvalue
identities themselves (via a from-scratch power-series exp), including
the truncation subtlety the lecture calls out: e^{kx} mod x^{N+1} is only
an APPROXIMATE eigenfunction once you truncate too early, which is
exactly why those algorithms pad the working degree to N+M before
truncating back down.
"""

# PDFBEGIN
from __future__ import annotations

import random

MOD = 998244353

def mul(a: int, b: int) -> int:
    return (a % MOD) * (b % MOD) % MOD

def mod_pow(base: int, exp: int) -> int:
    return pow(base % MOD, exp, MOD)

def mod_inv(a: int) -> int:
    return mod_pow(a, MOD - 2)

# ---------------------------------------------------------------------
# Problem 2: A_{i+1} = (1+x) A_i'.  Eigenfunctions (1+x)^k, eigenvalue k.
# ---------------------------------------------------------------------


def binomial_table(n: int) -> list[list[int]]:
    """Pascal's triangle mod p, up to row n."""
    c = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        c[i][0] = 1
        for j in range(1, i + 1):
            c[i][j] = (c[i - 1][j - 1] + (c[i - 1][j] if j <= i - 1 else 0)) % MOD
    return c

def taylor_shift(a: list[int], t: int) -> list[int]:
    """Taylor shift: given P of degree <= n, compute Q(x) = P(x + t), O(n^2)."""
    n = len(a) - 1
    c = binomial_table(n)
    t_pow = [1] * (n + 1)
    for i in range(1, n + 1):
        t_pow[i] = mul(t_pow[i - 1], t)
    q = [0] * (n + 1)
    for i in range(n + 1):
        s = 0
        for j in range(i, n + 1):
            s = (s + mul(a[j], mul(c[j][i], t_pow[j - i]))) % MOD
        q[i] = s
    return q

def solve_problem2_eigenfunction(a0: list[int], m: int) -> list[int]:
    """Solve Problem 2 via the eigenfunction method: O(N^2) here (Taylor
    shift implemented naively); O(N log N) with a fast Taylor shift.
    """
    # Expand A_0(x) = sum_k c_k (1+x)^k: substitute x -> x-1, i.e. shift by -1.
    c = taylor_shift(a0, MOD - 1)
    c = [mul(ck, mod_pow(k, m)) for k, ck in enumerate(c)]
    # Recombine: A_M(x) = sum_k c_k (1+x)^k = D(x+1), shift by +1.
    return taylor_shift(c, 1)

def solve_problem2_direct(a0: list[int], m: int) -> list[int]:
    """The definition, directly: simulate the DP table. O(N*M). Used only
    to cross-check the eigenfunction-based solution above.
    """
    n = len(a0) - 1
    dp = a0[:]
    for _ in range(m):
        nxt = [0] * (n + 1)
        for j in range(n + 1):
            dp_j1 = dp[j + 1] if j + 1 <= n else 0
            nxt[j] = (mul(j, dp[j]) + mul(j + 1, dp_j1)) % MOD
        dp = nxt
    return dp

# ---------------------------------------------------------------------
# Power series exp, from scratch (O(n^2)), used to numerically verify the
# eigenfunction claims in Problems 3 and 4 without building a full
# composition/reversion engine.
# ---------------------------------------------------------------------


def power_series_exp(f: list[int], n: int) -> list[int]:
    """exp(f) mod x^n, given f[0] = 0. Standard E' = f' E recurrence:
    k*E_k = sum_{i=1}^{k} i*f_i*E_{k-i}.
    """
    assert (f[0] if f else 0) == 0
    e = [0] * n
    e[0] = 1
    for k in range(1, n):
        s = 0
        for i in range(1, min(k, len(f) - 1) + 1):
            s = (s + mul(mul(i, f[i]), e[k - i])) % MOD
        e[k] = mul(s, mod_inv(k))
    return e

def derivative(a: list[int]) -> list[int]:
    return [mul(i, a[i]) for i in range(1, len(a))]

def poly_add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % MOD for i in range(n)]

def poly_scale(a: list[int], s: int) -> list[int]:
    return [mul(x, s) for x in a]

def x_pow_k_times_exp_neg_x(k: int, length: int) -> list[int]:
    """x^k * e^{-x}, as a power series truncated to `length` terms."""
    f = [0] * max(length, 1)
    if length > 1:
        f[1] = MOD - 1  # f(x) = -x
    e = power_series_exp(f, max(length - k, 1))
    out = [0] * length
    for i, v in enumerate(e):
        if k + i < length:
            out[k + i] = v
    return out

# PDFEND -- everything below is the test harness and is not printed in the PDF

def main() -> None:
    # --- Problem 2: eigenfunction method vs. direct DP simulation --------
    rng = random.Random(31415)
    p2_ok = True
    for _ in range(200):
        n = rng.randint(1, 8)
        m = rng.randint(0, 20)
        a0 = [rng.randint(0, MOD - 1) for _ in range(n + 1)]
        via_eigen = solve_problem2_eigenfunction(a0, m)
        via_direct = solve_problem2_direct(a0, m)
        if via_eigen != via_direct:
            p2_ok = False
            print(f"PROBLEM2 MISMATCH a0={a0} m={m}: eigen={via_eigen} direct={via_direct}")
    print(f"Problem 2 (A_i+1=(1+x)A_i'): eigenfunction method matches direct DP simulation: {'PASS' if p2_ok else 'FAIL'}")

    # --- Problem 3's eigenfunction claim: e^{kx} has eigenvalue (k+1) ----
    p3_exact_ok = True
    p3_caveat_demonstrated = False
    for trial in range(50):
        k = trial % 7
        n = 6
        pad = 3
        f = [0] * (n + pad + 2)
        if len(f) > 1:
            f[1] = k
        a_full = power_series_exp(f, n + pad + 1)  # e^{kx} mod x^{n+pad+1}
        a_deriv = derivative(a_full)  # exact for indices 0..n+pad-1
        lhs = poly_add(a_full[: n + pad], a_deriv[: n + pad])
        rhs = poly_scale(a_full[: n + pad], k + 1)
        if lhs != rhs:
            p3_exact_ok = False
            print(f"PROBLEM3 EIGEN MISMATCH k={k}: lhs={lhs} rhs={rhs}")

        # Under-padded version: truncate A to exactly n+1 terms, compute
        # A' from THAT truncation (losing the true top term), and check
        # the identity fails at degree n but holds for degrees 0..n-1.
        a_short = a_full[: n + 1]
        a_short_deriv = derivative(a_short)  # length n, indices 0..n-1
        lhs_short = poly_add(a_short, a_short_deriv)
        rhs_short = poly_scale(a_short, k + 1)
        matches_low = lhs_short[:n] == rhs_short[:n]
        matches_top = lhs_short[n] == rhs_short[n]
        if matches_low and not matches_top:
            p3_caveat_demonstrated = True
        elif not matches_low:
            print(f"PROBLEM3 CAVEAT: unexpected low-degree mismatch k={k}")
    print(f"Problem 3 eigenfunction identity A+A'=(k+1)A holds exactly with enough padding: {'PASS' if p3_exact_ok else 'FAIL'}")
    print(f"...and (as the lecture warns) breaks specifically at the top degree without padding: {'CONFIRMED' if p3_caveat_demonstrated else 'NOT OBSERVED'}")

    # --- Problem 4's eigenfunction claim: x^k e^{-x} has eigenvalue k ----
    p4_ok = True
    for trial in range(50):
        k = trial % 5
        n = 6
        pad = 3
        length = k + n + pad + 2
        a = x_pow_k_times_exp_neg_x(k, length)
        a_deriv = derivative(a)  # length-1, exact throughout
        # f(A) = x*(A + A'): compute (A+A') then shift up by one degree (multiply by x).
        s = poly_add(a[: length - 1], a_deriv)
        f_a = [0] * length
        for i, v in enumerate(s):
            if i + 1 < length:
                f_a[i + 1] = v
        want = poly_scale(a, k)
        check_len = k + n
        if f_a[:check_len] != want[:check_len]:
            p4_ok = False
            print(f"PROBLEM4 EIGEN MISMATCH k={k}: f(A)={f_a[:check_len]} k*A={want[:check_len]}")
    print(f"Problem 4 eigenfunction identity x(A+A')=kA holds (for x^k e^-x) with enough padding: {'PASS' if p4_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
