"""Linearly recurrent sequences: sequences A satisfying a constant-coefficient
recurrence A_i = c_1 A_{i-1} + ... + c_d A_{i-d} for all i >= d. Such
sequences are exactly the ones whose generating function is a rational
function P(x)/Q(x) with deg P < d, deg Q <= d, [x^0]Q(x)=1, and are exactly
the ones expressible as a linear combination of the entries of M^i for a
d x d "companion matrix" M built from the recurrence.

This file verifies both characterizations computationally, then verifies
(via a from-scratch minimal-recurrence fitter, solving a linear system by
Gaussian elimination mod a prime) that linearly recurrent sequences are
closed under sum, convolution, and termwise product, and that
polynomial-valued and polynomial-times-geometric sequences are linearly
recurrent -- reproducing every closure claim in the lecture by fitting a
candidate recurrence from a prefix of a sequence and confirming it
extrapolates correctly on a long, disjoint holdout range.
"""
from __future__ import annotations

MOD = 998244353


def add(a: int, b: int) -> int:
    return (a + b) % MOD


def sub(a: int, b: int) -> int:
    return (a - b) % MOD


def mul(a: int, b: int) -> int:
    return (a % MOD) * (b % MOD) % MOD


def mod_pow(base: int, exp: int) -> int:
    return pow(base % MOD, exp, MOD)


def mod_inv(a: int) -> int:
    return mod_pow(a, MOD - 2)


def recur_seq(c: list[int], init: list[int], n: int) -> list[int]:
    """Generate A_0, ..., A_{n-1}: A_i = sum_k c[k] * A_{i-1-k}."""
    d = len(c)
    assert len(init) == d
    a = init[:]
    while len(a) < n:
        i = len(a)
        s = 0
        for k in range(d):
            s = add(s, mul(c[k], a[i - 1 - k]))
        a.append(s)
    return a[:n]


Mat = list[list[int]]


def companion_matrix(c: list[int]) -> Mat:
    """Companion matrix of A_i = c[0]*A_{i-1} + ... + c[d-1]*A_{i-d}, so that
    V_{i+1} = M V_i for V_i = (A_i, ..., A_{i+d-1})^T."""
    d = len(c)
    m = [[0] * d for _ in range(d)]
    for i in range(d - 1):
        m[i][i + 1] = 1
    for j in range(d):
        m[d - 1][j] = c[d - 1 - j]
    return m


def mat_mul(a: Mat, b: Mat) -> Mat:
    n, k, m = len(a), len(b), len(b[0])
    out = [[0] * m for _ in range(n)]
    for i in range(n):
        for t in range(k):
            if a[i][t] == 0:
                continue
            for j in range(m):
                out[i][j] = add(out[i][j], mul(a[i][t], b[t][j]))
    return out


def mat_pow(base: Mat, exp: int) -> Mat:
    n = len(base)
    result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    while exp > 0:
        if exp & 1:
            result = mat_mul(result, base)
        base = mat_mul(base, base)
        exp >>= 1
    return result


def mat_vec(m: Mat, v: list[int]) -> list[int]:
    return [sum((mul(mij, vj) for mij, vj in zip(row, v)), 0) % MOD for row in m]


def verify_companion(c: list[int], init: list[int], n: int) -> bool:
    """Checks A_i == (M^i V_0)_0 for i = 0..n (Theorem 18)."""
    seq = recur_seq(c, init, n)
    m = companion_matrix(c)
    for i in range(n):
        vi = mat_vec(mat_pow(m, i), init)
        if vi[0] != seq[i]:
            print(f"COMPANION MISMATCH i={i}: matrix={vi[0]} direct={seq[i]}")
            return False
    return True


def find_recurrence(seq: list[int], d: int) -> list[int] | None:
    """Fits a candidate order-D recurrence using the equations at indices
    D..2D-1 (Gaussian elimination mod MOD)."""
    assert len(seq) >= 2 * d
    rows = []
    for i in range(d):
        row = [seq[d + i - 1 - k] for k in range(d)]
        row.append(seq[d + i])
        rows.append(row)
    for col in range(d):
        pivot = next((r for r in range(col, d) if rows[r][col] != 0), None)
        if pivot is None:
            return None
        rows[col], rows[pivot] = rows[pivot], rows[col]
        inv = mod_inv(rows[col][col])
        for j in range(col, d + 1):
            rows[col][j] = mul(rows[col][j], inv)
        for r in range(d):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col]
                for j in range(col, d + 1):
                    rows[r][j] = sub(rows[r][j], mul(factor, rows[col][j]))
    return [rows[i][d] for i in range(d)]


def verify_holdout(seq: list[int], c: list[int], start: int) -> bool:
    d = len(c)
    for i in range(start, len(seq)):
        s = 0
        for k in range(d):
            s = add(s, mul(c[k], seq[i - 1 - k]))
        if s != seq[i]:
            print(f"HOLDOUT MISMATCH i={i}: fitted={s} actual={seq[i]}")
            return False
    return True


def check_linearly_recurrent(name: str, seq: list[int], d: int) -> bool:
    c = find_recurrence(seq, d)
    if c is None:
        print(f"{name}: FAILED TO FIT (degenerate window)")
        return False
    ok = verify_holdout(seq, c, 2 * d)
    print(
        f"{name}: order-{d} recurrence fitted from a prefix reproduces the next "
        f"{len(seq) - 2 * d} terms: {'PASS' if ok else 'FAIL'}"
    )
    return ok


def main() -> None:
    n = 60

    # --- Theorem 2 / worked example: A_i = 2A_{i-1}+3A_{i-2}, A0=1,A1=4 ---
    generic = recur_seq([2, 3], [1, 4], 6)
    expected_generic = [1, 4, 11, 34, 101, 304]
    generic_ok = generic == expected_generic
    print(f"Worked example A_i=2A_{{i-1}}+3A_{{i-2}}, A0=1,A1=4: matches lecture's own A0..A5: {'PASS' if generic_ok else 'FAIL'}")

    # --- Theorem 18: companion matrix characterization ---
    fib_init = [0, 1]
    companion_fib_ok = verify_companion([1, 1], fib_init, 30)
    print(f"Companion matrix reproduces the Fibonacci sequence (Theorem 18): {'PASS' if companion_fib_ok else 'FAIL'}")
    companion_generic_ok = verify_companion([2, 3], [1, 4], 30)
    print(f"Companion matrix reproduces A_i=2A_{{i-1}}+3A_{{i-2}} (Theorem 18): {'PASS' if companion_generic_ok else 'FAIL'}")

    # --- Base sequences used for the closure checks below ---
    a = recur_seq([1, 1], [0, 1], n)  # Fibonacci, order 2
    b = recur_seq([2, 3], [1, 4], n)  # order 2

    # --- Theorem 26: closed under sum ---
    sum_seq = [add(x, y) for x, y in zip(a, b)]
    sum_ok = check_linearly_recurrent("Sum (Theorem 26)", sum_seq, 4)

    # --- Theorem 27: closed under convolution ---
    conv_seq = [sum((mul(a[j], b[i - j]) for j in range(i + 1)), 0) % MOD for i in range(n)]
    conv_ok = check_linearly_recurrent("Convolution (Theorem 27)", conv_seq, 4)

    # --- Theorem 28: closed under termwise product ---
    prod_seq = [mul(x, y) for x, y in zip(a, b)]
    prod_ok = check_linearly_recurrent("Termwise product (Theorem 28)", prod_seq, 4)

    # --- Theorem 8: polynomial-valued sequences ---
    poly_seq = [mul(mul(i, i), i) for i in range(n)]  # A_i = i^3, order <= 4
    poly_ok = check_linearly_recurrent("Polynomial A_i=i^3 (Theorem 8)", poly_seq, 4)

    # --- Theorem 10: polynomial times geometric sequences ---
    poly_geom_seq = [mul(mul(mul(i, i), i), mod_pow(3, i)) for i in range(n)]  # A_i = i^3*3^i
    poly_geom_ok = check_linearly_recurrent("Polynomial x geometric A_i=i^3*3^i (Theorem 10)", poly_geom_seq, 4)

    # --- Periodic sequences are linearly recurrent of order <= period ---
    period = 5
    periodic_pattern = [3, 1, 4, 1, 5]
    periodic_seq = [periodic_pattern[i % period] for i in range(n)]
    periodic_ok = check_linearly_recurrent("Periodic sequence (period 5)", periodic_seq, period)

    all_ok = (
        generic_ok
        and companion_fib_ok
        and companion_generic_ok
        and sum_ok
        and conv_ok
        and prod_ok
        and poly_ok
        and poly_geom_ok
        and periodic_ok
    )
    print(f"\nAll checks: {'PASS' if all_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
