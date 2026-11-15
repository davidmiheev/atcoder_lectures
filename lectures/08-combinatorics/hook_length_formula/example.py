"""The hook length formula: for a partition lambda of n, the number of
standard Young tableaux of shape lambda is

    f(lambda) = n! / prod over cells of (hook length of the cell).

The core below computes hook lengths, the count (exactly, and modulo a
prime), the same count from the determinant-derived product formula
(Frobenius), and a uniform random tableau generator built on the hook
walk of the first proof. The test harness after PDFEND cross-checks all
of them against each other and against definition-level brute force.
"""

# PDFBEGIN
from __future__ import annotations

from math import factorial
from typing import Callable

P = 998_244_353

# A partition: row lengths lam[0] >= lam[1] >= ... > 0 (row 0 is the top row).
Shape = list[int]


def conjugate(lam: Shape) -> list[int]:
    """Column lengths (the conjugate partition)."""
    width = lam[0] if lam else 0
    return [sum(1 for length in lam if length > c) for c in range(width)]


def hook_lengths(lam: Shape) -> list[list[int]]:
    """hook[r][c] = arm + leg + 1, where the arm is the number of cells to
    the right of (r, c) in its row and the leg the number of cells below it
    in its column."""
    conj = conjugate(lam)
    return [[(lam[r] - c - 1) + (conj[c] - r - 1) + 1 for c in range(lam[r])] for r in range(len(lam))]


def count_syt_exact(lam: Shape) -> int:
    """Exact f(lambda) by the hook length formula."""
    denominator = 1
    for row in hook_lengths(lam):
        for h in row:
            denominator *= h
    numerator = factorial(sum(lam))
    assert numerator % denominator == 0, "hook product must divide n!"
    return numerator // denominator


def count_syt_mod(lam: Shape) -> int:
    """f(lambda) modulo the prime P, for any n < P."""
    numerator = factorial(sum(lam)) % P
    denominator = 1
    for row in hook_lengths(lam):
        for h in row:
            denominator = denominator * h % P
    return numerator * pow(denominator, P - 2, P) % P


def count_syt_frobenius_mod(lam: Shape) -> int:
    """f(lambda) from Frobenius' product formula
        n! * prod_{i<j}(b_i - b_j) / prod_i b_i!,   b_i = lam_i + m - 1 - i,
    which the determinant proof produces directly (rows indexed from 0)."""
    m = len(lam)
    b = [lam[i] + m - 1 - i for i in range(m)]
    numerator = factorial(sum(lam)) % P
    denominator = 1
    for i in range(m):
        for j in range(i + 1, m):
            numerator = numerator * (b[i] - b[j]) % P
        denominator = denominator * (factorial(b[i]) % P) % P
    return numerator * pow(denominator, P - 2, P) % P


def sample_tableau(lam: Shape, rand_below: Callable[[int], int]) -> list[list[int]]:
    """A uniformly random standard Young tableau of shape `lam`, built by
    repeatedly running a hook walk from a uniformly random cell: the corner
    where the walk stops receives the largest remaining entry.
    `rand_below(k)` must return a uniform integer in range(k)."""
    cur = list(lam)
    tableau = [[0] * length for length in lam]
    for value in range(sum(lam), 0, -1):
        # Pick a uniformly random cell of the current shape.
        pick = rand_below(sum(cur))
        r = 0
        while pick >= cur[r]:
            pick -= cur[r]
            r += 1
        c = pick
        # Hook walk: jump to a uniformly random other cell of the hook until a corner.
        while True:
            arm = cur[r] - c - 1
            leg = sum(1 for length in cur[r + 1:] if length > c)
            if arm + leg == 0:
                break
            k = rand_below(arm + leg)
            if k < arm:
                c += 1 + k
            else:
                r += 1 + (k - arm)
        tableau[r][c] = value
        cur[r] -= 1
        if cur[r] == 0:
            cur.pop()
    return tableau

# PDFEND -- everything below is the test harness and is not printed in the PDF


def partitions(n: int) -> list[Shape]:
    """All partitions of n in non-increasing lexicographic order."""
    out: list[Shape] = []

    def go(remaining: int, max_part: int, current: Shape) -> None:
        if remaining == 0:
            out.append(list(current))
            return
        for part in range(min(remaining, max_part), 0, -1):
            current.append(part)
            go(remaining - part, part, current)
            current.pop()

    go(n, n, [])
    return out


def count_syt_recurrence(lam: Shape, memo: dict[tuple[int, ...], int]) -> int:
    """f(lambda) = sum of f(lambda - corner) over removable corners."""
    if not lam:
        return 1
    key = tuple(lam)
    if key in memo:
        return memo[key]
    total = 0
    for r in range(len(lam)):
        if r + 1 == len(lam) or lam[r + 1] < lam[r]:
            smaller = list(lam)
            smaller[r] -= 1
            if smaller[r] == 0:
                smaller.pop()
            total += count_syt_recurrence(smaller, memo)
    memo[key] = total
    return total


def count_syt_brute_force(lam: Shape) -> int:
    """Definition-level brute force: fill the cells in row-major order with a
    permutation of 1..n so that every row and every column increases."""
    n = sum(lam)
    grid = [[0] * length for length in lam]
    used = [False] * (n + 1)

    def go(r: int, c: int) -> int:
        if r == len(lam):
            return 1
        nr, nc = (r + 1, 0) if c + 1 == lam[r] else (r, c + 1)
        total = 0
        for v in range(1, n + 1):
            if used[v]:
                continue
            if c > 0 and grid[r][c - 1] > v:
                continue
            if r > 0 and grid[r - 1][c] > v:
                continue
            used[v] = True
            grid[r][c] = v
            total += go(nr, nc)
            used[v] = False
        return total

    return go(0, 0)


def count_syt_determinant_mod(lam: Shape) -> int:
    """n! * det[1/(b_j - a_i)!] modulo P, with a_i = m-1-i and b_j = lam_j + m-1-j,
    by Gaussian elimination: the determinant of the first proof, evaluated numerically."""
    m = len(lam)
    n = sum(lam)
    inv_fact = [pow(factorial(k) % P, P - 2, P) for k in range(lam[0] + m + 1)]
    a = [m - 1 - i for i in range(m)]
    b = [lam[j] + m - 1 - j for j in range(m)]
    mat = [[inv_fact[b[j] - a[i]] if b[j] >= a[i] else 0 for j in range(m)] for i in range(m)]
    det = 1
    for col in range(m):
        pivot = next((r for r in range(col, m) if mat[r][col] != 0), None)
        if pivot is None:
            return 0
        if pivot != col:
            mat[pivot], mat[col] = mat[col], mat[pivot]
            det = (P - det) % P
        det = det * mat[col][col] % P
        inv = pow(mat[col][col], P - 2, P)
        for r in range(col + 1, m):
            factor = mat[r][col] * inv % P
            for k in range(col, m):
                mat[r][k] = (mat[r][k] + P - factor * mat[col][k] % P) % P
    return det * (factorial(n) % P) % P


def walk_probabilities(lam: Shape, corner_cell: tuple[int, int]) -> list[list[float]]:
    """P(u, c) of the hook walk, as floating point, for every cell u and the
    corner c = corner_cell: probability that a walk started at u ends at c."""
    conj = conjugate(lam)
    hooks = hook_lengths(lam)
    prob = [[0.0] * length for length in lam]
    for r in range(len(lam) - 1, -1, -1):
        for c in range(lam[r] - 1, -1, -1):
            if hooks[r][c] == 1:
                prob[r][c] = 1.0 if (r, c) == corner_cell else 0.0
                continue
            total = 0.0
            for c2 in range(c + 1, lam[r]):
                total += prob[r][c2]
            for r2 in range(r + 1, conj[c]):
                total += prob[r2][c]
            prob[r][c] = total / (hooks[r][c] - 1)
    return prob


class Lcg:
    def __init__(self, seed: int) -> None:
        self.state = seed

    def below(self, bound: int) -> int:
        self.state = (self.state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        return (self.state >> 33) % bound


def product(values) -> int:
    out = 1
    for v in values:
        out *= v
    return out


def main() -> None:
    # --- A worked shape ---------------------------------------------------
    lam: Shape = [4, 2, 1]
    print(f"hook lengths of (4,2,1): {hook_lengths(lam)}")
    print(f"f(4,2,1) = {count_syt_exact(lam)} (7! / 144)")

    # --- Four independent counts agree on every partition of n <= 12 --------
    memo: dict[tuple[int, ...], int] = {}
    all_ok = True
    shapes_checked = 0
    for n in range(13):
        for part in partitions(n):
            shapes_checked += 1
            rec = count_syt_recurrence(part, memo)
            hook = count_syt_exact(part)
            hook_mod = count_syt_mod(part)
            if not (rec == hook and hook_mod == hook % P):
                all_ok = False
                print(f"MISMATCH (recurrence/hook) {part}: {rec} vs {hook}")
            if part:
                frob = count_syt_frobenius_mod(part)
                det = count_syt_determinant_mod(part)
                if frob != hook_mod or det != hook_mod:
                    all_ok = False
                    print(f"MISMATCH (Frobenius/determinant) {part}: {frob} {det} vs {hook_mod}")
    print(f"recurrence, hook formula, Frobenius product and determinant agree on all {shapes_checked} partitions of n <= 12: {'PASS' if all_ok else 'FAIL'}")

    # --- Definition-level brute force for n <= 8 ---------------------------
    brute_ok = True
    for n in range(9):
        for part in partitions(n):
            if count_syt_brute_force(part) != count_syt_exact(part):
                brute_ok = False
                print(f"BRUTE FORCE MISMATCH {part}")
    print(f"brute-force enumeration of tableaux matches the formula for all partitions of n <= 8: {'PASS' if brute_ok else 'FAIL'}")

    # --- sum over lambda of f^2 = n! (a consequence of the RSK correspondence)
    sq_ok = True
    for n in range(1, 13):
        total = sum(count_syt_exact(part) ** 2 for part in partitions(n))
        if total != factorial(n):
            sq_ok = False
    print(f"sum of f^2 over partitions of n equals n! for n <= 12: {'PASS' if sq_ok else 'FAIL'}")

    # --- Two-row rectangles give the Catalan numbers -------------------------
    catalan_ok = True
    catalan = 1
    for n in range(11):
        if n > 0:
            catalan = catalan * 2 * (2 * n - 1) // (n + 1)
        part = [] if n == 0 else [n, n]
        if count_syt_exact(part) != catalan:
            catalan_ok = False
    print(f"f(n,n) equals the n-th Catalan number for n <= 10: {'PASS' if catalan_ok else 'FAIL'}")

    # --- A large shape, modulo P ---------------------------------------------
    big: Shape = [500, 500]
    catalan_500 = factorial(1000) % P * pow(factorial(501) % P * (factorial(500) % P) % P, P - 2, P) % P
    ok = count_syt_mod(big) == catalan_500 and count_syt_frobenius_mod(big) == catalan_500
    print(f"f(500,500) mod p equals Catalan(500) mod p: {'PASS' if ok else 'FAIL'}")

    # --- The row-hook lemma: hooks of row i together with the gaps b_i - b_j are exactly {1..b_i}
    lemma_ok = True
    for n in range(1, 13):
        for part in partitions(n):
            m = len(part)
            b = [part[i] + m - 1 - i for i in range(m)]
            hooks = hook_lengths(part)
            for i in range(m):
                values = sorted(hooks[i] + [b[i] - b[j] for j in range(i + 1, m)])
                if values != list(range(1, b[i] + 1)):
                    lemma_ok = False
    print(f"row-hook lemma ({{hooks of row i}} + {{b_i - b_j}} = {{1..b_i}}) holds for all partitions of n <= 12: {'PASS' if lemma_ok else 'FAIL'}")

    # --- Hook walk reaches corner c with total probability H(lambda)/H(lambda - c)
    walk_ok = True
    for n in range(1, 11):
        for part in partitions(n):
            h_lam = float(product(h for row in hook_lengths(part) for h in row))
            for r in range(len(part)):
                if r + 1 == len(part) or part[r + 1] < part[r]:
                    corner = (r, part[r] - 1)
                    smaller = list(part)
                    smaller[r] -= 1
                    if smaller[r] == 0:
                        smaller.pop()
                    h_small = float(product(h for row in hook_lengths(smaller) for h in row))
                    total = sum(x for row in walk_probabilities(part, corner) for x in row)
                    if abs(total - h_lam / h_small) > 1e-9 * (1.0 + total):
                        walk_ok = False
    print(f"sum_u P(u, c) = H(lambda)/H(lambda - c) for every corner of every partition of n <= 10: {'PASS' if walk_ok else 'FAIL'}")

    # --- The hook-walk sampler is uniform --------------------------------------
    shape: Shape = [3, 2, 1]
    f = count_syt_exact(shape)
    trials = 16_000
    rng = Lcg(2026)
    counts: dict[tuple[tuple[int, ...], ...], int] = {}
    for _ in range(trials):
        t = tuple(tuple(row) for row in sample_tableau(shape, rng.below))
        counts[t] = counts.get(t, 0) + 1
    expected = trials / f
    sigma = (expected * (1.0 - 1.0 / f)) ** 0.5
    worst = max(abs(c - expected) for c in counts.values())
    ok = len(counts) == f and worst < 4.0 * sigma
    print(f"hook-walk sampler on (3,2,1): {len(counts)} distinct tableaux seen (f = {f}), every count within 4 sigma of uniform: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
