"""The transposition principle (Tellegen's principle): an algorithm for
"y = Ax" built out of (1) loading x into some scratch variables, (2) a
sequence of elementary steps z[a] += c * z[b], and (3) reading y out of
some scratch variables, can be mechanically transformed into an algorithm
for the TRANSPOSED problem "X = A^T Y" by reversing the step order,
swapping each step's a/b roles, and swapping the roles of the load/store
positions.

This file implements a generic interpreter for such programs (forward and
transposed), and: (1) validates the transposition principle itself
against the matrix identity it is a restatement of, over many random
programs; (2) applies it to the classic prefix-sum solution of range-sum
queries to mechanically derive a range-add solution, cross-checked
against both a hand-derived variant from the lecture and the well-known
"imos method"; (3) applies it to the fast subset zeta transform to
mechanically derive the fast superset zeta transform; and (4) reproduces
the lecture's own small worked example exactly.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field


@dataclass
class Program:
    """A "linear program": load x into some scratch slots, run a sequence
    of z[a] += c*z[b] steps, then read y out of some scratch slots."""

    k: int
    load_positions: list[int]  # length M: x[j] goes into z[load_positions[j]]
    ops: list[tuple[int, int, int]] = field(default_factory=list)  # (a, b, c): z[a] += c*z[b]
    store_positions: list[int] = field(default_factory=list)  # length N: y[i] = z[store_positions[i]]


def forward(prog: Program, x: list[int]) -> list[int]:
    z = [0] * prog.k
    for j, pos in enumerate(prog.load_positions):
        # Additive, not overwriting: if two x_j share a load position, the
        # matrix-multiplication reading of "z = Sx" (S a 0/1 selection
        # matrix) sums their contributions rather than keeping only the
        # last one. The lecture's own examples never hit this case (every
        # load position it uses is distinct), but summing is what makes
        # the transpose identity hold in general, and it costs nothing
        # when positions are already distinct.
        z[pos] += x[j]
    for a, b, c in prog.ops:
        z[a] += c * z[b]
    return [z[p] for p in prog.store_positions]


def transpose_run(prog: Program, y: list[int]) -> list[int]:
    """The transposed algorithm from Theorem 2: reverse step order, swap
    each step's (a,b) roles, and swap the roles of load/store positions."""
    z = [0] * prog.k
    for i, pos in enumerate(prog.store_positions):
        # Additive for the same reason as in `forward`'s load step, mirrored
        # here since the transposed algorithm loads Y through T^T.
        z[pos] += y[i]
    for a, b, c in reversed(prog.ops):
        z[b] += c * z[a]
    return [z[p] for p in prog.load_positions]


def compute_matrix(prog: Program) -> list[list[int]]:
    """Computes the N x M matrix A that `forward` implements, by running it
    on each standard basis vector. Used only to validate the transposition
    principle itself against A^T, not part of either algorithm."""
    m = len(prog.load_positions)
    n = len(prog.store_positions)
    a = [[0] * m for _ in range(n)]
    for j in range(m):
        x = [0] * m
        x[j] = 1
        y = forward(prog, x)
        for i in range(n):
            a[i][j] = y[i]
    return a


def mat_transpose_vec(a: list[list[int]], y: list[int]) -> list[int]:
    n = len(a)
    m = len(a[0]) if n > 0 else 0
    return [sum(a[i][j] * y[i] for i in range(n)) for j in range(m)]


def random_program(rng: random.Random) -> Program:
    k = rng.randint(2, 7)
    m = rng.randint(1, k)
    n = rng.randint(1, k)
    q = rng.randint(0, 14)
    load_positions = [rng.randint(0, k - 1) for _ in range(m)]
    store_positions = [rng.randint(0, k - 1) for _ in range(n)]
    ops = [(rng.randint(0, k - 1), rng.randint(0, k - 1), rng.randint(-5, 5)) for _ in range(q)]
    return Program(k=k, load_positions=load_positions, ops=ops, store_positions=store_positions)


# --- Section 5: range-sum queries (Problem 3A) and their transpose,
# range-add queries (Problem 3B), via the prefix-sum solution. ---


def interval_sum_program(n: int, queries: list[tuple[int, int]]) -> Program:
    """Builds a Program implementing Problem 3(A) exactly as in the
    lecture: scratch variables s_0..s_n (indices 0..n) and y_1..y_q
    (indices n+1..n+q), computing y_i = sum_{j=L_i}^{R_i} x_j via prefix
    sums."""
    q = len(queries)
    k = (n + 1) + q
    load_positions = [j + 1 for j in range(n)]
    ops: list[tuple[int, int, int]] = []
    for j in range(n):
        ops.append((j + 1, j, 1))
    for i, (l, r) in enumerate(queries):
        yi = n + 1 + i
        ops.append((yi, r + 1, 1))
        ops.append((yi, l, -1))
    store_positions = [n + 1 + i for i in range(q)]
    return Program(k=k, load_positions=load_positions, ops=ops, store_positions=store_positions)


def interval_sum_naive(queries: list[tuple[int, int]], x: list[int]) -> list[int]:
    return [sum(x[l : r + 1]) for l, r in queries]


def interval_add_naive(n: int, queries: list[tuple[int, int]], y: list[int]) -> list[int]:
    return [sum(yi for (l, r), yi in zip(queries, y) if l <= j <= r) for j in range(n)]


def transpose_hand_derived(n: int, queries: list[tuple[int, int]], y: list[int]) -> list[int]:
    """The hand-derived transpose from section 5.3 of the lecture."""
    s = [0] * (n + 1)
    for i, (l, r) in enumerate(queries):
        s[l] -= y[i]
        s[r + 1] += y[i]
    for j in range(n - 1, -1, -1):
        s[j] += s[j + 1]
    return [s[j + 1] for j in range(n)]


def transpose_imos(n: int, queries: list[tuple[int, int]], y: list[int]) -> list[int]:
    """The "imos method" variant from section 5.4 (different loop
    direction and sign convention, same result)."""
    s = [0] * (n + 1)
    for i, (l, r) in enumerate(queries):
        s[l] += y[i]
        s[r + 1] -= y[i]
    for j in range(n):
        s[j + 1] += s[j]
    return [s[j] for j in range(n)]


# --- Section 6: subset zeta transform (Problem 4A) and its transpose,
# the superset zeta transform (Problem 4B). ---


def subset_zeta_program(n: int) -> Program:
    size = 1 << n
    load_positions = list(range(size))
    ops: list[tuple[int, int, int]] = []
    for i in range(n):
        for t in range(size):
            if (t >> i) & 1:
                ops.append((t, t ^ (1 << i), 1))
    store_positions = list(range(size))
    return Program(k=size, load_positions=load_positions, ops=ops, store_positions=store_positions)


def subset_zeta_naive(n: int, a: list[int]) -> list[int]:
    size = 1 << n
    return [sum(a[s] for s in range(size) if s & t == s) for t in range(size)]


def superset_zeta_naive(n: int, a: list[int]) -> list[int]:
    size = 1 << n
    return [sum(a[s] for s in range(size) if s & t == t) for t in range(size)]


def main() -> None:
    rng = random.Random(2026092001)

    # --- Theorem 2 itself: transpose_run(Y) == A^T Y, for many random programs ---
    theorem_ok = True
    for _ in range(300):
        prog = random_program(rng)
        n = len(prog.store_positions)
        a = compute_matrix(prog)
        y = [rng.randint(-10, 10) for _ in range(n)]
        via_transpose = transpose_run(prog, y)
        via_matrix = mat_transpose_vec(a, y)
        if via_transpose != via_matrix:
            theorem_ok = False
            print("TRANSPOSITION PRINCIPLE MISMATCH")
    print(f"Transposition principle (transpose_run(Y) == A^T Y) over 300 random programs: {'PASS' if theorem_ok else 'FAIL'}")

    # --- Worked example (section 3.3): y1=2x1+3x2, y2=4x1+5x2, y3=6x1+7x2 ---
    example = Program(
        k=5,
        load_positions=[0, 1],
        ops=[(2, 0, 2), (2, 1, 3), (3, 0, 4), (3, 1, 5), (4, 0, 6), (4, 1, 7)],
        store_positions=[2, 3, 4],
    )
    example_ok = True
    for _ in range(50):
        x1 = rng.randint(-20, 20)
        x2 = rng.randint(-20, 20)
        y = forward(example, [x1, x2])
        expected_y = [2 * x1 + 3 * x2, 4 * x1 + 5 * x2, 6 * x1 + 7 * x2]
        y1 = rng.randint(-20, 20)
        y2 = rng.randint(-20, 20)
        y3 = rng.randint(-20, 20)
        x_t = transpose_run(example, [y1, y2, y3])
        expected_x = [2 * y1 + 4 * y2 + 6 * y3, 3 * y1 + 5 * y2 + 7 * y3]
        if y != expected_y or x_t != expected_x:
            example_ok = False
            print("WORKED EXAMPLE MISMATCH")
    print(f"Worked example (section 3.3): forward and transpose match the lecture's own formulas, over 50 random trials: {'PASS' if example_ok else 'FAIL'}")

    # --- Section 5: range-sum queries and their transpose ---
    interval_ok = True
    for _ in range(100):
        n = rng.randint(1, 20)
        q = rng.randint(0, 14)
        queries = []
        for _ in range(q):
            l = rng.randint(0, n - 1)
            r = rng.randint(l, n - 1)
            queries.append((l, r))
        prog = interval_sum_program(n, queries)

        x = [rng.randint(-20, 20) for _ in range(n)]
        y_via_prog = forward(prog, x)
        y_via_naive = interval_sum_naive(queries, x)
        if y_via_prog != y_via_naive:
            interval_ok = False
            print(f"INTERVAL SUM MISMATCH n={n} q={q}")

        y = [rng.randint(-20, 20) for _ in range(q)]
        x_via_transpose = transpose_run(prog, y)
        x_via_naive = interval_add_naive(n, queries, y)
        x_via_hand = transpose_hand_derived(n, queries, y)
        x_via_imos = transpose_imos(n, queries, y)
        if x_via_transpose != x_via_naive or x_via_transpose != x_via_hand or x_via_transpose != x_via_imos:
            interval_ok = False
            print(f"INTERVAL ADD MISMATCH n={n} q={q}")
    print(f"Range-sum queries (Problem 3A) match brute force, over 100 random trials: {'PASS' if interval_ok else 'FAIL'}")
    print(f"Mechanically-derived transpose == hand-derived (5.3) == imos method (5.4) == brute-force range-add, over 100 random trials: {'PASS' if interval_ok else 'FAIL'}")

    # --- Section 6: subset zeta transform and its transpose ---
    zeta_ok = True
    for _ in range(30):
        n = rng.randint(0, 6)
        size = 1 << n
        prog = subset_zeta_program(n)

        x = [rng.randint(-10, 10) for _ in range(size)]
        y_via_prog = forward(prog, x)
        y_via_naive = subset_zeta_naive(n, x)
        if y_via_prog != y_via_naive:
            zeta_ok = False
            print(f"SUBSET ZETA MISMATCH n={n}")

        y = [rng.randint(-10, 10) for _ in range(size)]
        x_via_transpose = transpose_run(prog, y)
        x_via_naive = superset_zeta_naive(n, y)
        if x_via_transpose != x_via_naive:
            zeta_ok = False
            print(f"SUPERSET ZETA (TRANSPOSE) MISMATCH n={n}")
    print(f"Subset zeta transform (Problem 4A) and its mechanically-derived transpose, the superset zeta transform (Problem 4B), over 30 random trials: {'PASS' if zeta_ok else 'FAIL'}")

    all_ok = theorem_ok and example_ok and interval_ok and zeta_ok
    print(f"\nAll checks: {'PASS' if all_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
