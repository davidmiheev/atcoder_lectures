"""Min-plus convolution of convex sequences: C_k = min_{i+j=k} (A_i + B_j).

In general this takes O(NM) time (and is believed to need essentially
that much). But if A and B are both convex ("convex below": nondecreasing
differences), it collapses to O(N+M): reinterpret A_i as "cheapest way to
pick i items from a bag priced by A's differences," and the convolution
becomes "cheapest way to pick k items from the two bags combined" -- i.e.
just merge the two (already-sorted) difference sequences.
"""
from __future__ import annotations

import random


def min_plus_convolution_naive(a: list[int], b: list[int]) -> list[int]:
    """General min-plus convolution, O(N*M): the definition, directly."""
    n, m = len(a) - 1, len(b) - 1
    c = [float("inf")] * (n + m + 1)
    for i in range(n + 1):
        for j in range(m + 1):
            s = a[i] + b[j]
            if s < c[i + j]:
                c[i + j] = s
    return c


def is_convex(a: list[int]) -> bool:
    """"Convex below": differences A_{i+1} - A_i are nondecreasing."""
    diffs = [a[i + 1] - a[i] for i in range(len(a) - 1)]
    return all(diffs[i] <= diffs[i + 1] for i in range(len(diffs) - 1))


def differences(a: list[int]) -> list[int]:
    return [a[i + 1] - a[i] for i in range(len(a) - 1)]


def merge_sorted(a: list[int], b: list[int]) -> list[int]:
    """Merge two already-sorted (nondecreasing) lists in O(n+m)."""
    result = []
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i] <= b[j]:
            result.append(a[i])
            i += 1
        else:
            result.append(b[j])
            j += 1
    result.extend(a[i:])
    result.extend(b[j:])
    return result


def min_plus_convolution_convex(a: list[int], b: list[int]) -> list[int]:
    """Min-plus convolution of two convex sequences, O(N+M): merge their
    (individually nondecreasing) difference sequences and re-accumulate.
    """
    assert is_convex(a) and is_convex(b)
    merged_diffs = merge_sorted(differences(a), differences(b))
    c = [a[0] + b[0]]
    for d in merged_diffs:
        c.append(c[-1] + d)
    return c


def random_convex_sequence(length: int, rng: random.Random) -> list[int]:
    """A random convex sequence of the given length: nondecreasing random
    differences, accumulated from a random starting value.
    """
    diffs = sorted(rng.randint(-5, 5) for _ in range(max(length - 1, 0)))
    seq = [rng.randint(-10, 10)]
    for d in diffs:
        seq.append(seq[-1] + d)
    return seq


def main() -> None:
    # --- The lecture's own worked examples ---------------------------------
    a1 = [5, 1, 4, 2, 2]
    b1 = [1, 3, 2]
    expected1 = [6, 2, 4, 3, 3, 4, 4]
    got1 = min_plus_convolution_naive(a1, b1)
    print(f"worked example 1 (general): {got1} (expected {expected1})")
    assert got1 == expected1

    a2 = [2, 0, 0, 1, 2]
    b2 = [0, 0, 1, 3]
    expected2 = [2, 0, 0, 0, 1, 2, 3, 5]
    assert is_convex(a2) and is_convex(b2)
    got2_fast = min_plus_convolution_convex(a2, b2)
    got2_naive = min_plus_convolution_naive(a2, b2)
    print(f"worked example 2 (convex):  {got2_fast} (expected {expected2})")
    assert got2_fast == expected2
    assert got2_naive == expected2

    # --- Random convex sequences: fast algorithm vs. naive definition -----
    rng = random.Random(2024)
    convolution_ok = True
    result_convex_ok = True
    for _ in range(300):
        n = rng.randint(1, 15)
        m = rng.randint(1, 15)
        a = random_convex_sequence(n, rng)
        b = random_convex_sequence(m, rng)
        fast = min_plus_convolution_convex(a, b)
        naive = min_plus_convolution_naive(a, b)
        if fast != naive:
            convolution_ok = False
            print(f"CONVOLUTION MISMATCH a={a} b={b}: fast={fast} naive={naive}")
        if not is_convex(fast):
            result_convex_ok = False
            print(f"RESULT NOT CONVEX a={a} b={b}: c={fast}")
    print(f"fast convex-convolution matches the naive O(NM) definition: {'PASS' if convolution_ok else 'FAIL'}")
    print(f"the convolution of two convex sequences is itself convex: {'PASS' if result_convex_ok else 'FAIL'}")

    # --- Problem 4: "choose exactly k elements, minimize the sum" ----------
    # The textbook O(N^2) DP is repeated min-plus convolution with (0, x).
    # Understanding that the dp table stays convex, and that convolving
    # with (0, x) just inserts x into the sorted difference multiset,
    # shows this dp is secretly just "sort A and take prefix sums."
    def slow_dp(a: list[int]) -> list[int]:
        dp = [0]
        for x in a:
            newdp = [float("inf")] * (len(dp) + 1)
            for i, v in enumerate(dp):
                newdp[i] = min(newdp[i], v)  # don't use x
                newdp[i + 1] = min(newdp[i + 1], v + x)  # use x
            dp = newdp
        return dp

    def fast_via_sort(a: list[int]) -> list[int]:
        dp = [0]
        for x in sorted(a):
            dp.append(dp[-1] + x)
        return dp

    dp_ok = True
    dp_stays_convex = True
    for _ in range(100):
        n = rng.randint(1, 10)
        a = [rng.randint(-5, 10) for _ in range(n)]
        slow = slow_dp(a)
        fast = fast_via_sort(a)
        if slow != fast:
            dp_ok = False
            print(f"DP MISMATCH a={a}: slow_dp={slow} fast_via_sort={fast}")
        # Confirm the dp table is convex after every prefix of A, as the
        # theory predicts (each step is a convex-convex convolution).
        dp = [0]
        for x in a:
            newdp = [float("inf")] * (len(dp) + 1)
            for i, v in enumerate(dp):
                newdp[i] = min(newdp[i], v)
                newdp[i + 1] = min(newdp[i + 1], v + x)
            dp = newdp
            if not is_convex(dp):
                dp_stays_convex = False
                print(f"DP NOT CONVEX a={a} partial dp={dp}")
    print(f"naive O(N^2) knapsack-min dp matches sort-and-prefix-sum: {'PASS' if dp_ok else 'FAIL'}")
    print(f"the dp table stays convex at every step (as the theorem predicts): {'PASS' if dp_stays_convex else 'FAIL'}")


if __name__ == "__main__":
    main()
