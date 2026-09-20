"""Zeta and Mobius transforms with respect to the subset/superset relation
on bitmasks. The subset zeta transform B_T = sum_{S subseteq T} A_S (a.k.a.
"sum over subsets", SOS) and its inverse, the subset Mobius transform, are
computed here three ways -- an O(4^n) definition-level double loop, an
O(3^n) submask-enumeration loop, and the O(n*2^n) fast zeta/Mobius
transform that flips one bit's role (fixed value -> wildcard) at a time --
and cross-checked against each other. The superset versions are analogous
with the inclusion direction reversed. Finally, OR- and AND-convolution
are built on top: C_U = sum_{S|T=U} A_S*B_T reduces to a pointwise product
after a subset zeta transform (since S|T=U implies, and nearly is implied
by, S,T subseteq U) followed by a subset Mobius transform back;
AND-convolution is the superset-relation mirror image. All of these are
cross-checked against brute-force O(4^n) definitions.
"""
from __future__ import annotations

import random


def subset_zeta_naive(n: int, a: list[int]) -> list[int]:
    """O(4^n): definition-level double loop."""
    size = 1 << n
    b = [0] * size
    for t in range(size):
        for s in range(size):
            if s & t == s:
                b[t] += a[s]
    return b


def subset_zeta_submask_enum(n: int, a: list[int]) -> list[int]:
    """O(3^n): enumerate submasks of T directly."""
    size = 1 << n
    b = [0] * size
    for t in range(size):
        s = t
        while True:
            b[t] += a[s]
            if s == 0:
                break
            s = (s - 1) & t
    return b


def subset_zeta_fast(n: int, a: list[int]) -> list[int]:
    """O(n*2^n): fast zeta transform (FZT)."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if (t >> i) & 1:
                a[t] += a[t ^ (1 << i)]
    return a


def subset_mobius_fast(n: int, a: list[int]) -> list[int]:
    """O(n*2^n): fast Mobius transform (FMT), the exact inverse of the fast
    zeta transform above."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if (t >> i) & 1:
                a[t] -= a[t ^ (1 << i)]
    return a


def superset_zeta_naive(n: int, a: list[int]) -> list[int]:
    """O(4^n): definition-level double loop for the superset zeta
    transform, B_T = sum_{S superseteq T} A_S."""
    size = 1 << n
    b = [0] * size
    for t in range(size):
        for s in range(size):
            if s & t == t:
                b[t] += a[s]
    return b


def superset_zeta_fast(n: int, a: list[int]) -> list[int]:
    """O(n*2^n): fast superset zeta transform."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if not (t >> i) & 1:
                a[t] += a[t ^ (1 << i)]
    return a


def superset_mobius_fast(n: int, a: list[int]) -> list[int]:
    """O(n*2^n): fast superset Mobius transform."""
    a = a[:]
    size = 1 << n
    for i in range(n):
        for t in range(size):
            if not (t >> i) & 1:
                a[t] -= a[t ^ (1 << i)]
    return a


def or_convolution(n: int, a: list[int], b: list[int]) -> list[int]:
    """C_U = sum_{S|T=U} A_S*B_T, via subset zeta -> pointwise product ->
    subset Mobius."""
    za = subset_zeta_fast(n, a)
    zb = subset_zeta_fast(n, b)
    d = [x * y for x, y in zip(za, zb)]
    return subset_mobius_fast(n, d)


def and_convolution(n: int, a: list[int], b: list[int]) -> list[int]:
    """C_U = sum_{S&T=U} A_S*B_T, via superset zeta -> pointwise product ->
    superset Mobius."""
    za = superset_zeta_fast(n, a)
    zb = superset_zeta_fast(n, b)
    d = [x * y for x, y in zip(za, zb)]
    return superset_mobius_fast(n, d)


def or_convolution_naive(n: int, a: list[int], b: list[int]) -> list[int]:
    size = 1 << n
    c = [0] * size
    for s in range(size):
        for t in range(size):
            c[s | t] += a[s] * b[t]
    return c


def and_convolution_naive(n: int, a: list[int], b: list[int]) -> list[int]:
    size = 1 << n
    c = [0] * size
    for s in range(size):
        for t in range(size):
            c[s & t] += a[s] * b[t]
    return c


def main() -> None:
    rng = random.Random(2026091801)

    # --- Worked example from the lecture: A = (3,1,4,1,5,9,2,6), n=3 ---
    a_example = [3, 1, 4, 1, 5, 9, 2, 6]
    expected_b = [3, 4, 7, 9, 8, 18, 14, 31]
    b_via_naive = subset_zeta_naive(3, a_example)
    b_via_submask = subset_zeta_submask_enum(3, a_example)
    b_via_fast = subset_zeta_fast(3, a_example)
    example_ok = b_via_naive == expected_b and b_via_submask == expected_b and b_via_fast == expected_b
    print(f"Worked example (n=3): all three zeta transform implementations match the lecture's own B: {'PASS' if example_ok else 'FAIL'}")

    # --- Round-trip: subset Mobius(subset zeta(A)) == A ---
    roundtrip_ok = True
    cross_impl_ok = True
    for _ in range(100):
        n = rng.randint(0, 8)
        size = 1 << n
        a = [rng.randint(-50, 50) for _ in range(size)]

        b_naive = subset_zeta_naive(n, a)
        b_submask = subset_zeta_submask_enum(n, a)
        b_fast = subset_zeta_fast(n, a)
        if b_naive != b_submask or b_naive != b_fast:
            cross_impl_ok = False
            print(f"SUBSET ZETA CROSS-IMPL MISMATCH n={n}")

        a_back = subset_mobius_fast(n, b_fast)
        if a_back != a:
            roundtrip_ok = False
            print(f"SUBSET ROUND-TRIP MISMATCH n={n}")

        sup_naive = superset_zeta_naive(n, a)
        sup_fast = superset_zeta_fast(n, a)
        if sup_naive != sup_fast:
            cross_impl_ok = False
            print(f"SUPERSET ZETA CROSS-IMPL MISMATCH n={n}")
        sup_back = superset_mobius_fast(n, sup_fast)
        if sup_back != a:
            roundtrip_ok = False
            print(f"SUPERSET ROUND-TRIP MISMATCH n={n}")

    print(f"Subset/superset zeta transform: naive (O(4^n)), submask-enum (O(3^n)), and fast (O(n*2^n)) all agree, over 100 random trials: {'PASS' if cross_impl_ok else 'FAIL'}")
    print(f"Mobius transform exactly inverts the zeta transform (subset and superset), over 100 random trials: {'PASS' if roundtrip_ok else 'FAIL'}")

    # --- OR- and AND-convolution vs. brute force ---
    or_ok = True
    and_ok = True
    for _ in range(100):
        n = rng.randint(0, 7)
        size = 1 << n
        a = [rng.randint(-10, 10) for _ in range(size)]
        b = [rng.randint(-10, 10) for _ in range(size)]

        or_fast = or_convolution(n, a, b)
        or_naive = or_convolution_naive(n, a, b)
        if or_fast != or_naive:
            or_ok = False
            print(f"OR-CONVOLUTION MISMATCH n={n}")

        and_fast = and_convolution(n, a, b)
        and_naive = and_convolution_naive(n, a, b)
        if and_fast != and_naive:
            and_ok = False
            print(f"AND-CONVOLUTION MISMATCH n={n}")

    print(f"OR-convolution (zeta -> pointwise product -> Mobius) matches brute force, over 100 random trials: {'PASS' if or_ok else 'FAIL'}")
    print(f"AND-convolution (superset zeta -> pointwise product -> superset Mobius) matches brute force, over 100 random trials: {'PASS' if and_ok else 'FAIL'}")

    all_ok = example_ok and cross_impl_ok and roundtrip_ok and or_ok and and_ok
    print(f"\nAll checks: {'PASS' if all_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
