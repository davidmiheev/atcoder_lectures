"""Parking functions: a length-n sequence a_1..a_n with values in [1,n] is a
parking function if, sending car i to look for the first free spot at or
after a_i (spots numbered 1..n), every car finds a spot. There are exactly
(n+1)^(n-1) of them, proved here via a circular extension trick.
"""
from __future__ import annotations

import itertools


def is_parking_function_sim(a: list[int], n: int) -> bool:
    """Definition 1, directly: simulate the parking procedure."""
    occupied = [False] * (n + 1)  # 1-indexed spots 1..n
    for ai in a:
        parked = False
        for j in range(ai, n + 1):
            if not occupied[j]:
                occupied[j] = True
                parked = True
                break
        if not parked:
            return False
    return True


def is_parking_function_condition(a: list[int], n: int) -> bool:
    """Proposition 3, condition 2: for every k, at least k of the a_i are
    <= k."""
    return all(sum(1 for ai in a if ai <= k) >= k for k in range(1, n + 1))


def circular_park_empty_slot(a: list[int], n: int) -> int:
    """The circular extension from the counting proof: n cars, spots
    0..n (n+1 of them), each a_i in 0..n, searching circularly from a_i.
    This never fails; returns the one spot left empty.
    """
    m = n + 1
    occupied = [False] * m
    for ai in a:
        j = ai
        while occupied[j]:
            j = (j + 1) % m
        occupied[j] = True
    return occupied.index(False)


def satisfies_generalized_condition(a: list[int], n: int, m: int) -> bool:
    """Proposition 5's generalized condition, for length-m sequences with
    values in [1,n]: for every k, at least (m-n+k) of the a_i are <= k.
    """
    return all(sum(1 for ai in a if ai <= k) >= (m - n + k) for k in range(1, n + 1))


def main() -> None:
    # --- The source's own worked examples (Problem 2) --------------------
    cases = [
        ([4, 1, 5, 1, 3], True),
        ([5, 3, 1, 4, 3], False),
        ([4, 1, 6, 5, 4, 1], False),
        ([3, 6, 1, 5, 3, 1], True),
    ]
    worked_ok = True
    for a, expected in cases:
        n = len(a)
        sim = is_parking_function_sim(a, n)
        cond = is_parking_function_condition(a, n)
        if sim != expected or cond != expected:
            worked_ok = False
            print(f"WORKED-EXAMPLE MISMATCH a={a}: sim={sim} cond={cond} expected={expected}")
    print(f"lecture's own worked examples (Problem 2) match: {'PASS' if worked_ok else 'FAIL'}")

    # --- Equivalence of the simulation and the threshold condition -------
    equiv_ok = True
    count_ok = True
    for n in range(1, 7):
        count = 0
        for a in itertools.product(range(1, n + 1), repeat=n):
            a = list(a)
            sim = is_parking_function_sim(a, n)
            cond = is_parking_function_condition(a, n)
            if sim != cond:
                equiv_ok = False
                print(f"EQUIV MISMATCH n={n} a={a}: sim={sim} cond={cond}")
            if sim:
                count += 1
        expected = (n + 1) ** (n - 1)
        if count != expected:
            count_ok = False
            print(f"COUNT MISMATCH n={n}: got {count} expected {expected}")
    print(f"simulation matches threshold condition, for n=1..=6: {'PASS' if equiv_ok else 'FAIL'}")
    print(f"parking function count matches (n+1)^(n-1), for n=1..=6: {'PASS' if count_ok else 'FAIL'}")

    # --- Circular extension: bijection + symmetry across slots -----------
    bijection_ok = True
    for n in range(1, 7):
        for a in itertools.product(range(1, n + 1), repeat=n):
            a = list(a)
            leaves_zero_empty = circular_park_empty_slot(a, n) == 0
            is_pf = is_parking_function_sim(a, n)
            if leaves_zero_empty != is_pf:
                bijection_ok = False
                print(f"BIJECTION MISMATCH n={n} a={a}: circular={leaves_zero_empty} linear={is_pf}")
    print(f"circular parking leaves spot 0 empty <=> a is an ordinary parking function, for n=1..=6: {'PASS' if bijection_ok else 'FAIL'}")

    symmetry_ok = True
    for n in range(1, 6):
        tally = [0] * (n + 1)
        for a in itertools.product(range(0, n + 1), repeat=n):
            tally[circular_park_empty_slot(list(a), n)] += 1
        expected = (n + 1) ** (n - 1)
        if any(c != expected for c in tally):
            symmetry_ok = False
            print(f"SYMMETRY MISMATCH n={n}: tally={tally} expected each={expected}")
    print(f"every spot is left empty equally often ((n+1)^(n-1) times), for n=1..=5: {'PASS' if symmetry_ok else 'FAIL'}")

    # --- Proposition 5: the generalized (length-m, values in [1,n]) count -
    prop5_ok = True
    for n in range(1, 7):
        for m in range(1, n + 1):
            count = 0
            for a in itertools.product(range(1, n + 1), repeat=m):
                if satisfies_generalized_condition(list(a), n, m):
                    count += 1
            expected = (n + 1 - m) * (n + 1) ** (m - 1)
            if count != expected:
                prop5_ok = False
                print(f"PROP5 MISMATCH n={n} m={m}: got {count} expected {expected}")
    print(f"Proposition 5 count (n+1-m)(n+1)^(m-1) matches brute force, for n,m<=6: {'PASS' if prop5_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
