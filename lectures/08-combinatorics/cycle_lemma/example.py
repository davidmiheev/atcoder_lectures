"""The Cycle Lemma: for an integer sequence A of length N with every
A_i <= 1 and nonnegative sum K, exactly K of its N cyclic shifts have
every prefix sum (after the empty prefix) strictly positive. Several
classical counting formulas -- Catalan, Fuss-Catalan, Narayana, and the
Catalan convolution identity -- follow from it via a "prepend +1" trick.
"""

# PDFBEGIN
from __future__ import annotations

import itertools

import math

import random

def catalan(n: int) -> int:
    return math.comb(2 * n, n) // (n + 1)

def prefix_always_nonneg(a: list[int]) -> bool:
    s = 0
    for x in a:
        s += x
        if s < 0:
            return False
    return True

def prefix_positive_after_first(a: list[int]) -> bool:
    s = 0
    for x in a:
        s += x
        if s <= 0:
            return False
    return True

def cycle_lemma_count(a: list[int]) -> int:
    """The Cycle Lemma count: how many of A's N cyclic shifts have every
    prefix sum (excluding S_0) strictly positive.
    """
    n = len(a)
    count = 0
    for j in range(n):
        shifted = a[j:] + a[:j]
        if prefix_positive_after_first(shifted):
            count += 1
    return count

# PDFEND -- everything below is the test harness and is not printed in the PDF

def all_sequences(total: int, num_plus: int, plus_val: int, minus_val: int) -> list[list[int]]:
    """Every length-`total` sequence with `plus_val` at `num_plus`
    positions (chosen every possible way) and `minus_val` everywhere
    else.
    """
    sequences = []
    for combo in itertools.combinations(range(total), num_plus):
        seq = [minus_val] * total
        for idx in combo:
            seq[idx] = plus_val
        sequences.append(seq)
    return sequences

def main() -> None:
    # --- Core theorem, on the lecture's own worked example ---------------
    example = [1, 0, -2, 1, -1, 1, 0, 0, 1, 1]
    k = sum(example)
    got = cycle_lemma_count(example)
    print(f"lecture's own example: sum K = {k}, cyclic shifts with all-positive prefix sums = {got}")
    assert k == 2
    assert got == 2

    # --- Core theorem, on random sequences with A_i <= 1, sum >= 0 -------
    rng = random.Random(99)
    core_ok = True
    for _ in range(500):
        n = rng.randint(1, 12)
        a = [rng.randint(-3, 1) for _ in range(n)]
        s = sum(a)
        if s < 0:
            deficit = -s
            idx = rng.randrange(n)
            a[idx] += deficit
            if a[idx] > 1:
                continue
        k = sum(a)
        got = cycle_lemma_count(a)
        if got != k:
            core_ok = False
            print(f"CYCLE LEMMA MISMATCH a={a}: K={k} got={got}")
    print(f"Cycle Lemma holds on random sequences (A_i <= 1, sum >= 0): {'PASS' if core_ok else 'FAIL'}")

    # --- Catalan numbers: both formulas from Theorem 6 --------------------
    catalan_ok = True
    for n in range(7):
        # Condition 1: correct bracket sequences of length 2n.
        count1 = 1 if n == 0 else sum(1 for a in all_sequences(2 * n, n, 1, -1) if prefix_always_nonneg(a))
        # Condition 2: length 2n+1, n+1 plus-ones, prefix positive after first.
        seqs_b = all_sequences(2 * n + 1, n + 1, 1, -1)
        count2 = sum(1 for a in seqs_b if prefix_positive_after_first(a))
        # Condition 4: same domain, but just those starting with +1.
        count4 = sum(1 for a in seqs_b if a[0] == 1)

        formula = catalan(n)
        total3 = math.comb(2 * n + 1, n)  # condition 3: all such sequences

        cycle_lemma_relation = count2 * (2 * n + 1) == total3
        second_relation = count2 * (n + 1) == count4
        if count1 != formula or count2 != formula or count4 != math.comb(2 * n, n) or not cycle_lemma_relation or not second_relation:
            catalan_ok = False
            print(
                f"CATALAN MISMATCH n={n}: count1={count1} count2={count2} count4={count4} "
                f"formula={formula} cycle_lemma_relation={cycle_lemma_relation} second_relation={second_relation}"
            )
    print(f"Catalan numbers: both formulas C_N = C(2N+1,N)/(2N+1) = C(2N,N)/(N+1) verified via Cycle Lemma: {'PASS' if catalan_ok else 'FAIL'}")

    # --- Fuss-Catalan numbers (Theorem 8) ----------------------------------
    fuss_ok = True
    for n in range(1, 5):
        for m in range(1, 4):
            total = n * m + n
            num_plus = n * m
            count = sum(1 for a in all_sequences(total, num_plus, 1, -m) if prefix_always_nonneg(a))
            formula = math.comb(total, n) // (n * m + 1)
            if count != formula:
                fuss_ok = False
                print(f"FUSS-CATALAN MISMATCH n={n} m={m}: count={count} formula={formula}")
    print(f"Fuss-Catalan numbers (Theorem 8) verified for small N,M: {'PASS' if fuss_ok else 'FAIL'}")

    # --- Narayana numbers (Theorem 9) --------------------------------------
    narayana_ok = True
    for n in range(1, 7):
        for m in range(1, n + 1):
            count = 0
            for a in all_sequences(2 * n, n, 1, -1):
                if not prefix_always_nonneg(a):
                    continue
                peaks = sum(1 for i in range(len(a) - 1) if a[i] == 1 and a[i + 1] == -1)
                if peaks == m:
                    count += 1
            formula = math.comb(n, m) * math.comb(n, m - 1) // n
            if count != formula:
                narayana_ok = False
                print(f"NARAYANA MISMATCH n={n} m={m}: count={count} formula={formula}")
    print(f"Narayana numbers (Theorem 9) verified for n<=6: {'PASS' if narayana_ok else 'FAIL'}")

    # --- Catalan convolution formula (Theorem 7) ---------------------------
    def compositions(total: int, parts: int) -> list[list[int]]:
        if parts == 1:
            return [[total]]
        result = []
        for first in range(total + 1):
            for rest in compositions(total - first, parts - 1):
                result.append([first] + rest)
        return result

    conv_ok = True
    for n in range(6):
        for m in range(1, 5):
            total_sum = sum(math.prod(catalan(x) for x in parts) for parts in compositions(n, m))
            formula = (m * math.comb(2 * n + m, n)) // (2 * n + m)
            if total_sum != formula:
                conv_ok = False
                print(f"CONVOLUTION MISMATCH n={n} m={m}: sum={total_sum} formula={formula}")
    print(f"Catalan convolution formula (Theorem 7) verified for small N,M: {'PASS' if conv_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
