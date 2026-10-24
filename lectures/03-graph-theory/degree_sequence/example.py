"""Degree sequences of simple undirected graphs: the Erdos-Gallai
realizability test and the Havel-Hakimi construction algorithm.
"""

# PDFBEGIN
from __future__ import annotations

import bisect

import itertools

def is_graphic_erdos_gallai(d: list[int]) -> bool:
    """Is `d` (any order) the degree sequence of some simple undirected
    graph? Checks the handshake-lemma parity condition and, for every
    prefix length k of the sorted-descending sequence, the Erdos-Gallai
    inequality sum(d[:k]) <= k(k-1) + sum_{i>=k} min(d[i], k).

    The right-hand sum is evaluated in O(log n) per k by binary-searching
    the (sorted, so monotone) boundary between "d[i] >= k" and "d[i] < k"
    in the suffix, then combining with prefix sums -- O(n log n) overall.
    """
    n = len(d)
    if any(x < 0 or x > n - 1 for x in d):
        return False
    if sum(d) % 2 != 0:
        return False

    sorted_desc = sorted(d, reverse=True)
    prefix = [0] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] + sorted_desc[i]

    for k in range(1, n + 1):
        lhs = prefix[k]
        suffix = sorted_desc[k:]  # already descending
        # bisect on the negated suffix to find how many entries are >= k
        neg_suffix = [-v for v in suffix]
        split = bisect.bisect_left(neg_suffix, -k)
        sum_ge = split * k
        sum_lt = (prefix[n] - prefix[k]) - sum(suffix[:split])
        rhs = k * (k - 1) + sum_ge + sum_lt
        if lhs > rhs:
            return False
    return True

def havel_hakimi(d: list[int]) -> list[tuple[int, int]] | None:
    """Construct a simple graph with degree sequence `d` (any order), or
    return None if none exists. Returns 0-indexed edges referring to the
    ORIGINAL indices of `d`.

    At each round, connects the current highest-degree vertex to the next
    `d1` highest-degree remaining vertices, then "retires" it.
    """
    n = len(d)
    if not is_graphic_erdos_gallai(d):
        return None
    remaining = list(zip(d, range(n)))
    edges: list[tuple[int, int]] = []

    while True:
        remaining.sort(key=lambda p: -p[0])
        d1, v1 = remaining[0]
        if d1 == 0:
            break
        if d1 >= len(remaining):
            return None
        for i in range(1, d1 + 1):
            deg, idx = remaining[i]
            if deg <= 0:
                return None
            edges.append((v1, idx))
            remaining[i] = (deg - 1, idx)
        remaining[0] = (0, v1)
    return edges

# PDFEND -- everything below is the test harness and is not printed in the PDF

def verify_degrees(n: int, edges: list[tuple[int, int]], want: list[int]) -> bool:
    got = [0] * n
    for a, b in edges:
        got[a] += 1
        got[b] += 1
    return got == want

def main() -> None:
    # The three worked examples from the lecture text.
    cases = [
        ("(4,4,3,2,2,1)", [4, 4, 3, 2, 2, 1]),
        ("(4,3,3,2,2,1)", [4, 3, 3, 2, 2, 1]),
        ("(4,4,4,2,1,1)", [4, 4, 4, 2, 1, 1]),
    ]

    for label, d in cases:
        graphic = is_graphic_erdos_gallai(d)
        print(f"{label}: graphic = {graphic}")
        edges = havel_hakimi(d)
        if edges is not None:
            ok = verify_degrees(len(d), edges, d)
            print(f"  constructed {len(edges)} edges: {edges}")
            print(f"  degree check: {'PASS' if ok else 'FAIL'}")
        else:
            print("  no construction (correctly not graphic)")

    # Brute-force cross-check: every possible 5-vertex graph's degree
    # sequence should be accepted by is_graphic_erdos_gallai.
    n = 5
    all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    mismatches = 0
    for mask in range(1 << len(all_pairs)):
        edges = [e for bit, e in enumerate(all_pairs) if mask & (1 << bit)]
        deg = [0] * n
        for a, b in edges:
            deg[a] += 1
            deg[b] += 1
        sorted_deg = sorted(deg, reverse=True)
        if not is_graphic_erdos_gallai(sorted_deg):
            mismatches += 1
            print(f"MISMATCH: degree sequence {sorted_deg} came from a real graph but is_graphic said false")

    print(
        f"\nexhaustive check over all {n}-vertex graphs: every realized degree "
        f"sequence passes is_graphic ({mismatches} mismatches)"
    )

if __name__ == "__main__":
    main()
