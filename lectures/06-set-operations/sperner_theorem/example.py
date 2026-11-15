"""Sperner's theorem: the largest antichain in the Boolean lattice B_n (all
subsets of {0,...,n-1}, ordered by inclusion) has exactly C(n, floor(n/2))
elements. Demonstrates both proofs from the lecture -- the symmetric chain
decomposition (constructive) and the LYM inequality (probabilistic) -- plus
an independent cross-check via Dilworth's theorem (max antichain = min
chain cover, computed by bipartite matching).
"""

# PDFBEGIN
from __future__ import annotations

import math

import random

from fractions import Fraction

def popcount(x: int) -> int:
    return bin(x).count("1")

# ---------------------------------------------------------------------
# Proof 2: symmetric chain decomposition
# ---------------------------------------------------------------------


def symmetric_chain_decomposition(n: int) -> list[list[int]]:
    """Partition B_n into disjoint symmetric chains, via the induction in
    the lecture's proof: extend each chain of B_{n-1} by adding bit (n-1)
    at the top, and peel off a length-1 leftover chain from the middle
    when the parent chain had length > 1.
    """
    if n == 0:
        return [[0]]
    prev = symmetric_chain_decomposition(n - 1)
    bit = 1 << (n - 1)
    result = []
    for chain in prev:
        m = len(chain)
        if m == 1:
            result.append([chain[0], chain[0] | bit])
        else:
            chain1 = chain + [chain[-1] | bit]
            chain2 = [s | bit for s in chain[: m - 1]]
            result.append(chain1)
            result.append(chain2)
    return result

# PDFEND -- everything below is the test harness and is not printed in the PDF

def verify_symmetric_chain_decomposition(n: int, chains: list[list[int]]) -> bool:
    total = 1 << n
    seen = [False] * total
    for chain in chains:
        for a, b in zip(chain, chain[1:]):
            if (a & b) != a or popcount(b) != popcount(a) + 1:
                return False
        lo, hi = popcount(chain[0]), popcount(chain[-1])
        if lo + hi != n:
            return False
        for s in chain:
            if seen[s]:
                return False
            seen[s] = True
    return all(seen)

def level_matching(chains: list[list[int]], k: int) -> list[tuple[int, int]]:
    """Corollary: a maximum matching between B_{n,k} and B_{n,k+1} (each
    edge a proper-subset pair), read directly off the chains that pass
    through both levels.
    """
    edges = []
    for chain in chains:
        for a, b in zip(chain, chain[1:]):
            if popcount(a) == k:
                edges.append((a, b))
    return edges

# ---------------------------------------------------------------------
# Proof 1: LYM inequality (exact rational arithmetic via Fraction)
# ---------------------------------------------------------------------


def lym_sum(n: int, antichain: list[int]) -> Fraction:
    """sum_{S in antichain} 1 / C(n, |S|)."""
    total = Fraction(0)
    for s in antichain:
        k = popcount(s)
        total += Fraction(1, math.comb(n, k))
    return total

def random_maximal_antichain(n: int, rng: random.Random) -> list[int]:
    """A random MAXIMAL (not necessarily maximum) antichain, built
    greedily: shuffle all subsets, add each one that's incomparable to
    everything already chosen.
    """
    order = list(range(1 << n))
    rng.shuffle(order)
    antichain: list[int] = []
    for s in order:
        if not any((s & t == s) or (s & t == t) for t in antichain):
            antichain.append(s)
    return antichain

# ---------------------------------------------------------------------
# Independent cross-check: Dilworth's theorem via bipartite matching.
# max antichain size = min chain cover = (#elements) - (max matching in
# the "u properly-contained-in v" bipartite graph).
# ---------------------------------------------------------------------


def max_bipartite_matching(adj: list[list[int]], n_right: int) -> int:
    match_right = [-1] * n_right

    def try_kuhn(u: int, visited: list[bool]) -> bool:
        for v in adj[u]:
            if visited[v]:
                continue
            visited[v] = True
            if match_right[v] == -1 or try_kuhn(match_right[v], visited):
                match_right[v] = u
                return True
        return False

    result = 0
    for u in range(len(adj)):
        visited = [False] * n_right
        if try_kuhn(u, visited):
            result += 1
    return result

def max_antichain_via_dilworth(n: int) -> int:
    total = 1 << n
    adj: list[list[int]] = [[] for _ in range(total)]
    for u in range(total):
        for v in range(total):
            if u != v and (u & v) == u:
                adj[u].append(v)
    matching = max_bipartite_matching(adj, total)
    return total - matching

def main() -> None:
    # --- Symmetric chain decomposition: validity + count -----------------
    decomp_ok = True
    for n in range(11):
        chains = symmetric_chain_decomposition(n)
        valid = verify_symmetric_chain_decomposition(n, chains)
        expected = math.comb(n, n // 2)
        if not valid or len(chains) != expected:
            decomp_ok = False
            print(f"DECOMP FAIL n={n}: valid={valid} chains={len(chains)} expected={expected}")
    print(f"symmetric chain decomposition valid, with C(n, floor(n/2)) chains, for n=0..=10: {'PASS' if decomp_ok else 'FAIL'}")

    # --- Level matching corollary -----------------------------------------
    matching_ok = True
    for n in range(1, 9):
        chains = symmetric_chain_decomposition(n)
        for k in range(n):
            edges = level_matching(chains, k)
            expected = min(math.comb(n, k), math.comb(n, k + 1))
            ss = {s for s, _ in edges}
            ts = {t for _, t in edges}
            injective = len(ss) == len(edges) and len(ts) == len(edges)
            proper = all((s & t == s) and s != t for s, t in edges)
            if len(edges) != expected or not injective or not proper:
                matching_ok = False
                print(f"MATCHING FAIL n={n} k={k}: got {len(edges)} expected {expected} injective={injective} proper={proper}")
    print(f"B_n,k <-> B_n,k+1 matching has size min(C(n,k),C(n,k+1)) and is injective: {'PASS' if matching_ok else 'FAIL'}")

    # --- LYM inequality on random maximal antichains ----------------------
    rng = random.Random(7)
    lym_ok = True
    for _ in range(300):
        n = rng.randint(1, 8)
        antichain = random_maximal_antichain(n, rng)
        s = lym_sum(n, antichain)
        if s > 1:
            lym_ok = False
            print(f"LYM FAIL n={n}: sum={s} > 1")
    print(f"LYM inequality (sum of 1/C(n,|S|) <= 1) over random maximal antichains: {'PASS' if lym_ok else 'FAIL'}")

    # Equality case: the full middle layer.
    n = 6
    middle = [s for s in range(1 << n) if popcount(s) == n // 2]
    s = lym_sum(n, middle)
    print(f"LYM sum for the full middle layer B_{{6,3}}: {s} (should be exactly 1)")

    # --- Independent cross-check via Dilworth's theorem --------------------
    dilworth_ok = True
    for n in range(6):
        via_dilworth = max_antichain_via_dilworth(n)
        via_formula = math.comb(n, n // 2)
        if via_dilworth != via_formula:
            dilworth_ok = False
            print(f"DILWORTH FAIL n={n}: matching-based={via_dilworth} formula={via_formula}")
    print(f"max antichain size via Dilworth's theorem (bipartite matching) matches C(n,floor(n/2)) for n=0..=5: {'PASS' if dilworth_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
