"""Pruefer code: a bijection between labeled trees on {0,...,n-1} and
sequences in {0,...,n-1}^(n-2).  Encoding repeatedly deletes the smallest
leaf and records its neighbor; decoding rebuilds the tree from the fact
that a vertex occurring c times in the code has degree c+1.  Consequences
(Cayley's formula, degree-sequence counts, forest and bipartite counts) are
checked by brute force in the harness below.
"""

# PDFBEGIN
from __future__ import annotations

import heapq


def encode_heap(n: int, edges: list[tuple[int, int]]) -> list[int]:
    """Pruefer code of a labeled tree, O(n log n) with a min-heap."""
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    deg = [len(adj[v]) for v in range(n)]
    leaves = [v for v in range(n) if deg[v] == 1]
    heapq.heapify(leaves)
    code = []
    for _ in range(n - 2):
        u = heapq.heappop(leaves)  # smallest current leaf
        v = next(w for w in adj[u] if deg[w] > 0)  # its only live neighbor
        code.append(v)
        deg[u] = 0
        deg[v] -= 1
        if deg[v] == 1:
            heapq.heappush(leaves, v)
    return code


def decode_heap(n: int, code: list[int]) -> list[tuple[int, int]]:
    """Inverse of encode_heap, O(n log n)."""
    left = [0] * n  # occurrences of each vertex in the unread part of code
    for v in code:
        left[v] += 1
    free = [v for v in range(n) if left[v] == 0]  # current leaves
    heapq.heapify(free)
    edges = []
    for v in code:
        u = heapq.heappop(free)
        edges.append((u, v))
        left[v] -= 1
        if left[v] == 0:
            heapq.heappush(free, v)
    a = heapq.heappop(free)
    b = heapq.heappop(free)
    edges.append((a, b))
    return edges


def encode_linear(n: int, edges: list[tuple[int, int]]) -> list[int]:
    """Pruefer code in O(n): root the tree at n-1 and sweep a pointer."""
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    parent = [-1] * n
    stack = [n - 1]
    while stack:  # orient every edge toward the root n-1
        v = stack.pop()
        for w in adj[v]:
            if w != parent[v]:
                parent[w] = v
                stack.append(w)
    deg = [len(adj[v]) for v in range(n)]
    ptr = 0
    while deg[ptr] != 1:
        ptr += 1
    leaf = ptr  # smallest current leaf
    code = []
    for _ in range(n - 2):
        p = parent[leaf]
        code.append(p)
        deg[p] -= 1
        if deg[p] == 1 and p < ptr:
            leaf = p  # p became a leaf smaller than everything unseen
        else:
            ptr += 1
            while deg[ptr] != 1:
                ptr += 1
            leaf = ptr
    return code


def decode_linear(n: int, code: list[int]) -> list[tuple[int, int]]:
    """Inverse of encode_linear, O(n)."""
    deg = [1] * n
    for v in code:
        deg[v] += 1
    ptr = 0
    while deg[ptr] != 1:
        ptr += 1
    leaf = ptr
    edges = []
    for v in code:
        edges.append((leaf, v))
        deg[v] -= 1
        if deg[v] == 1 and v < ptr:
            leaf = v
        else:
            ptr += 1
            while deg[ptr] != 1:
                ptr += 1
            leaf = ptr
    edges.append((leaf, n - 1))
    return edges

# PDFEND -- everything below is the test harness and is not printed in the PDF

import itertools
import math
import random
from collections import Counter


def canon(edges):
    return sorted((min(a, b), max(a, b)) for a, b in edges)


def is_tree(n, edges):
    if len(edges) != n - 1:
        return False
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        parent[ra] = rb
    return True


def all_trees(n):
    """Every labeled tree on n vertices, by brute force over edge subsets."""
    pairs = list(itertools.combinations(range(n), 2))
    for sub in itertools.combinations(pairs, n - 1):
        if is_tree(n, sub):
            yield list(sub)


def random_tree(n, rng):
    """Random labeled tree NOT built through Pruefer decoding (so that the
    round-trip tests are not circular): attach each vertex to an earlier
    one, then relabel at random."""
    perm = list(range(n))
    rng.shuffle(perm)
    return [(perm[i], perm[rng.randrange(i)]) for i in range(1, n)]


def bip_encode(n, edges, side):
    """Two-sequence Pruefer code of a spanning tree of a complete bipartite
    graph (side[v] in {0,1}): delete the smallest leaf until one vertex per
    side remains; the recorded neighbor goes to the sequence of ITS side."""
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    alive = set(range(n))
    seqs = [[], []]
    while len(alive) > 2:
        u = min(v for v in alive if len(adj[v]) == 1)
        (v,) = adj[u]
        seqs[side[v]].append(v)
        adj[v].discard(u)
        alive.discard(u)
        adj[u].clear()
    return tuple(seqs[0]), tuple(seqs[1])


def bip_decode(n, seqs, side):
    """Inverse of bip_encode: the smallest vertex absent from the unread part
    of its own side's sequence is the next leaf; its neighbor is the next
    unread entry of the OTHER side's sequence."""
    from collections import deque
    q = [deque(seqs[0]), deque(seqs[1])]
    left = [0] * n
    for s in (0, 1):
        for v in q[s]:
            left[v] += 1
    alive = set(range(n))
    edges = []
    while len(alive) > 2:
        u = min(v for v in alive if left[v] == 0)
        v = q[1 - side[u]].popleft()
        left[v] -= 1
        edges.append((u, v))
        alive.discard(u)
    a, b = sorted(alive)
    edges.append((a, b))
    return edges


def main():
    rng = random.Random(20240601)

    # --- the lecture's worked examples ------------------------------------
    ex_edges = [(0, 4), (1, 4), (4, 7), (2, 7), (7, 3), (3, 6), (5, 6), (6, 8)]
    ex_code = encode_heap(9, ex_edges)
    print("worked encode example (n=9):", ex_code)
    assert ex_code == [4, 4, 7, 7, 6, 3, 6] == encode_linear(9, ex_edges)
    assert canon(decode_heap(9, ex_code)) == canon(ex_edges)
    ex2 = [3, 3, 0, 5, 0, 3]
    t2 = decode_heap(8, ex2)
    print("worked decode example (n=8):", canon(t2))
    assert canon(t2) == [(0, 3), (0, 4), (0, 5), (1, 3), (2, 3), (3, 7), (5, 6)]
    assert encode_heap(8, t2) == ex2

    # --- bijection and Cayley's formula, n = 2..7 -------------------------
    for n in range(2, 8):
        trees = list(all_trees(n))
        codes = set()
        for t in trees:
            c = encode_heap(n, t)
            assert c == encode_linear(n, t)
            assert canon(decode_heap(n, c)) == canon(t)
            assert canon(decode_linear(n, c)) == canon(t)
            codes.add(tuple(c))
        expect = n ** (n - 2) if n >= 2 else 1
        assert len(trees) == expect == len(codes), (n, len(trees), expect)
        print(f"n={n}: {len(trees)} trees = n^(n-2), encode is a bijection")

    # --- decode every sequence, n <= 7 ------------------------------------
    for n in range(2, 8):
        for code in itertools.product(range(n), repeat=n - 2):
            code = list(code)
            t = decode_heap(n, code)
            assert is_tree(n, t)
            assert encode_heap(n, t) == code
            assert canon(t) == canon(decode_linear(n, code))
    print("decode gives a tree and encode inverts it for every sequence, n<=7")

    # --- degree lemma and linear/heap agreement on random big trees ------
    for _ in range(300):
        n = rng.randrange(2, 200)
        t = random_tree(n, rng)
        c = encode_heap(n, t)
        assert c == encode_linear(n, t)
        assert canon(decode_linear(n, c)) == canon(t)
        deg = Counter()
        for a, b in t:
            deg[a] += 1
            deg[b] += 1
        cnt = Counter(c)
        assert all(deg[v] == 1 + cnt[v] for v in range(n))
    print("degree lemma deg(v) = 1 + #occurrences, linear == heap: PASS")

    # --- degree sequences (multinomial formula) ---------------------------
    for n in range(3, 8):
        by_deg = Counter()
        for t in all_trees(n):
            d = [0] * n
            for a, b in t:
                d[a] += 1
                d[b] += 1
            by_deg[tuple(d)] += 1
        for d, cnt in by_deg.items():
            f = math.factorial(n - 2)
            for x in d:
                f //= math.factorial(x - 1)
            assert cnt == f, (n, d, cnt, f)
        # and every admissible degree sequence occurs
        admissible = [d for d in itertools.product(range(1, n), repeat=n)
                      if sum(d) == 2 * (n - 1)]
        assert set(admissible) == set(by_deg), n
    print("degree-sequence multinomial formula, n<=7: PASS")

    # --- trees where vertex 0 has degree d --------------------------------
    for n in range(3, 8):
        by = Counter()
        for t in all_trees(n):
            by[sum(1 for e in t if 0 in e)] += 1
        for d in range(1, n):
            assert by[d] == math.comb(n - 2, d - 1) * (n - 1) ** (n - 1 - d)
    print("trees with deg(0)=d count C(n-2,d-1)(n-1)^(n-1-d): PASS")

    # --- forests -> trees: N^(k-2) * prod(n_i) ----------------------------
    def components(n, edges):
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b in edges:
            parent[find(a)] = find(b)
        comp = {}
        for v in range(n):
            comp.setdefault(find(v), []).append(v)
        return list(comp.values())

    for n in range(3, 8):
        pairs = list(itertools.combinations(range(n), 2))
        checked = 0
        for m in range(0, n - 2):  # forests with k = n-m >= 3 components
            for fe in itertools.combinations(pairs, m):
                if len(components(n, fe)) != n - m:
                    continue  # has a cycle
                comps = components(n, fe)
                k = len(comps)
                pset = set(fe)
                ways = 0
                for add in itertools.combinations(
                        [p for p in pairs if p not in pset], k - 1):
                    if is_tree(n, list(fe) + list(add)):
                        ways += 1
                prod = 1
                for c in comps:
                    prod *= len(c)
                assert ways == n ** (k - 2) * prod
                checked += 1
                if checked >= 40:
                    break
            if checked >= 40:
                break
    print("forest-to-tree count N^(k-2) * prod n_i: PASS")

    # --- rooted forests with prescribed roots {0..k-1} --------------------
    for n in range(2, 8):
        pairs = list(itertools.combinations(range(n), 2))
        for k in range(1, n + 1):
            count = 0
            for m in range(n - k, n - k + 1):  # forest with k comps has n-k edges
                for fe in itertools.combinations(pairs, m):
                    comps = components(n, fe)
                    if len(comps) != k:
                        continue
                    if all(sum(1 for v in c if v < k) == 1 for c in comps):
                        count += 1
            assert count == k * n ** (n - k - 1), (n, k, count)
    print("rooted forests with root set {0..k-1}: k n^(n-k-1): PASS")

    # --- spanning trees of complete bipartite graphs ----------------------
    # brute-force count AND the two-sequence code, with random side labelings
    for a in range(1, 5):
        for b in range(1, 5):
            n = a + b
            if n < 3 or n > 7:
                continue
            labels = list(range(n))
            rng.shuffle(labels)
            part_a = set(labels[:a])
            side = [0 if v in part_a else 1 for v in range(n)]
            pairs = [(i, j) for i in part_a for j in range(n) if j not in part_a]
            trees = [list(sub) for sub in itertools.combinations(pairs, n - 1)
                     if is_tree(n, sub)]
            codes = set()
            for t in trees:
                c = bip_encode(n, t, side)
                assert len(c[0]) == b - 1 and len(c[1]) == a - 1
                assert canon(bip_decode(n, c, side)) == canon(t)
                codes.add(c)
            expect = a ** (b - 1) * b ** (a - 1)
            assert len(trees) == len(codes) == expect, (a, b, len(trees))
    print("spanning trees of K_{a,b} = a^(b-1) b^(a-1), two-sequence code "
          "is a bijection: PASS")

    # --- the specific numbers quoted in the lecture's worked examples -----
    # degree sequence (4,1,2,1,1,3,1,1) on 8 vertices: 6!/(3! 1! 2!) = 60
    want = (4, 1, 2, 1, 1, 3, 1, 1)
    hits = 0
    for code in itertools.product(range(8), repeat=6):
        cnt = Counter(code)
        if all(1 + cnt[v] == want[v] for v in range(8)):
            hits += 1
    assert hits == 60
    # forest {0}, {1,2}, {3,4,5,6} on 7 vertices: 7^(3-2) * 1*2*4 = 56
    fe = [(1, 2), (3, 4), (4, 5), (5, 6)]
    rest = [p for p in itertools.combinations(range(7), 2) if p not in fe]
    ways = sum(1 for add in itertools.combinations(rest, 2)
               if is_tree(7, fe + list(add)))
    assert ways == 56
    # K_{2,3} and K_{3,3}
    for a, b, want_cnt in ((2, 3, 12), (3, 3, 81)):
        pairs = [(i, j) for i in range(a) for j in range(a, a + b)]
        got = sum(1 for sub in itertools.combinations(pairs, a + b - 1)
                  if is_tree(a + b, sub))
        assert got == want_cnt
    print("numbers quoted in the notes (60, 56, 12, 81): PASS")

    print("ALL PASS")


if __name__ == "__main__":
    main()
