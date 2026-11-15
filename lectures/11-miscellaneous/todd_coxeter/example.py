"""The Todd-Coxeter algorithm for monoids. Given a finite alphabet and a
finite list of relations u = v between words, it builds the right Cayley
graph of the monoid they define, provided that monoid is finite. The graph
is a transition table T[v][x]: the element reached from element v by
multiplying on the right with the generator x, with the identity numbered 0.

The algorithm keeps a growing set of provisional vertices, creates a vertex
whenever a needed edge is undefined, and merges (with a union-find
structure) vertices that must stand for the same element. If the monoid is
finite it always stops; if it is infinite it never does, so the code takes
a cap on the number of vertices and reports when that cap is hit.
"""

# PDFBEGIN
from __future__ import annotations

from collections import deque

UNDEFINED = -1

Word = list[int]
Relation = tuple[Word, Word]


class Enumerator:
    def __init__(self, num_letters: int, relations: list[Relation]) -> None:
        self.k = num_letters
        self.relations = relations
        self.table: list[list[int]] = []   # table[v][x] = target vertex or UNDEFINED
        self.parent: list[int] = []        # union-find forest over vertices
        self.done: list[bool] = []         # class fully processed?
        self.created = 0                   # total vertices ever created

    # ---- union-find (the smaller index always becomes the representative) ----
    def find(self, v: int) -> int:
        while self.parent[v] != v:
            self.parent[v] = self.parent[self.parent[v]]
            v = self.parent[v]
        return v

    def add_vertex(self) -> int:
        v = len(self.table)
        self.table.append([UNDEFINED] * self.k)
        self.parent.append(v)
        self.done.append(False)
        self.created += 1
        return v

    def step(self, v: int, x: int) -> int:
        """Follow the edge labelled x from v, creating a vertex if it is undefined."""
        v = self.find(v)
        if self.table[v][x] == UNDEFINED:
            self.table[v][x] = self.add_vertex()
        return self.find(self.table[v][x])

    def walk(self, v: int, word: Word) -> int:
        for x in word:
            v = self.step(v, x)
        return self.find(v)

    def merge(self, p: int, q: int) -> None:
        """Declare vertices p and q to be the same element, and propagate."""
        pending = deque([(p, q)])
        while pending:
            p, q = pending.popleft()
            p, q = self.find(p), self.find(q)
            if p == q:
                continue
            if q < p:
                p, q = q, p                       # keep the smaller index as representative
            self.parent[q] = p
            self.done[p] = self.done[p] or self.done[q]
            for x in range(self.k):
                a, b = self.table[p][x], self.table[q][x]
                if a == UNDEFINED:
                    self.table[p][x] = b
                elif b != UNDEFINED:
                    pending.append((a, b))        # same element times x: targets coincide

    def run(self, max_vertices: int) -> list[list[int]] | None:
        """Returns the transition table (identity = 0, numbered in breadth-first
        order), or None if more than max_vertices vertices were needed."""
        self.add_vertex()
        p = 0
        while p < len(self.table):
            if len(self.table) > max_vertices:
                return None
            if not self.done[self.find(p)]:
                for u, v in self.relations:
                    self.merge(self.walk(p, u), self.walk(p, v))
                for x in range(self.k):
                    self.step(p, x)
                self.done[self.find(p)] = True
            p += 1
        return self.canonical_table()

    def canonical_table(self) -> list[list[int]]:
        """Renumber the surviving classes in breadth-first order from the identity."""
        root = self.find(0)
        label = {root: 0}
        queue = deque([root])
        order = []
        while queue:
            v = queue.popleft()
            order.append(v)
            for x in range(self.k):
                w = self.find(self.table[v][x])
                if w not in label:
                    label[w] = len(label)
                    queue.append(w)
        return [[label[self.find(self.table[v][x])] for x in range(self.k)] for v in order]


def todd_coxeter(num_letters: int, relations: list[Relation], max_vertices: int = 100_000) -> list[list[int]] | None:
    return Enumerator(num_letters, relations).run(max_vertices)

# PDFEND -- everything below is the test harness and is not printed in the PDF


def follow(table: list[list[int]], word: Word, start: int = 0) -> int:
    v = start
    for x in word:
        v = table[v][x]
    return v


def parse(text: str, letters: str) -> Word:
    return [letters.index(ch) for ch in text if ch != "e"]   # 'e' denotes the empty word


def closure_of_permutations(generators: list[tuple[int, ...]]) -> dict[tuple[int, ...], int]:
    """All permutations generated (as a monoid) by the given ones, by breadth-first search."""
    size = len(generators[0])
    identity = tuple(range(size))
    seen = {identity: 0}
    queue = deque([identity])
    while queue:
        g = queue.popleft()
        for h in generators:
            gh = tuple(g[h[i]] for i in range(size))   # (g h)(i) = g(h(i))
            if gh not in seen:
                seen[gh] = len(seen)
                queue.append(gh)
    return seen


def check_group_presentation(name: str, letters: str, relations: list[tuple[str, str]], generators: list[tuple[int, ...]], expected: int) -> bool:
    rels = [(parse(u, letters), parse(v, letters)) for u, v in relations]
    enum = Enumerator(len(letters), rels)
    table = enum.run(100_000)
    if table is None or len(table) != expected:
        print(f"  {name}: FAIL (got {None if table is None else len(table)}, expected {expected})")
        return False
    # The map word -> permutation must be a bijection from the enumerated elements
    # onto the permutation group, which proves the presentation defines that group.
    word_of = {0: []}
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for x in range(len(letters)):
            w = table[v][x]
            if w not in word_of:
                word_of[w] = word_of[v] + [x]
                queue.append(w)
    images = set()
    for v, word in word_of.items():
        perm = tuple(range(len(generators[0])))
        for x in word:
            perm = tuple(perm[generators[x][i]] for i in range(len(perm)))
        images.add(perm)
    group = closure_of_permutations(generators)
    ok = len(images) == expected == len(group) and images == set(group)
    print(f"  {name}: {len(table)} elements, {enum.created} vertices created, isomorphic to the permutation group: {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> None:
    all_ok = True

    # --- Worked example: a small monoid that is not a group ------------------
    rels = [(parse("aa", "ab"), parse("a", "ab")), (parse("bb", "ab"), parse("b", "ab")), (parse("ba", "ab"), parse("ab", "ab"))]
    enum = Enumerator(2, rels)
    table = enum.run(1000)
    print("<a,b | aa=a, bb=b, ba=ab>: transition table (rows = elements, columns = a, b):")
    for v, row in enumerate(table):
        print(f"  {v}: {row}")
    all_ok &= table == [[1, 2], [1, 3], [3, 2], [3, 3]]
    print(f"  4 elements ({enum.created} vertices created along the way): {'PASS' if len(table) == 4 else 'FAIL'}")

    # --- Groups given by presentations, checked against permutation groups --------
    print("\ngroups from presentations:")
    s3 = [(1, 0, 2), (0, 2, 1)]
    all_ok &= check_group_presentation("S3 = <a,b | aaa=e, bb=e, ba=aab>", "ab", [("aaa", "e"), ("bb", "e"), ("ba", "aab")], [(1, 2, 0), (1, 0, 2)], 6)
    all_ok &= check_group_presentation("S3 = <a,b | aa=e, bb=e, aba=bab>", "ab", [("aa", "e"), ("bb", "e"), ("aba", "bab")], s3, 6)
    for n in range(3, 9):
        rot = tuple((i + 1) % n for i in range(n))
        refl = tuple((-i) % n for i in range(n))
        all_ok &= check_group_presentation(f"dihedral D{n}", "ab", [("a" * n, "e"), ("bb", "e"), ("ba", "a" * (n - 1) + "b")], [rot, refl], 2 * n)
    for n, expected in ((4, 24), (5, 120)):
        letters = "abcdefgh"[: n - 1]
        rels_text = [(c + c, "e") for c in letters]
        for i in range(n - 1):
            for j in range(i + 1, n - 1):
                rels_text.append(((letters[i] + letters[j]) * (3 if j == i + 1 else 2), "e"))
        gens = []
        for i in range(n - 1):
            perm = list(range(n))
            perm[i], perm[i + 1] = perm[i + 1], perm[i]
            gens.append(tuple(perm))
        all_ok &= check_group_presentation(f"symmetric group S{n} (Coxeter presentation)", letters, rels_text, gens, expected)

    # --- Infinite monoids never finish, so a cap is needed ---------------------
    bicyclic = Enumerator(2, [(parse("ab", "ab"), parse("e", "ab"))])
    capped = bicyclic.run(300) is None
    print(f"\n<a,b | ab=e> is infinite; the run is stopped by the vertex cap: {'PASS' if capped else 'FAIL'}")
    all_ok &= capped

    # --- Many provisional vertices, one final element ----------------------------
    for m in (10, 50, 200):
        e = Enumerator(1, [([0] * m, []), ([0], [])])
        t = e.run(10_000)
        print(f"<a | a^{m}=e, a=e>: {len(t)} element, {e.created} vertices created")
        all_ok &= t == [[0]]

    # --- Random presentations: soundness checks that need no known answer ----------
    state = 12345

    def rnd(bound: int) -> int:
        nonlocal state
        state = (state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        return (state >> 33) % bound

    finite = 0
    tested = 200
    problems = 0
    for _ in range(tested):
        k = 2 + rnd(2)
        rels = []
        for _ in range(1 + rnd(3)):
            u = [rnd(k) for _ in range(rnd(4))]
            v = [rnd(k) for _ in range(rnd(4))]
            rels.append((u, v))
        table = todd_coxeter(k, rels, 400)
        if table is None:
            continue
        finite += 1
        n = len(table)
        # (1) every relation holds at every element
        for p in range(n):
            for u, v in rels:
                if follow(table, u, p) != follow(table, v, p):
                    problems += 1
        # (2) the multiplication read off the graph is associative
        rep = {0: []}
        queue = deque([0])
        while queue:
            v = queue.popleft()
            for x in range(k):
                w = table[v][x]
                if w not in rep:
                    rep[w] = rep[v] + [x]
                    queue.append(w)
        mult = [[follow(table, rep[h], g) for h in range(n)] for g in range(n)]
        for g in range(n):
            for h in range(n):
                for i in range(n):
                    if mult[mult[g][h]][i] != mult[g][mult[h][i]]:
                        problems += 1
        # (3) words connected by rewriting with the relations end at the same element
        for length in range(0, 6):
            for code in range(k ** length):
                word = [(code // k ** j) % k for j in range(length)]
                for u, v in rels:
                    for pos in range(length - len(u) + 1):
                        if word[pos: pos + len(u)] == u:
                            other = word[:pos] + v + word[pos + len(u):]
                            if follow(table, word) != follow(table, other):
                                problems += 1
    print(f"\n{finite} of {tested} random presentations define finite monoids (cap 400 vertices); relations hold, multiplication is associative, rewriting is respected: {'PASS' if problems == 0 else 'FAIL'}")
    all_ok &= problems == 0

    print(f"\nAll checks: {'PASS' if all_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
