"""Matchings on a path and their applications.

Vertices 0..n-1 lie on a path and carry weights x_i.  A matching of size k
covers 2k vertices; its value is the total weight of the covered vertices.
f(k) is the best value.  Everything in this file revolves around one fact:
f(k+1) - f(k) is obtained by picking the best pair of *consecutive
unmatched* vertices, so f is concave and its increments ("gains") can be
produced by a heap, by a stack in O(n), or queried on sub-intervals.
"""

# PDFBEGIN
from __future__ import annotations

import heapq
import random
from collections import defaultdict

NEG = -(1 << 62)  # "minus infinity" for max structures
INF = 1 << 62

# ---------------------------------------------------------------- Section 4
def greedy_pairs(x: list[int]) -> list[tuple[int, int, int]]:
    """Heap + linked list greedy.  Returns the list of (gain, i, j): in
    step t the best pair of neighbouring unmatched vertices is (i, j).
    O(n log n).  The gains are nonincreasing, and f(k) = gains[0..k) summed."""
    n = len(x)
    prv = list(range(-1, n - 1))
    nxt = list(range(1, n + 1))
    alive = [True] * n
    heap = [(-(x[i] + x[i + 1]), i, i + 1) for i in range(n - 1)]
    heapq.heapify(heap)
    out = []
    while len(out) < n // 2:
        neg, i, j = heapq.heappop(heap)
        if not (alive[i] and alive[j] and nxt[i] == j):
            continue  # stale entry
        out.append((-neg, i, j))
        alive[i] = alive[j] = False
        a, d = prv[i], nxt[j]
        if a >= 0:
            nxt[a] = d
        if d < n:
            prv[d] = a
        if a >= 0 and d < n:
            heapq.heappush(heap, (-(x[a] + x[d]), a, d))
    return out

# ---------------------------------------------------------------- Section 5
def gains_linear(x: list[int]) -> list[int]:
    """All floor(n/2) gains (as a multiset, unordered) in O(n) time:
    delete local maxima with a stack.  Sort descending to get f(k)."""
    st: list[int] = []
    gains: list[int] = []
    for v in x:
        st.append(v)
        while len(st) >= 3:
            b, c, d = st[-3], st[-2], st[-1]
            left_ok = len(st) < 4 or c >= st[-4]   # b+c >= a+b
            if left_ok and b >= d:                 # b+c >= c+d
                gains.append(b + c)
                del st[-3:-1]
            else:
                break
    while len(st) >= 2:                            # the top pair is a local max
        gains.append(st[-2] + st[-1])
        del st[-2:]
    return gains

def sum_of_largest(vals: list[int], k: int) -> int:
    """Sum of the k largest values in expected O(n) (quickselect)."""
    vals = list(vals)
    total, lo_k = 0, k
    while lo_k > 0 and vals:
        p = random.choice(vals)
        big = [v for v in vals if v > p]
        eq = vals.count(p)
        if len(big) >= lo_k:
            vals = big
        elif len(big) + eq >= lo_k:
            return total + sum(big) + (lo_k - len(big)) * p
        else:
            total += sum(big) + eq * p
            lo_k -= len(big) + eq
            vals = [v for v in vals if v < p]
    return total

# ---------------------------------------------------------------- Section 2-3
def edge_to_vertex_weights(w: list[int]) -> list[int]:
    """x_0 = 0, x_{i+1} = w_i - x_i, so that w_i = x_i + x_{i+1}."""
    x = [0]
    for wi in w:
        x.append(wi - x[-1])
    return x

def max_nonadjacent(b: list[int], k: int) -> int:
    """Largest sum of k pairwise non-adjacent entries of b, O(n)."""
    x = edge_to_vertex_weights(b)
    return sum_of_largest(gains_linear(x), k)

def trade_path(prices: list[int]) -> list[int]:
    """Path B_0, S_1, B_1, S_2, ..., B_{N-1}, S_N with weights +P / -P."""
    x = []
    for i in range(len(prices) - 1):
        x += [prices[i], -prices[i + 1]]
    return x

def max_profit_k_trades(prices: list[int], k: int) -> int:
    """Exactly k trades (buy on day l, sell on day r > l, next buy on day
    >= r), prices P_0..P_N.  profit(k) = -(sum of the k smallest gains)."""
    g = gains_linear(trade_path(prices))
    return sum_of_largest([-v for v in g], k)

def max_k_subarrays(a: list[int], k: int) -> int:
    """Largest total of k disjoint non-empty subarrays of a."""
    p = [0]
    for v in a:
        p.append(p[-1] + v)
    return max_profit_k_trades(p, k)

# PDFEND -- everything below (range-query structures and the test harness) is not printed in the PDF

# ---------------------------------------------------------------- Section 6
class Fenwick:
    def __init__(self, n: int):
        self.n, self.t = n, [0] * (n + 1)

    def add(self, i: int, v: int) -> None:
        i += 1
        while i <= self.n:
            self.t[i] += v
            i += i & -i

    def prefix(self, i: int) -> int:  # sum over [0, i)
        r = 0
        while i > 0:
            r += self.t[i]
            i -= i & -i
        return r

    def kth(self, k: int) -> int:  # smallest idx with prefix(idx+1) > k
        pos, step = 0, 1 << self.n.bit_length()
        while step:
            if pos + step <= self.n and self.t[pos + step] <= k:
                pos += step
                k -= self.t[pos]
            step >>= 1
        return pos

class WaveletMatrix:
    """k-th smallest in a range of a static array of small non-negative ints."""
    def __init__(self, a: list[int]):
        self.n, self.bits = len(a), max(1, max(a, default=0).bit_length())
        self.rank0, self.zeros = [], []
        cur = a[:]
        for bit in range(self.bits - 1, -1, -1):
            r0 = [0]
            for v in cur:
                r0.append(r0[-1] + (0 if (v >> bit) & 1 else 1))
            self.rank0.append(r0)
            self.zeros.append(r0[-1])
            cur = [v for v in cur if not (v >> bit) & 1] + [v for v in cur if (v >> bit) & 1]

    def kth_smallest(self, l: int, r: int, k: int) -> int:  # 0-indexed k
        res = 0
        for lvl, bit in enumerate(range(self.bits - 1, -1, -1)):
            r0 = self.rank0[lvl]
            zl, zr = r0[l], r0[r]
            if k < zr - zl:
                l, r = zl, zr
            else:
                k -= zr - zl
                res |= 1 << bit
                l, r = self.zeros[lvl] + l - zl, self.zeros[lvl] + r - zr
        return res

class SparseMin:
    """Range argmin of values vals[i] (INF = excluded), ties -> any index."""
    def __init__(self, vals: list[int]):
        self.t = [[(v, i) for i, v in enumerate(vals)]]
        j = 1
        while (1 << j) <= len(vals):
            p = self.t[-1]
            self.t.append([min(p[i], p[i + (1 << (j - 1))]) for i in range(len(vals) - (1 << j) + 1)])
            j += 1

    def argmin(self, l: int, r: int) -> int:  # over [l, r), r > l
        j = (r - l).bit_length() - 1
        return min(self.t[j][l], self.t[j][r - (1 << j)])[1]

class MaxSeg:
    """Point update, range max of (value, index)."""
    def __init__(self, n: int):
        self.sz = 1
        while self.sz < max(n, 1):
            self.sz *= 2
        self.t = [(NEG, -1)] * (2 * self.sz)

    def set(self, i: int, v: int) -> None:
        i += self.sz
        self.t[i] = (v, i - self.sz)
        i >>= 1
        while i:
            self.t[i] = max(self.t[2 * i], self.t[2 * i + 1])
            i >>= 1

    def query(self, l: int, r: int) -> tuple[int, int]:  # over [l, r)
        res = (NEG, -1)
        l += self.sz
        r += self.sz
        while l < r:
            if l & 1:
                res = max(res, self.t[l]); l += 1
            if r & 1:
                r -= 1; res = max(res, self.t[r])
            l >>= 1
            r >>= 1
        return res

class RangeMatching:
    """Offline answers to (L, R, K): best size-K matching using only the
    vertices L..R-1.  O((n + q) log n)."""
    def __init__(self, x: list[int]):
        n = len(x)
        self.n, self.x = n, x
        self.pairs = greedy_pairs(x)
        self.time = [len(self.pairs)] * n          # step at which a vertex is matched
        for t, (_, i, j) in enumerate(self.pairs):
            self.time[i] = self.time[j] = t
        self.wm = WaveletMatrix(self.time)
        self.pre = [0]
        for v in x:
            self.pre.append(self.pre[-1] + v)
        self.pmin = [SparseMin([x[i] if i % 2 == p else INF for i in range(n)]) for p in (0, 1)]
        self.alive = Fenwick(n)                    # 1 for every unmatched vertex
        self.dead_w = Fenwick(n)                   # weight of every matched vertex
        self.seg = MaxSeg(n)                       # key i: x_i + x_next(i)
        for i in range(n):
            self.alive.add(i, 1)
        for i in range(n - 1):
            self.seg.set(i, x[i] + x[i + 1])

    # --- the current global state: vertices with time < t are matched
    def _next(self, i: int) -> int:               # first unmatched index >= i, or n
        c = self.alive.prefix(i)
        return self.alive.kth(c) if c < self.alive.prefix(self.n) else self.n

    def _prev(self, i: int) -> int:               # last unmatched index <= i, or -1
        c = self.alive.prefix(i + 1)
        return self.alive.kth(c - 1) if c else -1

    def _refresh(self, i: int) -> None:
        if i < 0:
            return
        j = self._next(i + 1) if self._is_alive(i) else self.n
        self.seg.set(i, self.x[i] + self.x[j] if j < self.n else NEG)

    def _is_alive(self, i: int) -> bool:
        return self.alive.prefix(i + 1) - self.alive.prefix(i) == 1

    def _flip(self, i: int, to_alive: bool) -> None:
        d = 1 if to_alive else -1
        self.alive.add(i, d)
        self.dead_w.add(i, -d * self.x[i])
        self._refresh(i)
        self._refresh(self._prev(i - 1))

    def _answer(self, L: int, R: int, K: int) -> int:
        a = self.alive
        left_cross = (L - a.prefix(L)) & 1                      # a matching edge crosses L
        right_cross = ((self.n - R) - (a.prefix(self.n) - a.prefix(R))) & 1
        u, v = self._next(L), self._prev(R - 1)                 # first / last unmatched in range
        log = []                                                # (vertex, became unmatched?)
        if left_cross:                                          # release b_L
            log.append((self.pmin[L % 2].argmin(L, u), True))
        if right_cross:                                         # release b_R
            log.append((self.pmin[(R - 1) % 2].argmin(v + 1, R), True))
        for b, _ in log:
            self._flip(b, True)
        s = (R - L) - (a.prefix(R) - a.prefix(L))              # matched vertices inside
        assert s % 2 == 0 and s <= 2 * K
        for _ in range((2 * K - s) // 2):                       # greedy on what is left
            p = self._prev(R - 1)
            _, i = self.seg.query(L, p)
            j = self._next(i + 1)
            self._flip(i, False)
            self._flip(j, False)
            log += [(i, False), (j, False)]
        ans = self.dead_w.prefix(R) - self.dead_w.prefix(L)
        for b, was_released in reversed(log):                   # undo everything
            self._flip(b, not was_released)
        return ans

    def answer(self, queries: list[tuple[int, int, int]]) -> list[int]:
        res = [0] * len(queries)
        by_time = defaultdict(list)
        for idx, (L, R, K) in enumerate(queries):
            if K == 0:
                res[idx] = 0
            elif 2 * K == R - L:
                res[idx] = self.pre[R] - self.pre[L]
            else:  # time at which the range gets its (2K)-th matched vertex
                by_time[self.wm.kth_smallest(L, R, 2 * K - 1)].append(idx)
        for t in range(len(self.pairs) + 1):
            for idx in by_time[t]:
                res[idx] = self._answer(*queries[idx])
            if t < len(self.pairs):
                _, i, j = self.pairs[t]
                self._flip(i, False); self._flip(j, False)
        return res


# ------------------------------------------------------------ brute force
def brute_f(x: list[int]) -> list[int]:
    """f(0..n//2) by the obvious DP over prefixes: dp[i][k]."""
    n = len(x)
    dp = [[NEG] * (n // 2 + 1) for _ in range(n + 1)]
    dp[0][0] = 0
    for i in range(1, n + 1):
        for k in range(n // 2 + 1):
            best = dp[i - 1][k]
            if i >= 2 and k >= 1 and dp[i - 2][k - 1] > NEG:
                best = max(best, dp[i - 2][k - 1] + x[i - 1] + x[i - 2])
            dp[i][k] = best
    return dp[n]

def brute_nonadjacent(b: list[int], k: int) -> int:
    n = len(b)
    dp = [[NEG] * (k + 1) for _ in range(n + 2)]
    for i in range(n + 2):
        dp[i][0] = 0
    for i in range(1, n + 1):
        for j in range(1, k + 1):
            dp[i + 1][j] = dp[i][j]
            if dp[i - 1][j - 1] > NEG:
                dp[i + 1][j] = max(dp[i + 1][j], dp[i - 1][j - 1] + b[i - 1])
    return dp[n + 1][k]

def brute_subarrays(a: list[int], k: int) -> int:
    """Exactly k disjoint non-empty subarrays: states (used, inside)."""
    n = len(a)
    dp = [[NEG, NEG] for _ in range(k + 1)]   # dp[j][in] = j subarrays started
    dp[0][0] = 0
    for v in a:
        nd = [[NEG, NEG] for _ in range(k + 1)]
        for j in range(k + 1):
            for inside in (0, 1):
                cur = dp[j][inside]
                if cur == NEG:
                    continue
                nd[j][0] = max(nd[j][0], cur)                      # stay outside / close
                if inside:
                    nd[j][1] = max(nd[j][1], cur + v)              # extend
                if j < k:
                    nd[j + 1][1] = max(nd[j + 1][1], cur + v)      # start a new one
        dp = nd
    return max(dp[k])

def brute_range(x: list[int], L: int, R: int, K: int) -> int:
    return brute_f(x[L:R])[K]

def main() -> None:
    random.seed(12345)
    checked = {"heap": 0, "stack": 0, "reduce": 0, "range": 0}

    # Sections 4 and 5: the three ways of computing f agree with the DP.
    for _ in range(3000):
        n = random.randint(0, 14)
        hi = random.choice([1, 2, 5, 50])           # tiny ranges force many ties
        x = [random.randint(-hi, hi) for _ in range(n)]
        want = brute_f(x)
        pairs = greedy_pairs(x)
        gains = [g for g, _, _ in pairs]
        assert all(gains[i] >= gains[i + 1] for i in range(len(gains) - 1)), x
        f = [0]
        for g in gains:
            f.append(f[-1] + g)
        assert f == want, (x, f, want)
        checked["heap"] += 1
        lin = sorted(gains_linear(x), reverse=True)
        assert lin == gains, (x, lin, gains)
        for k in range(n // 2 + 1):
            assert sum_of_largest(gains_linear(x), k) == want[k]
        checked["stack"] += 1

    # Section 3: the reductions.
    for _ in range(2000):
        n = random.randint(1, 12)
        hi = random.choice([2, 6, 40])
        b = [random.randint(-hi, hi) for _ in range(n)]
        for k in range((n + 1) // 2 + 1):
            assert max_nonadjacent(b, k) == brute_nonadjacent(b, k), (b, k)
        a = [random.randint(-hi, hi) for _ in range(n)]
        for k in range(n + 1):
            assert max_k_subarrays(a, k) == brute_subarrays(a, k), (a, k)
        p = [random.randint(0, hi) for _ in range(n + 1)]
        d = [p[i + 1] - p[i] for i in range(n)]
        for k in range(n + 1):
            assert max_profit_k_trades(p, k) == brute_subarrays(d, k), (p, k)
        checked["reduce"] += 1
    w = [4, -7, 2, 9, 1]
    x = edge_to_vertex_weights(w)
    assert [x[i] + x[i + 1] for i in range(len(w))] == w

    # Section 6: range queries against the brute force on the slice.
    for _ in range(600):
        n = random.randint(2, 16)
        hi = random.choice([1, 3, 9, 100])
        x = [random.randint(-hi, hi) for _ in range(n)]
        queries = []
        for L in range(n):
            for R in range(L + 1, n + 1):
                for K in range((R - L) // 2 + 1):
                    queries.append((L, R, K))
        random.shuffle(queries)
        got = RangeMatching(x).answer(queries)
        for q, g in zip(queries, got):
            assert g == brute_range(x, *q), (x, q, g)
            checked["range"] += 1

    # A larger instance, just to see the structures cope.
    n = 3000
    x = [random.randint(-10**6, 10**6) for _ in range(n)]
    queries = []
    for _ in range(300):
        L = random.randint(0, n - 2)
        R = random.randint(L + 2, n)
        queries.append((L, R, random.randint(1, (R - L) // 2)))
    got = RangeMatching(x).answer(queries)
    for _ in range(5):
        i = random.randrange(len(queries))
        assert got[i] == brute_range(x, *queries[i])
    print("all checks passed:", checked)

if __name__ == "__main__":
    main()
