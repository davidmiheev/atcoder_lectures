// Matchings on a path and their applications.
//
// Vertices 0..n-1 lie on a path and carry weights x_i. A matching of size k
// covers 2k vertices; its value is the total weight of the covered vertices.
// f(k) is the best value. Everything here revolves around one fact:
// f(k+1) - f(k) is obtained by picking the best pair of *consecutive
// unmatched* vertices, so f is concave and its increments ("gains") can be
// produced by a heap, by a stack in O(n), or queried on sub-intervals.

// PDFBEGIN
use std::cmp::Reverse;
use std::collections::BinaryHeap;

const NEG: i64 = -(1 << 62); // "minus infinity" for max structures
const INF: i64 = 1 << 62;

// ---------------------------------------------------------------- Section 4
/// Heap + linked list greedy. Returns (gain, i, j) per step: the best pair of
/// neighbouring unmatched vertices is (i, j). O(n log n). The gains are
/// nonincreasing and f(k) is the sum of the first k of them.
fn greedy_pairs(x: &[i64]) -> Vec<(i64, usize, usize)> {
    let n = x.len();
    let mut prv: Vec<isize> = (0..n as isize).map(|i| i - 1).collect();
    let mut nxt: Vec<usize> = (1..=n).collect();
    let mut alive = vec![true; n];
    let mut heap: BinaryHeap<(i64, Reverse<usize>, usize)> = (0..n.saturating_sub(1))
        .map(|i| (x[i] + x[i + 1], Reverse(i), i + 1))
        .collect();
    let mut out = Vec::new();
    while out.len() < n / 2 {
        let (g, Reverse(i), j) = heap.pop().unwrap();
        if !(alive[i] && alive[j] && nxt[i] == j) {
            continue; // stale entry
        }
        out.push((g, i, j));
        alive[i] = false;
        alive[j] = false;
        let (a, d) = (prv[i], nxt[j]);
        if a >= 0 {
            nxt[a as usize] = d;
        }
        if d < n {
            prv[d] = a;
        }
        if a >= 0 && d < n {
            heap.push((x[a as usize] + x[d], Reverse(a as usize), d));
        }
    }
    out
}

// ---------------------------------------------------------------- Section 5
/// All floor(n/2) gains (as an unordered multiset) in O(n): delete local
/// maxima with a stack. Sort descending to obtain f(k).
fn gains_linear(x: &[i64]) -> Vec<i64> {
    let mut st: Vec<i64> = Vec::new();
    let mut gains = Vec::new();
    for &v in x {
        st.push(v);
        while st.len() >= 3 {
            let m = st.len();
            let (b, c, d) = (st[m - 3], st[m - 2], st[m - 1]);
            let left_ok = m < 4 || c >= st[m - 4]; // b+c >= a+b
            if left_ok && b >= d {                 // b+c >= c+d
                gains.push(b + c);
                st.truncate(m - 3);
                st.push(d);
            } else {
                break;
            }
        }
    }
    while st.len() >= 2 {                          // the top pair is a local max
        let m = st.len();
        gains.push(st[m - 2] + st[m - 1]);
        st.truncate(m - 2);
    }
    gains
}

/// Sum of the k largest values in O(n) (introselect from the standard library).
fn sum_of_largest(mut vals: Vec<i64>, k: usize) -> i64 {
    if k == 0 {
        return 0;
    }
    if k < vals.len() {
        vals.select_nth_unstable_by(k - 1, |a, b| b.cmp(a));
    }
    vals[..k.min(vals.len())].iter().sum()
}

// ---------------------------------------------------------------- Section 2-3
/// x_0 = 0, x_{i+1} = w_i - x_i, so that w_i = x_i + x_{i+1}.
fn edge_to_vertex_weights(w: &[i64]) -> Vec<i64> {
    let mut x = vec![0];
    for &wi in w {
        let last = *x.last().unwrap();
        x.push(wi - last);
    }
    x
}

/// Largest sum of k pairwise non-adjacent entries of b, O(n).
fn max_nonadjacent(b: &[i64], k: usize) -> i64 {
    sum_of_largest(gains_linear(&edge_to_vertex_weights(b)), k)
}

/// Path B_0, S_1, B_1, S_2, ..., B_{N-1}, S_N with weights +P / -P.
fn trade_path(prices: &[i64]) -> Vec<i64> {
    let mut x = Vec::new();
    for i in 0..prices.len() - 1 {
        x.push(prices[i]);
        x.push(-prices[i + 1]);
    }
    x
}

/// Exactly k trades on prices P_0..P_N: profit(k) = -(sum of the k smallest gains).
fn max_profit_k_trades(prices: &[i64], k: usize) -> i64 {
    let g: Vec<i64> = gains_linear(&trade_path(prices)).into_iter().map(|v| -v).collect();
    sum_of_largest(g, k)
}

/// Largest total of k disjoint non-empty subarrays of a.
fn max_k_subarrays(a: &[i64], k: usize) -> i64 {
    let mut p = vec![0];
    for &v in a {
        let last = *p.last().unwrap();
        p.push(last + v);
    }
    max_profit_k_trades(&p, k)
}

// ---------------------------------------------------------------- Section 6
struct Fenwick {
    n: usize,
    t: Vec<i64>,
}
impl Fenwick {
    fn new(n: usize) -> Self {
        Fenwick { n, t: vec![0; n + 1] }
    }
    fn add(&mut self, i: usize, v: i64) {
        let mut i = i + 1;
        while i <= self.n {
            self.t[i] += v;
            i += i & i.wrapping_neg();
        }
    }
    fn prefix(&self, i: usize) -> i64 {
        let (mut i, mut r) = (i, 0);
        while i > 0 {
            r += self.t[i];
            i -= i & i.wrapping_neg();
        }
        r
    }
    /// Smallest idx with prefix(idx + 1) > k.
    fn kth(&self, mut k: i64) -> usize {
        let mut pos = 0;
        let mut step = 1usize << (usize::BITS - self.n.leading_zeros());
        while step > 0 {
            if pos + step <= self.n && self.t[pos + step] <= k {
                pos += step;
                k -= self.t[pos];
            }
            step >>= 1;
        }
        pos
    }
}

/// k-th smallest in a range of a static array of small non-negative ints.
struct WaveletMatrix {
    bits: usize,
    rank0: Vec<Vec<usize>>,
    zeros: Vec<usize>,
}
impl WaveletMatrix {
    fn new(a: &[usize]) -> Self {
        let mx = a.iter().copied().max().unwrap_or(0);
        let bits = ((usize::BITS - mx.leading_zeros()) as usize).max(1);
        let (mut rank0, mut zeros) = (Vec::new(), Vec::new());
        let mut cur = a.to_vec();
        for bit in (0..bits).rev() {
            let mut r0 = vec![0];
            for &v in &cur {
                r0.push(r0.last().unwrap() + usize::from((v >> bit) & 1 == 0));
            }
            zeros.push(*r0.last().unwrap());
            rank0.push(r0);
            let mut next: Vec<usize> = cur.iter().copied().filter(|v| (v >> bit) & 1 == 0).collect();
            next.extend(cur.iter().copied().filter(|v| (v >> bit) & 1 == 1));
            cur = next;
        }
        WaveletMatrix { bits, rank0, zeros }
    }
    fn kth_smallest(&self, mut l: usize, mut r: usize, mut k: usize) -> usize {
        let mut res = 0;
        for lvl in 0..self.bits {
            let bit = self.bits - 1 - lvl;
            let (zl, zr) = (self.rank0[lvl][l], self.rank0[lvl][r]);
            if k < zr - zl {
                l = zl;
                r = zr;
            } else {
                k -= zr - zl;
                res |= 1 << bit;
                l = self.zeros[lvl] + l - zl;
                r = self.zeros[lvl] + r - zr;
            }
        }
        res
    }
}

/// Range argmin of vals[i] (INF = excluded).
struct SparseMin {
    t: Vec<Vec<(i64, usize)>>,
}
impl SparseMin {
    fn new(vals: &[i64]) -> Self {
        let mut t = vec![vals.iter().copied().enumerate().map(|(i, v)| (v, i)).collect::<Vec<_>>()];
        let mut j = 1;
        while (1 << j) <= vals.len() {
            let p = &t[j - 1];
            let row = (0..=vals.len() - (1 << j)).map(|i| p[i].min(p[i + (1 << (j - 1))])).collect();
            t.push(row);
            j += 1;
        }
        SparseMin { t }
    }
    fn argmin(&self, l: usize, r: usize) -> usize {
        let j = (usize::BITS - (r - l).leading_zeros() - 1) as usize;
        self.t[j][l].min(self.t[j][r - (1 << j)]).1
    }
}

/// Point update, range max of (value, index).
struct MaxSeg {
    sz: usize,
    t: Vec<(i64, usize)>,
}
impl MaxSeg {
    fn new(n: usize) -> Self {
        let sz = n.max(1).next_power_of_two();
        MaxSeg { sz, t: vec![(NEG, 0); 2 * sz] }
    }
    fn set(&mut self, i: usize, v: i64) {
        let mut i = i + self.sz;
        self.t[i] = (v, i - self.sz);
        i >>= 1;
        while i > 0 {
            self.t[i] = self.t[2 * i].max(self.t[2 * i + 1]);
            i >>= 1;
        }
    }
    fn query(&self, l: usize, r: usize) -> (i64, usize) {
        let (mut res, mut l, mut r) = ((NEG, 0), l + self.sz, r + self.sz);
        while l < r {
            if l & 1 == 1 {
                res = res.max(self.t[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                res = res.max(self.t[r]);
            }
            l >>= 1;
            r >>= 1;
        }
        res
    }
}

/// Offline answers to (L, R, K): best size-K matching using only the
/// vertices L..R-1. O((n + q) log n).
struct RangeMatching {
    n: usize,
    x: Vec<i64>,
    pairs: Vec<(i64, usize, usize)>,
    wm: WaveletMatrix,
    pre: Vec<i64>,
    pmin: [SparseMin; 2],
    alive: Fenwick,  // 1 for every unmatched vertex
    dead_w: Fenwick, // weight of every matched vertex
    seg: MaxSeg,     // key i: x_i + x_next(i)
}

impl RangeMatching {
    fn new(x: &[i64]) -> Self {
        let n = x.len();
        let pairs = greedy_pairs(x);
        let mut time = vec![pairs.len(); n]; // step at which a vertex is matched
        for (t, &(_, i, j)) in pairs.iter().enumerate() {
            time[i] = t;
            time[j] = t;
        }
        let mut pre = vec![0];
        for &v in x {
            pre.push(pre.last().unwrap() + v);
        }
        let par = |p: usize| -> Vec<i64> { (0..n).map(|i| if i % 2 == p { x[i] } else { INF }).collect() };
        let mut rm = RangeMatching {
            n,
            x: x.to_vec(),
            wm: WaveletMatrix::new(&time),
            pairs,
            pre,
            pmin: [SparseMin::new(&par(0)), SparseMin::new(&par(1))],
            alive: Fenwick::new(n),
            dead_w: Fenwick::new(n),
            seg: MaxSeg::new(n),
        };
        for i in 0..n {
            rm.alive.add(i, 1);
        }
        for i in 0..n.saturating_sub(1) {
            rm.seg.set(i, x[i] + x[i + 1]);
        }
        rm
    }

    // the current global state: vertices with time < t are matched
    fn next(&self, i: usize) -> usize {
        // first unmatched index >= i, or n
        let c = self.alive.prefix(i);
        if c < self.alive.prefix(self.n) { self.alive.kth(c) } else { self.n }
    }
    fn prev(&self, i: usize) -> Option<usize> {
        // last unmatched index <= i
        let c = self.alive.prefix(i + 1);
        if c > 0 { Some(self.alive.kth(c - 1)) } else { None }
    }
    fn is_alive(&self, i: usize) -> bool {
        self.alive.prefix(i + 1) - self.alive.prefix(i) == 1
    }
    fn refresh(&mut self, i: usize) {
        let j = if self.is_alive(i) { self.next(i + 1) } else { self.n };
        let v = if j < self.n { self.x[i] + self.x[j] } else { NEG };
        self.seg.set(i, v);
    }
    fn flip(&mut self, i: usize, to_alive: bool) {
        let d = if to_alive { 1 } else { -1 };
        self.alive.add(i, d);
        self.dead_w.add(i, -d * self.x[i]);
        self.refresh(i);
        if i > 0 {
            if let Some(a) = self.prev(i - 1) {
                self.refresh(a);
            }
        }
    }

    fn answer_one(&mut self, l: usize, r: usize, k: usize) -> i64 {
        let left_cross = (l as i64 - self.alive.prefix(l)) & 1 == 1; // a matching edge crosses l
        let right_cross =
            ((self.n - r) as i64 - (self.alive.prefix(self.n) - self.alive.prefix(r))) & 1 == 1;
        let (u, v) = (self.next(l), self.prev(r - 1).unwrap()); // first / last unmatched in range
        let mut log: Vec<(usize, bool)> = Vec::new(); // (vertex, became unmatched?)
        if left_cross {
            log.push((self.pmin[l % 2].argmin(l, u), true)); // release b_L
        }
        if right_cross {
            log.push((self.pmin[(r - 1) % 2].argmin(v + 1, r), true)); // release b_R
        }
        for &(b, _) in &log {
            self.flip(b, true);
        }
        let s = (r - l) as i64 - (self.alive.prefix(r) - self.alive.prefix(l));
        assert!(s % 2 == 0 && s <= 2 * k as i64);
        for _ in 0..(2 * k as i64 - s) / 2 {
            // greedy on what is left
            let p = self.prev(r - 1).unwrap();
            let (_, i) = self.seg.query(l, p);
            let j = self.next(i + 1);
            self.flip(i, false);
            self.flip(j, false);
            log.push((i, false));
            log.push((j, false));
        }
        let ans = self.dead_w.prefix(r) - self.dead_w.prefix(l);
        for &(b, was_released) in log.iter().rev() {
            self.flip(b, !was_released); // undo everything
        }
        ans
    }

    fn answer(&mut self, queries: &[(usize, usize, usize)]) -> Vec<i64> {
        let mut res = vec![0; queries.len()];
        let mut by_time = vec![Vec::new(); self.pairs.len() + 1];
        for (idx, &(l, r, k)) in queries.iter().enumerate() {
            if k == 0 {
                res[idx] = 0;
            } else if 2 * k == r - l {
                res[idx] = self.pre[r] - self.pre[l];
            } else {
                // time at which the range gets its (2K)-th matched vertex
                by_time[self.wm.kth_smallest(l, r, 2 * k - 1)].push(idx);
            }
        }
        for t in 0..=self.pairs.len() {
            for &idx in &by_time[t] {
                let (l, r, k) = queries[idx];
                res[idx] = self.answer_one(l, r, k);
            }
            if t < self.pairs.len() {
                let (_, i, j) = self.pairs[t];
                self.flip(i, false);
                self.flip(j, false);
            }
        }
        res
    }
}
// PDFEND -- everything below is the test harness and is not printed in the PDF

// ------------------------------------------------------------ brute force
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
}

fn brute_f(x: &[i64]) -> Vec<i64> {
    let n = x.len();
    let mut dp = vec![vec![NEG; n / 2 + 1]; n + 1];
    dp[0][0] = 0;
    for i in 1..=n {
        for k in 0..=n / 2 {
            let mut best = dp[i - 1][k];
            if i >= 2 && k >= 1 && dp[i - 2][k - 1] > NEG {
                best = best.max(dp[i - 2][k - 1] + x[i - 1] + x[i - 2]);
            }
            dp[i][k] = best;
        }
    }
    dp[n].clone()
}

fn brute_nonadjacent(b: &[i64], k: usize) -> i64 {
    let n = b.len();
    let mut dp = vec![vec![NEG; k + 1]; n + 2];
    for row in dp.iter_mut() {
        row[0] = 0;
    }
    for i in 1..=n {
        for j in 1..=k {
            dp[i + 1][j] = dp[i][j];
            if dp[i - 1][j - 1] > NEG {
                dp[i + 1][j] = dp[i + 1][j].max(dp[i - 1][j - 1] + b[i - 1]);
            }
        }
    }
    dp[n + 1][k]
}

fn brute_subarrays(a: &[i64], k: usize) -> i64 {
    let mut dp = vec![[NEG, NEG]; k + 1];
    dp[0][0] = 0;
    for &v in a {
        let mut nd = vec![[NEG, NEG]; k + 1];
        for j in 0..=k {
            for inside in 0..2 {
                let cur = dp[j][inside];
                if cur == NEG {
                    continue;
                }
                nd[j][0] = nd[j][0].max(cur);
                if inside == 1 {
                    nd[j][1] = nd[j][1].max(cur + v);
                }
                if j < k {
                    nd[j + 1][1] = nd[j + 1][1].max(cur + v);
                }
            }
        }
        dp = nd;
    }
    dp[k][0].max(dp[k][1])
}

fn main() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let (mut heap_n, mut reduce_n, mut range_n) = (0, 0, 0);

    // Sections 4 and 5: the heap, the stack and the DP agree.
    for _ in 0..3000 {
        let n = rng.range(0, 14) as usize;
        let hi = [1, 2, 5, 50][rng.range(0, 3) as usize]; // tiny ranges force ties
        let x: Vec<i64> = (0..n).map(|_| rng.range(-hi, hi)).collect();
        let want = brute_f(&x);
        let gains: Vec<i64> = greedy_pairs(&x).iter().map(|p| p.0).collect();
        assert!(gains.windows(2).all(|w| w[0] >= w[1]), "{:?}", x);
        let mut f = vec![0];
        for g in &gains {
            f.push(f.last().unwrap() + g);
        }
        assert_eq!(f, want, "{:?}", x);
        let mut lin = gains_linear(&x);
        lin.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(lin, gains, "{:?}", x);
        for k in 0..=n / 2 {
            assert_eq!(sum_of_largest(gains_linear(&x), k), want[k]);
        }
        heap_n += 1;
    }

    // Section 3: the reductions.
    for _ in 0..2000 {
        let n = rng.range(1, 12) as usize;
        let hi = [2, 6, 40][rng.range(0, 2) as usize];
        let b: Vec<i64> = (0..n).map(|_| rng.range(-hi, hi)).collect();
        for k in 0..=(n + 1) / 2 {
            assert_eq!(max_nonadjacent(&b, k), brute_nonadjacent(&b, k), "{:?} {}", b, k);
        }
        let a: Vec<i64> = (0..n).map(|_| rng.range(-hi, hi)).collect();
        for k in 0..=n {
            assert_eq!(max_k_subarrays(&a, k), brute_subarrays(&a, k), "{:?} {}", a, k);
        }
        let p: Vec<i64> = (0..=n).map(|_| rng.range(0, hi)).collect();
        let d: Vec<i64> = (0..n).map(|i| p[i + 1] - p[i]).collect();
        for k in 0..=n {
            assert_eq!(max_profit_k_trades(&p, k), brute_subarrays(&d, k), "{:?} {}", p, k);
        }
        reduce_n += 1;
    }
    let w = [4, -7, 2, 9, 1];
    let x = edge_to_vertex_weights(&w);
    assert!((0..w.len()).all(|i| x[i] + x[i + 1] == w[i]));

    // Section 6: range queries against the brute force on the slice.
    for _ in 0..600 {
        let n = rng.range(2, 16) as usize;
        let hi = [1, 3, 9, 100][rng.range(0, 3) as usize];
        let x: Vec<i64> = (0..n).map(|_| rng.range(-hi, hi)).collect();
        let mut queries = Vec::new();
        for l in 0..n {
            for r in l + 1..=n {
                for k in 0..=(r - l) / 2 {
                    queries.push((l, r, k));
                }
            }
        }
        for i in (1..queries.len()).rev() {
            let j = rng.range(0, i as i64) as usize;
            queries.swap(i, j);
        }
        let got = RangeMatching::new(&x).answer(&queries);
        for (q, g) in queries.iter().zip(got) {
            assert_eq!(g, brute_f(&x[q.0..q.1])[q.2], "{:?} {:?}", x, q);
            range_n += 1;
        }
    }

    // A larger instance, just to see the structures cope.
    let n = 3000;
    let x: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000, 1_000_000)).collect();
    let queries: Vec<(usize, usize, usize)> = (0..300)
        .map(|_| {
            let l = rng.range(0, n as i64 - 2) as usize;
            let r = rng.range(l as i64 + 2, n as i64) as usize;
            (l, r, rng.range(1, ((r - l) / 2) as i64) as usize)
        })
        .collect();
    let got = RangeMatching::new(&x).answer(&queries);
    for _ in 0..5 {
        let i = rng.range(0, queries.len() as i64 - 1) as usize;
        let q = queries[i];
        assert_eq!(got[i], brute_f(&x[q.0..q.1])[q.2]);
    }
    println!(
        "all checks passed: heap/stack {}, reductions {}, range queries {}",
        heap_n, reduce_n, range_n
    );
}
