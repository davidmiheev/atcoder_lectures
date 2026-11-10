# Matchings on a Path and Their Applications

Difficulty: 8/9 · Category: Graph Theory

Original notes on weighted matchings on a path graph and the optimization
problems that reduce to them: choosing exactly K pairwise non-adjacent
array entries, choosing K disjoint subarrays with maximum total, and
making exactly K stock trades. Covers the reductions (including why K
trades correspond to matchings of size N−K on a 2N-vertex buy/sell path
and why profit(K) is minus the sum of the K smallest gains), a
self-contained proof that augmenting along the best augmenting path yields
optimal matchings of every size, the resulting O(n log n) heap greedy,
concavity of the optimum in K, the peeling theorem (any local maximum
among neighbouring candidate pairs may be taken first), a stack-based O(n)
algorithm with selection for a single K, and, as an advanced topic,
offline static range queries in O((n+q) log n) using the global greedy's
deletion times, a "release lemma" bounding what must be undone at the
interval boundary, a wavelet matrix, Fenwick trees, sparse tables and a
segment tree.

Every lemma and theorem is proved in full, with the release lemma stated
and proved for arbitrary weights (no distinctness assumption).

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
contain the heap greedy, the stack algorithm, the three reductions and the
range-query structure, and cross-check them against dynamic programming on
thousands of random inputs (including tie-heavy ones with weights in
{−1,0,1}) and on about 10^5 exhaustively enumerated range queries.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
