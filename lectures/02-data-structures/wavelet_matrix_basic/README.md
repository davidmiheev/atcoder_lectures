# Wavelet Matrix (Basics)

Difficulty: 6/9 · Category: Data Structures

Original notes on the wavelet matrix: how a height-H binary trie over a
static array's values can be queried against *any* index range [L, R) in
O(H) time, by flattening the (conceptual) wavelet tree level-by-level into
one bit-plane and rank table per height. Covers the O(1) child-range
formula, O(N H) construction via MSB-first stable partitioning, and the
two flagship queries (range k-th smallest, range count in [a, b)) built on
top of it.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex), with reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
— both are cross-checked against brute force over every possible
`(l, r, k)` and `(l, r, a, b)` query on a 9-element example array.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
