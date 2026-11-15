# Persistent Segment Tree

Difficulty: 7/9 · Category: Data Structures

Original notes on making a pointer-based segment tree persistent: path
copying (clone the O(log N) nodes on the update path instead of mutating
them, so both the old and new versions stay valid), extending the idea to
lazy range updates, and the same-range-copy operation across versions.
Covers the classic application — a chain of frequency-array versions
answering range count-<=x and, via simultaneous descent of two same-shaped
trees, range k-th-smallest in O(log N) — plus a note on extending it to
tree-path queries with inclusion-exclusion via LCA.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
cover both a branching multi-version array (point update, range min) and
the chained frequency-tree (range count-<=x, range k-th-smallest), the
latter cross-checked against brute force over every possible query.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
