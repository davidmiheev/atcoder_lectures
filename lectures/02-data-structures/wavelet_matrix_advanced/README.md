# Wavelet Matrix (Advanced)

Difficulty: 7/9 · Category: Data Structures

Original notes extending the basic wavelet matrix (see the companion
`wavelet_matrix_basic` lecture) to: XOR-shifted range k-th smallest, range
sums over an attached array (and group-valued aggregates in general),
non-invertible aggregates like max via a segment tree per level (with
point updates via segment tree or Fenwick tree per level), binary search
against a running sum, and 2D rectangle counting via coordinate
compression, plus offline handling of points that are inserted and
deleted over time.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement and brute-force-cross-check three representative techniques:
XOR-shifted k-th smallest, generalized range sum via an attached array, and
2D rectangle counting by coordinate compression. The max/segment-tree,
point-update, binary-search, and offline-dynamic-points extensions are
covered in the notes with complexity analysis but left as an implementation
exercise, matching how the source material itself treats them.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
