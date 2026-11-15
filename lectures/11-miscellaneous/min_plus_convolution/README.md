# Min-Plus Convolution of Convex Sequences

Difficulty: 4/9 · Category: Miscellaneous

Original notes on the min-plus convolution C_k = min_{i+j=k}(A_i + B_j),
generally O(NM) (and believed hard to beat), collapsing to O(N+M) when
both sequences are convex. Covers the "bag of items priced by the
difference sequence" reinterpretation that makes the fast algorithm (and
why the result stays convex) almost obvious, plus a DP-speedup
application: recognizing that a naive O(N^2) "choose exactly k elements
minimizing their sum" DP is secretly a chain of convex-convolutions with
(0, x), which collapses to plain sorting.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
match both of the lecture's own worked examples exactly, cross-check the
fast O(N+M) algorithm against the naive O(NM) definition on random convex
sequences (confirming the result stays convex), and verify the DP
matches sort-and-prefix-sum while staying convex at every step.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
