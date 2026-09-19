# Degree Sequences of Simple Undirected Graphs

Difficulty: 6/9 · Category: Graph Theory

Original notes on the Erdős–Gallai theorem (the necessary-and-sufficient
condition for a nonnegative integer sequence to be the degree sequence of
a simple graph) and the Havel–Hakimi algorithm (a greedy construction
realizing such a sequence, which doubles as a realizability test). Covers
the handshake lemma, a necessity proof via counting edges within/across a
prefix of vertices, a proof sketch of the Havel–Hakimi reduction via an
edge-rewiring exchange argument, and complexity (O(n) realizability
checking, and what it takes to get Havel–Hakimi down to O(n + sum d_i)).

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex), with reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
— both verified against the lecture's own three worked examples and
cross-checked exhaustively against every possible 5-vertex graph.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
