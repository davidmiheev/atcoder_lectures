# Prüfer Code and Cayley's Formula

Difficulty: 6/9 · Category: Combinatorics

Original notes on the Prüfer code, a bijection between labeled trees on
{0,…,n−1} and sequences in {0,…,n−1}^(n−2) obtained by repeatedly deleting
the smallest leaf and recording its neighbor. Covers the degree lemma (a
vertex of degree d occurs d−1 times in the code), the proof that encoding
is a bijection (by induction, with decoding read off from the proof),
heap-based O(n log n) and pointer-sweep O(n) implementations, and the
counting results that follow: Cayley's formula n^(n−2), trees with a
prescribed degree sequence and the degree generating function, the number
of ways to join a forest into a tree (n^(k−2) · ∏ sizes), rooted forests
with a prescribed root set (k·n^(n−k−1)), and spanning trees of complete
bipartite graphs (a^(b−1) · b^(a−1)) via a two-sequence variant of the code.

Every lemma and theorem in the notes is proved in full (including the
correctness of the O(n) pointer-sweep algorithm and the bipartite
two-sequence code), relying only on elementary tree facts that are
proved at the start.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
contain both encoders/decoders and cross-check them against brute-force
enumeration of all labeled trees for small n, against each other on random
larger trees, and against every counting formula in the notes (the Python
file also checks the bipartite two-sequence code).

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
