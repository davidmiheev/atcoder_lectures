# Sperner's Theorem

Difficulty: 5/9 · Category: Set Operations & Functions

Original notes on Sperner's theorem: the maximum antichain in the Boolean
lattice B_n has exactly C(n, floor(n/2)) elements. Gives two independent
proofs — the LYM inequality (a probabilistic argument over random maximal
chains, with the exact equality-case characterization) and an explicit
symmetric chain decomposition (which also yields a maximum matching
between consecutive levels of the lattice) — plus a note on the
generalization to products of chains.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
build and verify the symmetric chain decomposition and level matching for
n up to 8–10, check the LYM inequality (with exact rational arithmetic)
on hundreds of random maximal antichains including the exact equality
case, and cross-check the whole theorem independently via Dilworth's
theorem (max antichain = min chain cover, computed by bipartite matching).

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
