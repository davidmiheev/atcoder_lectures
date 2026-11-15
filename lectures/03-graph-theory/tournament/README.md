# Tournaments

Difficulty: 5/9 · Category: Graph Theory

Original notes on tournaments (complete directed graphs from round-robin
results): the equivalence of transitivity, being a DAG, and having no
3-cycle; why strongly connected components form a strict total order and
can be decomposed from in-degrees alone; Rédei's theorem guaranteeing a
Hamiltonian path (and its link to comparison sorting); Camion's theorem
guaranteeing a Hamiltonian cycle when strongly connected; and Landau's
theorem on 2-kings.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
cover all five results — transitivity testing, in-degree-only SCC
decomposition, Hamiltonian path/cycle construction, and 2-king — each
validated against hundreds of random tournaments (the SCC method is also
cross-checked against a reachability-based reference implementation).

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
