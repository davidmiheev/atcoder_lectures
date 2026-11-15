# Zeta & Möbius Transform (Subset Relations)

Difficulty: 4/9 · Category: Set Operations & Functions

An original lecture on the zeta and Möbius transforms with respect to the
subset/superset relation on bitmasks — the subset zeta transform being
the "sum over subsets" (SOS) operation familiar from competitive
programming. Develops $O(4^n)$, $O(3^n)$, and $O(n\cdot 2^n)$
("fast zeta transform") algorithms for it and its inverse (the fast
Möbius transform), the superset mirror versions, and derives
$O(n\cdot 2^n)$ OR-convolution and AND-convolution algorithms by
composing a zeta transform, a pointwise product, and a Möbius transform.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement subset/superset zeta and Möbius transforms three ways each
(reproducing the lecture's own worked example exactly), cross-check the
three implementations against each other and confirm the Möbius
transform is an exact inverse of the zeta transform over 100 random
trials, and implement OR-/AND-convolution via the zeta/pointwise-product/
Möbius pipeline, cross-checked against an $O(4^n)$ brute-force
convolution over 100 further random trials each.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
