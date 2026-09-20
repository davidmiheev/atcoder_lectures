# Transposition Principle Introduction

Difficulty: 6/9 · Category: Linear Algebra

An original introductory lecture on the transposition principle
(Tellegen's principle): given an algorithm for $y=Ax$ built from
elementary steps $z_a\mathrel{+}=c\cdot z_b$, mechanically reversing the
step order and swapping each step's two roles (plus swapping load/store
positions) produces an algorithm for the *transposed* problem $X=A^TY$.
Proves this two ways (a matrix identity, and a duality-identity
invariant), then applies it to the classic prefix-sum solution of
range-sum queries to mechanically derive a range-add solution — cross-
checked against the well-known "imos method" — plus sketches further
transposition pairs (subset↔superset zeta transform, divisor↔multiple
zeta transform, point-update/range-sum↔range-update/point-query "dual"
Fenwick trees).

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement a generic interpreter for the elementary-step programs described
in the theorem, validate the transposition principle itself against a
directly-computed $A^TY$ over 300 random programs, reproduce the lecture's
own small worked example exactly, and confirm the mechanically-derived
range-add solution agrees with the lecture's hand-derived version, the
imos method, and brute force over 100 random trials — plus a bonus check
that transposing the fast subset zeta transform mechanically yields the
fast superset zeta transform.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
