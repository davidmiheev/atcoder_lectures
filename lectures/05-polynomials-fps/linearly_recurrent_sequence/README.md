# Linearly Recurrent Sequences

Difficulty: 6/9 · Category: Polynomials & Formal Power Series

An original lecture on sequences satisfying a constant-coefficient linear
recurrence $A_i = c_1A_{i-1} + \cdots + c_dA_{i-d}$ (the Fibonacci
sequence being the standard example), and their two other
characterizations: a rational generating function, and a linear
combination of the entries of $M^i$ for a fixed "companion matrix" $M$.
Uses these to show the class is closed under sum, convolution, and
termwise product, and works through many examples — prefix sums, forced
recurrences, interleaved sequences, polynomial-valued and
polynomial-times-geometric sequences, chained/simultaneous recurrences,
and DP state-transition counting problems recast as matrix powers.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
generate sequences from explicit recurrences, build companion matrices and
confirm matrix powers reproduce them exactly (cross-checked against the
lecture's own worked values for $A_i=2A_{i-1}+3A_{i-2}$), and verify every
closure claim with a reusable technique: fit a candidate recurrence to a
prefix of a sequence by solving a linear system, then confirm it correctly
extrapolates dozens of further terms.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
