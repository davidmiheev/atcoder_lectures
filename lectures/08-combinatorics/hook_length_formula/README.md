# Hook Length Formula

Difficulty: 7/9 · Category: Combinatorics

An original lecture on the hook length formula: the number of standard Young
tableaux of shape λ ⊢ n is n! divided by the product of the hook lengths of the
cells. Two independent proofs are given. The first is the probabilistic hook
walk of Greene–Nijenhuis–Wilf, built on a weighted lattice-path identity; it
also yields a uniform random tableau generator. The second reads a tableau as
non-colliding walkers on a line, turns the count into a determinant with a
sign-reversing involution, evaluates it with a Vandermonde determinant
(Frobenius' formula), and converts factorials to hook lengths with a
boundary-path lemma. Figures illustrate hooks, a hook walk, the walker
cancellation and the boundary path.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
compute hook lengths, the exact count and the count modulo 998244353,
Frobenius' product formula and a hook-walk tableau sampler. The full files
check them against the corner recurrence, a numerically evaluated determinant
and brute-force enumeration (all partitions up to n = 12), verify Σ (f^λ)² = n!,
the Catalan special case, the row-hook lemma, the hook-walk corner
probabilities, and the uniformity of the sampler.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
