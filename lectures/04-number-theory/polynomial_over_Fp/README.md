# Polynomials over F_p

Difficulty: 5/9 · Category: Number Theory

An original lecture on the algebra of polynomials with coefficients in
$\mathbb{F}_p$: polynomial long division, the remainder and factor
theorems, the bound that a nonzero degree-$n$ polynomial has at most $n$
roots, the resulting factorization $x^p - x = \prod_{a \in
\mathbb{F}_p}(x-a)$ and Wilson's theorem, and a worked comparison showing
exactly which of these fail once the modulus is composite (where
$\mathbb{F}_p$'s field structure — every nonzero element invertible — is
no longer available). Closes with the Schwartz–Zippel lemma, the
multivariate generalization used throughout randomized-algorithm
correctness proofs.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement polynomial long division from scratch over a small prime field
and cross-check the remainder theorem, factor theorem, and root-count
bound against brute-force evaluation over 300 random trials each; verify
the $x^p-x$ factorization and Wilson's theorem by direct expansion; and
verify Schwartz–Zippel by exhaustive evaluation of random bivariate
polynomials over the full field.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
