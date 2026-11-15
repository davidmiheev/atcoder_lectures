# DP Acceleration via Eigenfunctions

Difficulty: 9/9 · Category: Polynomials & Formal Power Series

An original lecture on speeding up polynomial-valued DP recurrences
$A_{i+1} = f(A_i)$, where $f$ is a linear map built from differentiation
and multiplication by $x$. When $f$'s eigenfunctions are known in closed
form, $M$ iterations collapse to a single change of basis, one
exponentiation of each eigenvalue, and a change back — turning an $O(M)$
loop into $O(\log M)$. Works through three problems of increasing
difficulty: $A_{i+1}=(1+x)A_i'$ (eigenbasis $(1+x)^k$, solved via Taylor
shift), $A_{i+1}=A_i+A_i'$ (eigenbasis $e^{kx}$, with a truncation caveat
that needs $N+M$ padding), and $A_{i+1}=x(A_i+A_i')+b_iA_i$ (eigenbasis
$x^k e^{-x}$, an approach used in one solution to AtCoder ABC272 Ex).

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
fully implement and cross-check the first problem (eigenfunction method vs.
a direct $O(NM)$ DP simulation, 200 random trials) and numerically verify
the eigenfunction identities behind the other two — including a concrete
demonstration of the truncation caveat — via a from-scratch power series
$\exp$.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
