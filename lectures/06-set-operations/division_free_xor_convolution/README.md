# XOR Convolution (Division-Free)

Difficulty: 8/9 · Category: Set Operations & Functions

An original lecture on computing the XOR convolution of two length-$2^N$
arrays using only ring operations (add, subtract, multiply) — no
division, unlike the classical Walsh–Hadamard-transform solution, which
divides by $2$ or a power of it. Originally posed as
[UOJ Long Round 3, Problem D](https://uoj.ac/contest/103/problem/1015).
The technique: recast XOR convolution as multilinear polynomial
multiplication modulo $x_i^2-1$, change basis via $x_i=y_i-1$ to turn
this into multiplication modulo $y_i^2-2y_i$, and solve that with a
"ranked OR-convolution" — exactly the subset-convolution technique, kept
as a full polynomial in an extra set-size-tracking variable so the final
$2^{|S|+|T|-|U|}$ factor needs only repeated multiplication by $2$, never
division. Overall $O(N^2\cdot 2^N)$.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement every stage separately and cross-check it: the basis change is
checked to round-trip exactly, the ranked OR-convolution routine is
checked against an $O(4^N)$ brute-force evaluation, and the full pipeline
is checked against an $O(4^N)$ brute-force XOR convolution — all
arithmetic carried out modulo $2^{20}$ specifically because $2$ is a zero
divisor there (so the Walsh–Hadamard approach is simply unavailable),
demonstrating the algorithm works where the classical one cannot.

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
