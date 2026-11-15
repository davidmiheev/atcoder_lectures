# Kth Term of a Linearly Recurrent Sequence

Difficulty: 6/9 · Category: Polynomials & Formal Power Series

Original notes on two O(M(d) log K)-time algorithms for the K-th term of
a degree-d linear recurrence, beating naive O(d^3 log K) matrix power:
**Fiduccia's algorithm** (reduces to computing x^K modulo the
recurrence's characteristic polynomial, with a cheap O(d) "multiply by x"
step) and **Bostan–Mori's algorithm** (works directly with the rational
generating function P(x)/Q(x), halving the target index each round via
the Q(x)Q(-x) even/odd-split trick). Covers the theorem linking recurrence
unrolling to polynomial remainders, the even-part lemma behind Bostan–
Mori, and the known constant-factor comparison between the two methods.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex), with reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
— both implement Fiduccia and Bostan–Mori with naive O(d^2) polynomial
arithmetic, reproduce the lecture's own worked example (A_5 = 16) exactly,
and cross-check against direct recurrence unrolling and against each
other, including at K = 10^12.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
