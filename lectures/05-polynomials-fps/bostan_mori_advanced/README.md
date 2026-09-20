# Bostan–Mori's Algorithm (Advanced)

Difficulty: 8/9 · Category: Polynomials & Formal Power Series

Original notes extending the companion Kth-term lecture: the MSB-first
variant of Bostan–Mori that extracts a whole interval of coefficients
[x^[L,R]] P(x)/Q(x) in O(M(d) log R), derived from the same
Q(x)Q(-x)=V(x^2) identity but halving the interval's right end from the
top down instead of halving a single target index from the bottom up.
Covers why this is the transpose of ordinary (LSB-first) Bostan–Mori, a
byproduct O(M(d)) formal power series inverse, two applications (shifting
a linearly recurrent sequence's generating function, and computing
x^K mod Gamma(x) via interval extraction instead of binary
exponentiation), and — in prose, with complexity analysis — the
constant-factor speedup from keeping polynomials in DFT representation
across Bostan–Mori's iterations.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement interval coefficient extraction (cross-checked against
brute-force power series inversion, including negative L), the shift
application, and x^K mod Gamma(x) via Bostan–Mori (cross-checked against
Fiduccia's binary exponentiation) — all verified over hundreds of random
cases. The FFT/DFT-representation speedup is not implemented, since it
requires a full NTT engine and is a constant-factor refinement layered on
top of the core algorithm, consistent with how the source material itself
treats it.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
