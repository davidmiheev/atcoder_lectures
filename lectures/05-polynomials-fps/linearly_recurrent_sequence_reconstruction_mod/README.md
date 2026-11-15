# Linear Recurrence Recovery (mod m)

Difficulty: 9/9 · Category: Polynomials & Formal Power Series

Original notes on the Reeds–Sloane algorithm: recovering a minimal-order
linear recurrence for a sequence over `Z/mZ` for an arbitrary (not
necessarily prime) modulus, by reducing to each prime-power factor of `m`
via CRT and, within a prime power, tracking one candidate recurrence per
p-adic valuation level so that every cancellation step only ever needs to
divide by a unit. Covers why plain Berlekamp–Massey breaks over a
composite modulus, the non-uniqueness of the minimal recurrence (with a
worked counterexample mod 4), and why the predicted continuation is still
unique whenever `2L <= N`.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex), with reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
— both implement the full factor → recover → CRT-recombine pipeline and are
checked against a mod-12 Fibonacci-like sequence, the classic mod-4
non-uniqueness witness, and a mod-5 (prime) sanity check.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>,
and from the original paper: J. A. Reeds, N. J. A. Sloane, "Shift-Register
Synthesis (Modulo m)", SIAM Journal on Computing 14(3), 1985.
