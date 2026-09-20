# Parking Function

Difficulty: 4/9 · Category: Combinatorics

Original notes on parking functions: a length-n sequence with values in
[1,n] where n cars, each taking the first free spot at or after their
preference, all successfully park. Covers four equivalent
characterizations (direct simulation, a threshold count condition, a
sorted-sequence bound, and a Hall's-theorem-style matching condition) and
the count (n+1)^(n-1), proved via a circular extension of the parking lot
where every sequence succeeds and, by rotational symmetry, every spot is
left empty equally often — plus the generalization to shorter sequences.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
reproduce the lecture's own worked examples exactly, cross-check the
simulation against the threshold condition and the circular-extension
bijection by brute force through n=6, verify the slot-symmetry claim
underlying the counting proof, and confirm the generalized count formula.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
