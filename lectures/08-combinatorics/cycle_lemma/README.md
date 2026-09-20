# Cycle Lemma

Difficulty: 5/9 · Category: Combinatorics

Original notes on the Cycle Lemma: for an integer sequence with every
term at most 1 and nonnegative sum K, exactly K of its cyclic shifts have
every prefix sum (after the first) strictly positive. Covers the proof
idea (rotating the prefix-sum path to its minimum, then cutting at the
last occurrence of each height) and derives four classical counting
results from it via the same "prepend +1" trick: Catalan numbers (two
equivalent closed forms), the Catalan convolution formula, Fuss–Catalan
numbers, and Narayana numbers.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
verify the core lemma exactly against the lecture's own worked example
(and against random sequences), then brute-force-check every corollary
against its closed form for small parameters.

Not a translation — topic and technical content taken from the AtCoder
Algorithm Lectures index: <https://info.atcoder.jp/entry/algorithm_lectures/index>.
