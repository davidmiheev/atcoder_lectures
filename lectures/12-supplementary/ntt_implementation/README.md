# The Number Theoretic Transform

Difficulty: 6/9 · Category: Supplementary (original, not from AtCoder's index)

An original supplementary lecture — not drawn from AtCoder's Algorithm
Lectures index — implementing the Number Theoretic Transform from
scratch: a Discrete Fourier Transform carried out in Z/pZ (p = 998244353)
instead of C, giving the same O(n log n) fast-convolution trick as the
classical FFT with exact integer arithmetic. Added because several
lectures elsewhere in this series (the Kth-term and advanced Bostan–Mori
lectures) refer to NTT-based polynomial multiplication without building
the transform itself.

Covers primitive roots of unity modulo a prime, the forward/inverse
transform and convolution theorem, and the iterative in-place
Cooley–Tukey algorithm (bit-reversal permutation + butterfly passes).

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). Reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
implement the iterative NTT and NTT-based polynomial multiplication,
cross-checked against the naive O(nm) convolution on hundreds of random
polynomial pairs, plus a forward/inverse round-trip check and a
hand-verifiable worked example.
