# Competitive Programming Lecture Notes

This repository is an **original**, English-language set of lecture notes on
competitive-programming topics, organized after the table of contents of
AtCoder's ["AtCoder Algorithm Lectures"](https://info.atcoder.jp/entry/algorithm_lectures/index)
series.

**Important:** this is not a translation of AtCoder's lecture text, and it is
not affiliated with or endorsed by AtCoder Inc. AtCoder's own lecture pages
(in Japanese, with video) remain the canonical source; each note here links
back to the matching entry in their index. Everything in this repo —
explanations, worked examples, and code — was written independently, from the
topic name and general algorithmic literature, not by translating AtCoder's
material.

Every lecture folder under [`lectures/`](lectures/) contains:

- `notes.tex` / `notes.pdf` — the write-up, typeset with the shared style in
  [`tex/lecturestyle.sty`](tex/lecturestyle.sty)
- `example.rs` — a self-contained Rust reference implementation
- `example.py` — a self-contained Python reference implementation (both
  files put the core algorithm first, then a `PDFEND` comment, then the
  randomized test harness; the PDF prints only the part between the
  `PDFBEGIN` and `PDFEND` marker comments)
- `README.md` — a one-paragraph summary and link to the source AtCoder entry

## Building

```
scripts/build_pdfs.sh     # compiles every notes.tex -> notes.pdf (needs latexmk + xelatex)
scripts/test_examples.sh  # compiles every example.rs and runs every example.py as a smoke test
```

## Index

### 1. Beginner
| Lecture | Difficulty |
|---|---|
| [Competitive Programming Basic Terminology](lectures/01-beginner/basic_terms/) | 1/9 |
| [Debugging with Random Tests](lectures/01-beginner/randomtest/) | 1/9 |

### 2. Data Structures
| Lecture | Difficulty |
|---|---|
| [Monoid Design for Range Queries](lectures/02-data-structures/range_query_monoid_design/) | 6/9 |
| [Sparse Table](lectures/02-data-structures/sparse_table/) | 3/9 |
| [Disjoint Sparse Table](lectures/02-data-structures/disjoint_sparse_table/) | 4/9 |
| [Sqrt Tree](lectures/02-data-structures/sqrt_tree/) | 5/9 |
| [Static Range Product Queries](lectures/02-data-structures/static_range_product/) | 7/9 |
| [Union-Find Complexity Bounds](lectures/02-data-structures/unionfind_complexity/) | 8/9 |
| [Wavelet Matrix (Basics)](lectures/02-data-structures/wavelet_matrix_basic/) | 6/9 |
| [Wavelet Matrix (Advanced)](lectures/02-data-structures/wavelet_matrix_advanced/) | 7/9 |
| [Persistent Segment Tree](lectures/02-data-structures/persistent_segment_tree/) | 7/9 |

### 3. Graph Theory
| Lecture | Difficulty |
|---|---|
| [Graph Theory Terminology](lectures/03-graph-theory/graph_terms/) | 1/9 |
| [Complete Graph Edge Coloring](lectures/03-graph-theory/complete_graph_edge_coloring/) | 3/9 |
| [Degree Sequences of Simple Undirected Graphs](lectures/03-graph-theory/degree_sequence/) | 6/9 |
| [Euler Tour Technique](lectures/03-graph-theory/euler_tour_technique/) | 4/9 |
| [Heavy-Light Decomposition](lectures/03-graph-theory/heavy_light_decomposition/) | 4/9 |
| [Linear Time LCA](lectures/03-graph-theory/linear_time_lca/) | 6/9 |
| [Tournament](lectures/03-graph-theory/tournament/) | 5/9 |

### 4. Number Theory
| Lecture | Difficulty |
|---|---|
| [Modular Arithmetic Basics](lectures/04-number-theory/modular_arithmetic/) | 1/9 |
| [Euclidean Algorithm](lectures/04-number-theory/euclidean_algorithm/) | 2/9 |
| [Prime Number Properties](lectures/04-number-theory/prime/) | 3/9 |
| [Prime Sieve](lectures/04-number-theory/prime_sieve/) | 3/9 |
| [Polynomials over F_p](lectures/04-number-theory/polynomial_over_Fp/) | 5/9 |
| [Primitive Root](lectures/04-number-theory/primitive_root/) | 6/9 |
| [Chinese Remainder Theorem](lectures/04-number-theory/chinese_remainder_theorem/) | 4/9 |
| [Prime-Related Complexity Analysis](lectures/04-number-theory/prime_related_complexity/) | 7/9 |

### 5. Polynomials & Formal Power Series
| Lecture | Difficulty |
|---|---|
| [Linearly Recurrent Sequences](lectures/05-polynomials-fps/linearly_recurrent_sequence/) | 6/9 |
| [Kth Term of a Linearly Recurrent Sequence](lectures/05-polynomials-fps/linearly_recurrent_sequence_kth_term/) | 6/9 |
| [Bostan–Mori Algorithm (Advanced)](lectures/05-polynomials-fps/bostan_mori_advanced/) | 8/9 |
| [Linear Recurrence Recovery](lectures/05-polynomials-fps/linearly_recurrent_sequence_reconstruction/) | 7/9 |
| [Linear Recurrence Recovery (mod m)](lectures/05-polynomials-fps/linearly_recurrent_sequence_reconstruction_mod/) | 9/9 |
| [Power Sum Formula (Faulhaber)](lectures/05-polynomials-fps/faulhaber_formula/) | 6/9 |
| [DP Acceleration via Eigenfunctions](lectures/05-polynomials-fps/eigenfunction/) | 9/9 |

### 6. Set Operations & Functions
| Lecture | Difficulty |
|---|---|
| [Zeta & Möbius Transform (Subset Relations)](lectures/06-set-operations/subset_zeta_mobius/) | 4/9 |
| [XOR Convolution (Division-Free)](lectures/06-set-operations/division_free_xor_convolution/) | 8/9 |
| [Sperner's Theorem](lectures/06-set-operations/sperner_theorem/) | 5/9 |

### 7. Linear Algebra
| Lecture | Difficulty |
|---|---|
| [Transposition Principle Introduction](lectures/07-linear-algebra/transposition_principle/) | 6/9 |

### 8. Combinatorics
| Lecture | Difficulty |
|---|---|
| [Parking Function](lectures/08-combinatorics/parking_function/) | 4/9 |
| [Cycle Lemma](lectures/08-combinatorics/cycle_lemma/) | 5/9 |
| [Hook Length Formula](lectures/08-combinatorics/hook_length_formula/) | 7/9 |

### 9. Strings
| Lecture | Difficulty |
|---|---|
| [String Rearrangement and Lexicographically Minimal Ordering](lectures/09-strings/lex_min_string_concatenation/) | 5/9 |

### 10. Game Theory
| Lecture | Difficulty |
|---|---|
| [Graph Games to Unvisited Vertices](lectures/10-game-theory/undirected_vertex_geography/) | 5/9 |

### 11. Miscellaneous
| Lecture | Difficulty |
|---|---|
| [Algebraic Structure Terminology](lectures/11-miscellaneous/algebra_terms/) | 6/9 |
| [Ackermann & Inverse Ackermann Functions](lectures/11-miscellaneous/ackermann_function/) | 7/9 |
| [Min-Plus Convolution of Convex Sequences](lectures/11-miscellaneous/min_plus_convolution/) | 4/9 |
| [Color-Coding](lectures/11-miscellaneous/color_coding/) | 5/9 |
| [Todd–Coxeter Algorithm](lectures/11-miscellaneous/todd_coxeter/) | 7/9 |

### 12. Supplementary (original, not from AtCoder's index)
| Lecture | Difficulty |
|---|---|
| [The Number Theoretic Transform](lectures/12-supplementary/ntt_implementation/) | 6/9 |

These lectures are not part of AtCoder's Algorithm Lectures table of
contents; they fill in prerequisite material that other lectures in this
series assume (e.g. an actual NTT implementation, referenced but not built
by the advanced Bostan–Mori lecture).

## Source

Topic list and difficulty ratings from the AtCoder Algorithm Lectures index:
<https://info.atcoder.jp/entry/algorithm_lectures/index>
