# Todd–Coxeter Algorithm

Difficulty: 7/9 · Category: Miscellaneous

An original lecture on the Todd–Coxeter algorithm for monoids. It recalls
generators and right Cayley graphs (the transition table of a monoid),
explains how a finite list of relations between words defines a monoid, and
presents the algorithm that builds the Cayley graph of that monoid: grow a
graph of provisional vertices, create a vertex wherever an edge is missing,
and merge (with union-find) vertices that must be the same element. Proofs
cover correctness when the algorithm stops and termination whenever the
monoid is finite; the notes also discuss cost, the dependence on the order of
work, infinite monoids (the run needs a cap) and the word problem. A full
trace of a small example and three figures illustrate the algorithm.

See [`notes.pdf`](notes.pdf) / [`notes.tex`](notes.tex). The reference
implementations in [`example.rs`](example.rs) and [`example.py`](example.py)
return the transition table with the identity numbered 0 and the remaining
vertices numbered breadth-first, and take a vertex cap. The full files check
them on presentations of S3, the dihedral groups D3–D8 and the Coxeter
presentations of S4 and S5 (each element is mapped to a permutation and the
map is verified to be a bijection onto the permutation group), on a monoid that
is not a group, on an infinite monoid, and on 200 random presentations with no
known answer (relations hold, multiplication is associative, rewriting is
respected).

This is not a translation of AtCoder's lecture text; see the repository
[README](../../../README.md) for the project's scope and sourcing policy.
