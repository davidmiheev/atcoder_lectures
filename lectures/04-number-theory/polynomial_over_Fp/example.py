"""Polynomials over F_p: the division/remainder/factor theorems that hold
because F_p is a field (every nonzero element is invertible), the bound
"a degree-n polynomial has at most n roots," the resulting factorization
x^p - x = prod_{a in F_p} (x - a) (hence Wilson's theorem), and the
Schwartz-Zippel lemma bounding how often a nonzero multivariate
polynomial vanishes at a random point.

This file verifies all of these computationally: polynomial long division
against a brute-force check, the root-count bound against naive
root-finding by evaluation, x^p - x's factorization by directly expanding
prod (x-a) and comparing coefficients, Wilson's theorem by direct
computation, and Schwartz-Zippel by measuring vanishing frequency of
random nonzero polynomials against the d/p bound.
"""

# PDFBEGIN
from __future__ import annotations

import random

P = 17  # small prime, so exhaustive checks over all of F_p are cheap

Poly = list[int]

def add(a: int, b: int) -> int:
    return (a + b) % P

def sub(a: int, b: int) -> int:
    return (a - b) % P

def mul(a: int, b: int) -> int:
    return (a % P) * (b % P) % P

def mod_pow(base: int, exp: int, m: int) -> int:
    return pow(base % m, exp, m)

def mod_inv(a: int) -> int:
    return mod_pow(a, P - 2, P)

def trim(a: Poly) -> Poly:
    a = a[:]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def degree(a: Poly) -> int:
    a = trim(a)
    if len(a) == 1 and a[0] == 0:
        return -1
    return len(a) - 1

def eval_poly(a: Poly, x: int) -> int:
    acc = 0
    for c in reversed(a):
        acc = add(mul(acc, x), c)
    return acc

def poly_divmod(f: Poly, g: Poly) -> tuple[Poly, Poly]:
    """Polynomial long division: returns (quotient, remainder) with
    deg(remainder) < deg(divisor). Requires g != 0."""
    assert degree(g) >= 0, "division by the zero polynomial"
    r = trim(f)
    dg = degree(g)
    lead_inv = mod_inv(g[-1])
    q = [0]
    while degree(r) >= dg:
        dr = degree(r)
        shift = dr - dg
        coeff = mul(r[-1], lead_inv)
        if len(q) < shift + 1:
            q += [0] * (shift + 1 - len(q))
        q[shift] = coeff
        for i, gi in enumerate(g):
            r[shift + i] = sub(r[shift + i], mul(coeff, gi))
        r = trim(r)
        if degree(r) < 0:
            break
    return trim(q), trim(r)

def poly_mul(a: Poly, b: Poly) -> Poly:
    if degree(a) < 0 or degree(b) < 0:
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = add(out[i + j], mul(ai, bj))
    return trim(out)

def roots_by_evaluation(f: Poly) -> list[int]:
    """All roots of f in F_p, found by brute-force evaluation (P is small)."""
    return [a for a in range(P) if eval_poly(f, a) == 0]

# PDFEND -- everything below is the test harness and is not printed in the PDF

def main() -> None:
    rng = random.Random(20260320)

    # --- Remainder theorem: f(x) mod (x - a) == f(a), for many random f, a ---
    remainder_ok = True
    for _ in range(300):
        deg = rng.randint(0, 6)
        f = trim([rng.randint(0, P - 1) for _ in range(deg + 1)])
        a = rng.randint(0, P - 1)
        divisor = [sub(0, a), 1]  # x - a
        _, r = poly_divmod(f, divisor)
        r_const = 0 if degree(r) < 0 else r[0]
        if r_const != eval_poly(f, a):
            remainder_ok = False
            print(f"REMAINDER MISMATCH f={f} a={a}: r={r_const} f(a)={eval_poly(f, a)}")
    print(f"Remainder theorem (f(x) mod (x-a) == f(a)) over 300 random trials: {'PASS' if remainder_ok else 'FAIL'}")

    # --- Factor theorem + root count bound ---
    factor_ok = True
    root_bound_ok = True
    for _ in range(300):
        deg = rng.randint(1, 6)
        f = [rng.randint(0, P - 1) for _ in range(deg + 1)]
        if f[deg] == 0:
            f[deg] = 1  # ensure the intended degree is exact
        roots = roots_by_evaluation(f)
        if len(roots) > degree(f):
            root_bound_ok = False
            print(f"ROOT COUNT BOUND VIOLATED f={f} roots={roots}")
        prod: Poly = [1]
        for r in roots:
            prod = poly_mul(prod, [sub(0, r), 1])
        if degree(prod) >= 0:
            q, rem = poly_divmod(f, prod)
            if degree(rem) >= 0:
                factor_ok = False
                print(f"FACTOR THEOREM: prod(x-root) does not divide f; f={f} roots={roots}")
            reexpanded = poly_mul(q, prod)
            if reexpanded != trim(f):
                factor_ok = False
                print(f"FACTOR THEOREM: re-expansion mismatch f={f}")
    print(f"Root count bound (#roots <= degree) over 300 random trials: {'PASS' if root_bound_ok else 'FAIL'}")
    print(f"Factor theorem (prod(x-root_i) divides f, re-expands to f) over 300 random trials: {'PASS' if factor_ok else 'FAIL'}")

    # --- x^p - x = prod_{a in F_p} (x - a), checked by direct expansion ---
    rhs: Poly = [1]
    for a in range(P):
        rhs = poly_mul(rhs, [sub(0, a), 1])
    lhs = [0] * (P + 1)
    lhs[P] = 1
    lhs[1] = sub(lhs[1], 1)  # x^p - x
    xp_x_ok = trim(lhs) == trim(rhs)
    print(f"x^p - x == prod_{{a in F_p}}(x-a), checked by direct expansion (p={P}): {'PASS' if xp_x_ok else 'FAIL'}")

    # --- Wilson's theorem ---
    factorial = 1
    for k in range(1, P):
        factorial = mul(factorial, k)
    coeff_x1 = rhs[1] if len(rhs) > 1 else 0
    wilson_ok = factorial == P - 1 and coeff_x1 == P - 1
    print(f"Wilson's theorem ((p-1)! == -1 mod p, matching the x^1 coefficient above): {'PASS' if wilson_ok else 'FAIL'}")

    # --- Schwartz-Zippel lemma ---
    sz_ok = True
    for _ in range(30):
        d = rng.randint(1, 5)
        terms: list[tuple[int, int, int]] = []
        for i in range(d + 1):
            for j in range(d - i + 1):
                c = rng.randint(0, P - 1)
                if c != 0:
                    terms.append((i, j, c))
        if not terms:
            continue

        def eval_f(x: int, y: int) -> int:
            acc = 0
            for i, j, c in terms:
                acc = add(acc, mul(c, mul(mod_pow(x, i, P), mod_pow(y, j, P))))
            return acc

        zero_count = 0
        total = P * P
        for x in range(P):
            for y in range(P):
                if eval_f(x, y) == 0:
                    zero_count += 1
        observed = zero_count / total
        bound = d / P
        if observed > bound + 1e-9:
            sz_ok = False
            print(f"SCHWARTZ-ZIPPEL VIOLATED d={d}: observed={observed:.4f} bound={bound:.4f}")
    print(f"Schwartz-Zippel lemma (P[f=0] <= d/p, checked exhaustively over F_p x F_p): {'PASS' if sz_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
