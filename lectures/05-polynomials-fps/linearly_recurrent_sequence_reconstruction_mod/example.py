"""Reeds-Sloane algorithm: recover a minimal-order linear recurrence for a
sequence over Z/mZ, for arbitrary (not necessarily prime) modulus m.

Strategy: factor m = prod p_j^{e_j}, solve the recovery problem modulo each
prime power separately (the hard part -- ordinary Berlekamp-Massey does not
work since Z/p^eZ is not a field), then recombine coefficient-wise with the
Chinese Remainder Theorem.
"""

# PDFBEGIN
from __future__ import annotations

def ext_gcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return a, 1, 0
    g, x1, y1 = ext_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1

def mod_inv(a: int, m: int) -> int:
    """Inverse of a unit `a` modulo `m` (requires gcd(a, m) = 1)."""
    g, x, _ = ext_gcd(a % m, m)
    assert g == 1, f"{a} is not invertible mod {m}"
    return x % m

def factorize(m: int) -> list[tuple[int, int]]:
    factors = []
    p = 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            factors.append((p, e))
        p += 1
    if m > 1:
        factors.append((m, 1))
    return factors

def valuation_unit(x: int, p: int, e: int, pe: int) -> tuple[int, int]:
    """Decompose a residue x in [0, p^e) as x = t * p^u, t a unit mod p^e.

    Returns the sentinel (e, 1) for x = 0 ("infinite valuation").
    """
    x %= pe
    if x == 0:
        return e, 1
    u = 0
    while x % p == 0:
        x //= p
        u += 1
    return u, x

def reeds_sloane_prime_power(s: list[int], p: int, e: int) -> tuple[int, list[int]]:
    """Reeds-Sloane recovery of the minimal (L, Q(x)) with Q(0) = 1, over
    Z/p^eZ. Returns (L, Q) where Q has L + 1 coefficients, Q[0] = 1.
    """
    pe = p**e
    s = [x % pe for x in s]
    n_terms = len(s)

    # Slot j tracks (Q_j, L_j) with Q_j(0) = p^j, the minimal-order
    # recurrence at "valuation level" j seen so far.
    q: list[list[int]] = [[pow(p, j, pe)] for j in range(e)]
    length = [0] * e
    # Auxiliary "last improvement" data per slot: (B, L_B, n_B, t_B) or None.
    aux: list[tuple[list[int], int, int, int] | None] = [None] * e

    for n in range(n_terms):
        delta = [0] * e
        for j in range(e):
            d = 0
            for i, qi in enumerate(q[j]):
                if i > n:
                    break
                d = (d + qi * s[n - i]) % pe
            delta[j] = d

        u = [0] * e
        t = [0] * e
        for j in range(e):
            u[j], t[j] = valuation_unit(delta[j], p, e, pe)

        q_new = [row[:] for row in q]
        l_new = length[:]

        for j in range(e):
            if u[j] == e:
                continue  # discrepancy is exactly 0: nothing to fix
            k = e - 1 - u[j]
            if aux[k] is None:
                # No donor available yet: record the failure, keep Q_j.
                l_new[j] = n + 1
            else:
                bk, lbk, nbk, tbk = aux[k]
                c = t[j] * mod_inv(tbk, pe) % pe
                shift = n - nbk
                l_new[j] = max(length[j], lbk + shift)
                poly = [0] * (l_new[j] + 1)
                for i, qi in enumerate(q[j]):
                    poly[i] = qi
                for i, bi in enumerate(bk):
                    poly[shift + i] = (poly[shift + i] - c * bi) % pe
                q_new[j] = poly
            if len(q_new[j]) != l_new[j] + 1:
                q_new[j] = q_new[j] + [0] * (l_new[j] + 1 - len(q_new[j]))

        # Save auxiliary data using the OLD (pre-update) donor slot k.
        aux_new = aux[:]
        for j in range(e):
            if length[j] < l_new[j]:
                k = e - 1 - u[j]
                aux_new[j] = (q[k][:], length[k], n, t[k])

        q, length, aux = q_new, l_new, aux_new

    return length[0], q[0]

def crt_pair(r1: int, m1: int, r2: int, m2: int) -> tuple[int, int]:
    g, x, _ = ext_gcd(m1, m2)
    assert g == 1, "moduli must be coprime"
    m = m1 * m2
    t = ((r2 - r1) % m2) * x % m2
    return (r1 + m1 * t) % m, m

def reeds_sloane(s: list[int], m: int) -> tuple[int, list[int]]:
    """Full Reeds-Sloane recovery for an arbitrary modulus m: solve modulo
    each prime power factor, then recombine every coefficient via CRT.
    """
    factors = factorize(m)
    per_factor = []
    for p, e in factors:
        l, q = reeds_sloane_prime_power(s, p, e)
        per_factor.append((p**e, l, q))

    l = max(lj for _, lj, _ in per_factor)

    q = [0] * (l + 1)
    for i in range(l + 1):
        r, cur_mod = 0, 1
        for pe, lj, qj in per_factor:
            coeff = qj[i] if i <= lj else 0
            r, cur_mod = crt_pair(r, cur_mod, coeff, pe)
        q[i] = r
    return l, q

def extrapolate(known: list[int], q: list[int], m: int, count: int) -> list[int]:
    """Extend a sequence using S_i = -sum_{j=1}^{L} Q[j] * S[i-j] (mod m),
    given Q[0] = 1. Returns `count` further terms after the known prefix.
    """
    s = known[:]
    l = len(q) - 1
    for _ in range(count):
        n = len(s)
        acc = 0
        for j in range(1, l + 1):
            if j <= n:
                acc = (acc + q[j] * s[n - j]) % m
        s.append((-acc) % m)
    return s[len(known):]

# PDFEND -- everything below is the test harness and is not printed in the PDF

def main() -> None:
    # --- Example 1: composite, non-prime-power modulus ---------------
    # Fibonacci-like sequence mod 12 = 4 * 3, satisfying S_i = S_{i-1} + S_{i-2}.
    m = 12
    fib = [1, 1]
    for i in range(2, 20):
        fib.append((fib[i - 1] + fib[i - 2]) % m)
    known = fib[:10]
    l, q = reeds_sloane(known, m)
    print(f"mod 12 Fibonacci-like: recovered L = {l}, Q = {q}")
    predicted = extrapolate(known, q, m, len(fib) - len(known))
    print(f"predicted tail: {predicted}")
    print(f"actual tail:    {fib[len(known):]}")
    print(f"match: {predicted == fib[len(known):]}\n")

    # --- Example 2: the classic non-uniqueness witness ----------------
    # S = (2,2,2,...) mod 4 satisfies both S_i = S_{i-1} and S_i = 3 S_{i-1}.
    s2 = [2] * 8
    l2, q2 = reeds_sloane(s2, 4)
    print(f"mod 4 constant sequence: recovered L = {l2}, Q = {q2}")
    more = extrapolate(s2, q2, 4, 4)
    print(f"predicted continuation: {more} (all should be 2)\n")

    # --- Example 3: prime modulus, sanity check against a field --------
    # Over Z/5Z this must reduce to ordinary Berlekamp-Massey behaviour.
    m3 = 5
    trib = [1, 1, 2]
    for i in range(3, 12):
        trib.append((trib[i - 1] + trib[i - 2] + trib[i - 3]) % m3)
    known3 = trib[:8]
    l3, q3 = reeds_sloane(known3, m3)
    print(f"mod 5 tribonacci-like: recovered L = {l3}, Q = {q3}")
    pred3 = extrapolate(known3, q3, m3, len(trib) - len(known3))
    print(f"match: {pred3 == trib[len(known3):]}")

if __name__ == "__main__":
    main()
