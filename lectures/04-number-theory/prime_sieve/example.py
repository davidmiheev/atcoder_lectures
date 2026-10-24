"""Prime sieves: Eratosthenes, the linear sieve (with smallest prime factors),
bulk multiplicative functions, and a segmented sieve for a high window.
"""

# PDFBEGIN
from __future__ import annotations

import math

def sieve_of_eratosthenes(n: int) -> list[int]:
    """Primes up to n, striking multiples of p starting at p*p."""
    if n < 2:
        return []
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0

    p = 2
    while p * p <= n:
        if is_prime[p]:
            # Slice assignment strikes the whole arithmetic progression at once.
            step_count = len(range(p * p, n + 1, p))
            is_prime[p * p : n + 1 : p] = bytearray(step_count)
        p += 1

    return [i for i in range(2, n + 1) if is_prime[i]]

def linear_sieve(n: int) -> tuple[list[int], list[int]]:
    """Every composite is written exactly once, by its smallest prime factor.

    Returns (primes, spf) with spf[0] = spf[1] = 0.
    """
    spf = [0] * (n + 1)
    primes: list[int] = []

    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i] = i
            primes.append(i)
        for p in primes:
            # The guard p <= spf[i] is what makes each composite unique.
            if p > spf[i] or p * i > n:
                break
            spf[p * i] = p

    return primes, spf

def factorize(n: int, spf: list[int]) -> list[tuple[int, int]]:
    """Factor n into (prime, exponent) pairs via repeated spf lookups."""
    factors: list[tuple[int, int]] = []
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        factors.append((p, e))
    return factors

def totient_table(n: int, spf: list[int]) -> list[int]:
    """Euler's totient for 1..n, using phi(p^e) = p^(e-1) (p - 1)."""
    phi = [0] * (n + 1)
    if n >= 1:
        phi[1] = 1
    for i in range(2, n + 1):
        p = spf[i]
        m, e = i, 0
        while m % p == 0:
            m //= p
            e += 1
        phi[i] = phi[m] * (p - 1) * p ** (e - 1)
    return phi

def divisor_count_table(n: int, spf: list[int]) -> list[int]:
    """Divisor count for 1..n, using d(p^e) = e + 1."""
    d = [0] * (n + 1)
    if n >= 1:
        d[1] = 1
    for i in range(2, n + 1):
        p = spf[i]
        m, e = i, 0
        while m % p == 0:
            m //= p
            e += 1
        d[i] = d[m] * (e + 1)
    return d

def segmented_sieve(lo: int, hi: int) -> list[int]:
    """Primes in [lo, hi], for hi far beyond what a full sieve could hold."""
    if hi < 2:
        return []
    lo = max(lo, 2)
    if lo > hi:
        return []

    base = sieve_of_eratosthenes(math.isqrt(hi))
    is_prime = bytearray([1]) * (hi - lo + 1)

    for p in base:
        # First multiple of p inside the window, but never below p*p.
        start = max(((lo + p - 1) // p) * p, p * p)
        for m in range(start, hi + 1, p):
            is_prime[m - lo] = 0

    return [lo + i for i, ok in enumerate(is_prime) if ok]

# PDFEND -- everything below is the test harness and is not printed in the PDF

def main() -> None:
    n = 50

    primes = sieve_of_eratosthenes(n)
    print(f"primes up to {n}: {primes}")

    lin_primes, spf = linear_sieve(n)
    print(f"linear sieve agrees: {lin_primes == primes}")
    print(f"spf[2..12] = {spf[2:13]}")

    for m in (360, 97, 1024):
        _, spf_m = linear_sieve(m)
        pretty = " * ".join(f"{p}^{e}" for p, e in factorize(m, spf_m))
        print(f"{m} = {pretty}")

    phi = totient_table(n, spf)
    print(f"phi(1..12) = {phi[1:13]}")

    d = divisor_count_table(n, spf)
    print(f"d(1..12)   = {d[1:13]}")

    lo, hi = 1_000_000_000_000, 1_000_000_000_100
    print(f"primes in [{lo}, {hi}]: {segmented_sieve(lo, hi)}")

if __name__ == "__main__":
    main()
