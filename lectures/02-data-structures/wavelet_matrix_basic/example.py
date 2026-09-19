"""Wavelet Matrix: answers range k-th-smallest and range count-in-[a,b)
queries on a static array of small non-negative integers, by re-encoding
the array as one bit-plane per level plus a rank (prefix-count) table.
"""
from __future__ import annotations


class WaveletMatrix:
    def __init__(self, height: int, a: list[int]) -> None:
        self.height = height
        self.n = n = len(a)
        a = a[:]
        # data[h][i] = number of 1-bits among the first i elements of the
        # level-h bit-plane (a "rank" prefix-sum table).
        self.data: list[list[int]] = [[] for _ in range(height)]

        # Process levels from the most significant bit down. At each level,
        # stably partition the current order into "bit = 0" then "bit = 1",
        # and record the bit-plane's prefix sum before overwriting `a`.
        for h in range(height - 1, -1, -1):
            dir_bits = [(a[i] >> h) & 1 for i in range(n)]
            left = [a[i] for i in range(n) if dir_bits[i] == 0]
            right = [a[i] for i in range(n) if dir_bits[i] == 1]
            a = left + right

            prefix = [0] * (n + 1)
            for i in range(n):
                prefix[i + 1] = prefix[i] + dir_bits[i]
            self.data[h] = prefix

    def _subtree_range(self, h: int, l: int, r: int) -> tuple[int, int, int, int]:
        """Given a node at height h+1 covering [l, r), return the ranges its
        left and right children cover at height h."""
        dat = self.data[h]
        a1, b1 = dat[l], dat[r]
        a0, b0 = l - a1, r - b1
        c0 = self.n - dat[-1]  # total zeros at this level
        return a0, b0, c0 + a1, c0 + b1

    def kth_smallest(self, l: int, r: int, k: int) -> int:
        """0-indexed k-th smallest value among A[l:r]."""
        value = 0
        for h in range(self.height - 1, -1, -1):
            l0, r0, l1, r1 = self._subtree_range(h, l, r)
            left_size = r0 - l0
            if k < left_size:
                l, r = l0, r0
            else:
                k -= left_size
                value |= 1 << h
                l, r = l1, r1
        return value

    def count_less(self, l: int, r: int, x: int) -> int:
        """Count of A[l:r] strictly less than x."""
        limit = 1 << self.height
        if x <= 0:
            return 0
        if x >= limit:
            return r - l
        count = 0
        for h in range(self.height - 1, -1, -1):
            bit = (x >> h) & 1
            l0, r0, l1, r1 = self._subtree_range(h, l, r)
            if bit == 1:
                count += r0 - l0  # everything landing in the 0-branch is < x
                l, r = l1, r1
            else:
                l, r = l0, r0
        return count

    def range_count(self, l: int, r: int, a: int, b: int) -> int:
        """Count of A[l:r] with a <= value < b."""
        return self.count_less(l, r, b) - self.count_less(l, r, a)


def main() -> None:
    a = [3, 1, 2, 1, 1, 0, 3, 2, 2]
    n = len(a)
    height = 2  # values fit in 0..4

    wm = WaveletMatrix(height, a)

    # Cross-check against brute force on every (l, r, k) and (l, r, a, b).
    all_ok = True
    for l in range(n):
        for r in range(l + 1, n + 1):
            sorted_slice = sorted(a[l:r])
            for k in range(len(sorted_slice)):
                got = wm.kth_smallest(l, r, k)
                if got != sorted_slice[k]:
                    all_ok = False
                    print(f"MISMATCH kth_smallest({l},{r},{k}): got {got}, want {sorted_slice[k]}")
            for lo in range(5):
                for hi in range(lo, 5):
                    expected = sum(1 for v in a[l:r] if lo <= v < hi)
                    got = wm.range_count(l, r, lo, hi)
                    if got != expected:
                        all_ok = False
                        print(f"MISMATCH range_count({l},{r},{lo},{hi}): got {got}, want {expected}")

    print(f"brute-force cross-check over all queries: {'PASS' if all_ok else 'FAIL'}")

    print(f"\nA = {a}")
    print(f"kth_smallest(2, 7, 0) = {wm.kth_smallest(2, 7, 0)}  (median-ish of A[2:7])")
    print(f"kth_smallest(2, 7, 4) = {wm.kth_smallest(2, 7, 4)}  (largest of A[2:7])")
    print(f"range_count(0, 9, 1, 3) = {wm.range_count(0, 9, 1, 3)}  (values in [1,3) over all of A)")


if __name__ == "__main__":
    main()
