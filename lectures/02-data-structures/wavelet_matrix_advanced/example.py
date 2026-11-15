"""Wavelet Matrix, advanced queries: XOR-shifted k-th smallest, range sums
over an attached value array, and 2D rectangle counting via coordinate
compression. Builds on the plain range k-th-smallest / range-count wavelet
matrix from the companion "basics" lecture.
"""

# PDFBEGIN
from __future__ import annotations

import bisect

class WaveletMatrixSum:
    def __init__(self, height: int, a: list[int], x: list[int]) -> None:
        self.height = height
        self.n = n = len(a)
        a = a[:]
        x = x[:]

        # Plain range sum of X over the ORIGINAL (unpermuted) index order.
        # Used to answer threshold >= 2^height, where every element qualifies.
        self.total_prefix = [0] * (n + 1)
        for i in range(n):
            self.total_prefix[i + 1] = self.total_prefix[i] + x[i]

        # dat[h][i]  = number of 1-bits among the first i elements of the
        #              level-h bit-plane (same rank table as the basic matrix).
        # sum0[h][i] = sum of X over the first i elements of level h's
        #              *pre-partition* order, restricted to those with bit 0.
        self.dat: list[list[int]] = [[] for _ in range(height)]
        self.sum0: list[list[int]] = [[] for _ in range(height)]

        for h in range(height - 1, -1, -1):
            dir_bits = [(a[i] >> h) & 1 for i in range(n)]
            dat_h = [0] * (n + 1)
            sum0_h = [0] * (n + 1)
            for i in range(n):
                dat_h[i + 1] = dat_h[i] + dir_bits[i]
                sum0_h[i + 1] = sum0_h[i] + (0 if dir_bits[i] else x[i])
            self.dat[h] = dat_h
            self.sum0[h] = sum0_h

            left_a = [a[i] for i in range(n) if dir_bits[i] == 0]
            right_a = [a[i] for i in range(n) if dir_bits[i] == 1]
            left_x = [x[i] for i in range(n) if dir_bits[i] == 0]
            right_x = [x[i] for i in range(n) if dir_bits[i] == 1]
            a = left_a + right_a
            x = left_x + right_x

    def _subtree_range(self, h: int, l: int, r: int) -> tuple[int, int, int, int]:
        dat = self.dat[h]
        a1, b1 = dat[l], dat[r]
        a0, b0 = l - a1, r - b1
        c0 = self.n - dat[-1]
        return a0, b0, c0 + a1, c0 + b1

    def kth_smallest(self, l: int, r: int, k: int) -> int:
        """Plain range k-th smallest (as in the basics lecture)."""
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

    def kth_smallest_xor(self, l: int, r: int, k: int, mask: int) -> int:
        """k-th smallest of (A_i XOR mask) over i in [l, r). Bit h of
        A_i XOR mask equals bit h of A_i, flipped whenever mask's bit h is
        1 -- so we just swap which child counts as "left" at levels where
        mask has a 1 bit.
        """
        value = 0
        for h in range(self.height - 1, -1, -1):
            l0, r0, l1, r1 = self._subtree_range(h, l, r)
            if (mask >> h) & 1:
                l0, r0, l1, r1 = l1, r1, l0, r0
            left_size = r0 - l0
            if k < left_size:
                l, r = l0, r0
            else:
                k -= left_size
                value |= 1 << h
                l, r = l1, r1
        return value

    def sum_less(self, l: int, r: int, threshold: int) -> int:
        """Sum of X_i over i in [l, r) with A_i < threshold."""
        limit = 1 << self.height
        if threshold <= 0:
            return 0
        if threshold >= limit:
            return self.total_prefix[r] - self.total_prefix[l]
        total = 0
        for h in range(self.height - 1, -1, -1):
            bit = (threshold >> h) & 1
            l0, r0, l1, r1 = self._subtree_range(h, l, r)
            if bit:
                total += self.sum0[h][r] - self.sum0[h][l]
                l, r = l1, r1
            else:
                l, r = l0, r0
        return total

    def range_sum(self, l: int, r: int, a: int, b: int) -> int:
        """Sum of X_i over i in [l, r) with a <= A_i < b."""
        return self.sum_less(l, r, b) - self.sum_less(l, r, a)

    def range_count(self, l: int, r: int, a: int, b: int) -> int:
        """Count of i in [l, r) with a <= A_i < b (range_sum with X = all ones)."""
        return self.range_sum(l, r, a, b)

class RectangleCounter:
    """2D rectangle counting: given N points, answer "how many points have
    L <= x < R and D <= y < U" by coordinate-compressing y, sorting by x,
    and reusing the 1D wavelet matrix range-count query.
    """

    def __init__(self, points: list[tuple[int, int]]) -> None:
        pts = sorted(points, key=lambda p: p[0])
        self.sorted_x = [p[0] for p in pts]

        self.sorted_unique_y = sorted(set(p[1] for p in pts))
        ranks = [bisect.bisect_left(self.sorted_unique_y, p[1]) for p in pts]

        n = len(pts)
        height = 1
        while (1 << height) < max(len(self.sorted_unique_y), 1):
            height += 1
        ones = [1] * n
        self.wm = WaveletMatrixSum(height, ranks, ones)

    def query(self, l: int, r: int, d: int, u: int) -> int:
        il = bisect.bisect_left(self.sorted_x, l)
        ir = bisect.bisect_left(self.sorted_x, r)
        dl = bisect.bisect_left(self.sorted_unique_y, d)
        du = bisect.bisect_left(self.sorted_unique_y, u)
        return self.wm.range_count(il, ir, dl, du)

# PDFEND -- everything below is the test harness and is not printed in the PDF

def main() -> None:
    # --- XOR-shifted k-th smallest -------------------------------------
    a = [3, 1, 2, 1, 1, 0, 3, 2, 2]
    n = len(a)
    height = 2
    wm = WaveletMatrixSum(height, a, a[:])

    xor_ok = True
    for l in range(n):
        for r in range(l + 1, n + 1):
            for mask in range(4):
                xored = sorted(v ^ mask for v in a[l:r])
                for k in range(len(xored)):
                    got = wm.kth_smallest_xor(l, r, k, mask)
                    if got != xored[k]:
                        xor_ok = False
                        print(f"XOR MISMATCH l={l} r={r} k={k} mask={mask}: got {got} want {xored[k]}")
    print(f"XOR k-th-smallest cross-check: {'PASS' if xor_ok else 'FAIL'}")
    print(f"sanity: kth_smallest(2,7,0) = {wm.kth_smallest(2, 7, 0)}, "
          f"kth_smallest_xor(2,7,0,0) = {wm.kth_smallest_xor(2, 7, 0, 0)}")

    # --- Range sum over an attached array -------------------------------
    x_vals = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    wm_sum = WaveletMatrixSum(height, a, x_vals)
    sum_ok = True
    for l in range(n):
        for r in range(l + 1, n + 1):
            for lo in range(5):
                for hi in range(lo, 5):
                    expected = sum(x_vals[i] for i in range(l, r) if lo <= a[i] < hi)
                    got = wm_sum.range_sum(l, r, lo, hi)
                    if got != expected:
                        sum_ok = False
                        print(f"SUM MISMATCH l={l} r={r} [{lo},{hi}): got {got} want {expected}")
    print(f"range-sum cross-check: {'PASS' if sum_ok else 'FAIL'}")

    # --- 2D rectangle counting -------------------------------------------
    points = [(1, 5), (2, 3), (2, 8), (4, 1), (5, 5), (5, 9), (7, 2), (8, 8)]
    rc = RectangleCounter(points)
    rect_ok = True
    for l, r, d, u in [(0, 6, 0, 6), (2, 8, 2, 9), (1, 9, 0, 10), (4, 5, 0, 10)]:
        expected = sum(1 for x, y in points if l <= x < r and d <= y < u)
        got = rc.query(l, r, d, u)
        if got != expected:
            rect_ok = False
            print(f"RECT MISMATCH ({l},{r},{d},{u}): got {got} want {expected}")
        print(f"rect [{l},{r}) x [{d},{u}) -> {got} points")
    print(f"rectangle-count cross-check: {'PASS' if rect_ok else 'FAIL'}")

if __name__ == "__main__":
    main()
