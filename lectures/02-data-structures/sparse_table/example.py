"""Sparse table for idempotent range queries (here: range minimum).

Build: O(n log n). Query: O(1).
"""
from __future__ import annotations


class SparseTable:
    def __init__(self, a: list[int]) -> None:
        n = len(a)
        self.log2floor = [0] * (n + 1)
        for i in range(2, n + 1):
            self.log2floor[i] = self.log2floor[i // 2] + 1

        self.table = [a[:]]
        k = 1
        while (1 << k) <= n:
            half = 1 << (k - 1)
            length = 1 << k
            prev = self.table[k - 1]
            row = [min(prev[i], prev[i + half]) for i in range(n - length + 1)]
            self.table.append(row)
            k += 1

    def query(self, l: int, r: int) -> int:
        """Minimum of a[l:r], 0 <= l < r <= n."""
        k = self.log2floor[r - l]
        window = 1 << k
        return min(self.table[k][l], self.table[k][r - window])


def main() -> None:
    a = [5, 2, 4, 7, 1, 3, 6]
    st = SparseTable(a)

    for l, r in [(1, 6), (0, 7), (3, 4), (4, 6)]:
        print(f"min(a[{l}:{r}]) = {st.query(l, r)}")


if __name__ == "__main__":
    main()
