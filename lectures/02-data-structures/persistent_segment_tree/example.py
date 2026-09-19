"""Persistent segment tree: path copying turns an O(log N)-time update into
an O(log N)-time, O(log N)-space *new version*, while every earlier version
stays intact and queryable. Two flavors are demonstrated:

 1. PersistentArray: a branching multi-version array with point update and
    range-min query (the direct generalization of an ordinary pointer-based
    segment tree).
 2. PersistentFrequencyTree: a chain of versions S_0, S_1, ..., S_N where
    S_i is the frequency array of A_1..A_i, giving range "count <= x" and
    range k-th-smallest queries by comparing two versions.

Both use an arena of nodes addressed by index, with index 0 reserved as the
"empty subtree" sentinel -- the standard way to avoid materializing every
empty leaf.
"""
from __future__ import annotations

NIL = 0
INF = float("inf")


class PersistentArray:
    """Branching multi-version array: point update, range min."""

    def __init__(self) -> None:
        # node[0] is the NIL sentinel: an "empty" subtree contributes the
        # monoid identity (here, +infinity for min).
        self.left = [NIL]
        self.right = [NIL]
        self.val = [INF]

    def _new_node(self, left: int, right: int, val: float) -> int:
        self.left.append(left)
        self.right.append(right)
        self.val.append(val)
        return len(self.val) - 1

    def build(self, lo: int, hi: int, a: list[int]) -> int:
        """Build the initial version from a 1-indexed array."""
        if lo == hi:
            return self._new_node(NIL, NIL, a[lo])
        mid = (lo + hi) // 2
        l = self.build(lo, mid, a)
        r = self.build(mid + 1, hi, a)
        return self._new_node(l, r, min(self.val[l], self.val[r]))

    def update(self, root: int, lo: int, hi: int, pos: int, x: int) -> int:
        """Return the root of a NEW version with position `pos` set to `x`,
        built from an existing version's root `root`. Only the O(log N)
        nodes on the path to `pos` are copied; everything else is shared.
        """
        if lo == hi:
            return self._new_node(NIL, NIL, x)
        mid = (lo + hi) // 2
        old_l, old_r = self.left[root], self.right[root]
        if pos <= mid:
            new_l, new_r = self.update(old_l, lo, mid, pos, x), old_r
        else:
            new_l, new_r = old_l, self.update(old_r, mid + 1, hi, pos, x)
        return self._new_node(new_l, new_r, min(self.val[new_l], self.val[new_r]))

    def query(self, root: int, lo: int, hi: int, l: int, r: int) -> float:
        """Range min of [l, r] under the version rooted at `root`."""
        if r < lo or hi < l:
            return INF
        if l <= lo and hi <= r:
            return self.val[root]
        mid = (lo + hi) // 2
        return min(
            self.query(self.left[root], lo, mid, l, r),
            self.query(self.right[root], mid + 1, hi, l, r),
        )


class PersistentFrequencyTree:
    """Chained frequency-array versions: range count <= x, range k-th
    smallest.
    """

    def __init__(self) -> None:
        self.left = [NIL]
        self.right = [NIL]
        self.count = [0]

    def _new_node(self, left: int, right: int, count: int) -> int:
        self.left.append(left)
        self.right.append(right)
        self.count.append(count)
        return len(self.count) - 1

    def insert(self, prev: int, lo: int, hi: int, value: int) -> int:
        """New version = old version (root `prev`) with one more occurrence
        of `value` recorded.
        """
        base = self.count[prev]
        if lo == hi:
            return self._new_node(NIL, NIL, base + 1)
        mid = (lo + hi) // 2
        old_l = self.left[prev] if prev != NIL else NIL
        old_r = self.right[prev] if prev != NIL else NIL
        if value <= mid:
            new_l, new_r = self.insert(old_l, lo, mid, value), old_r
        else:
            new_l, new_r = old_l, self.insert(old_r, mid + 1, hi, value)
        return self._new_node(new_l, new_r, base + 1)

    def count_less(self, root: int, lo: int, hi: int, threshold: int) -> int:
        """Count of values < threshold recorded in the version rooted at
        `root`."""
        if root == NIL or threshold <= lo:
            return 0
        if hi < threshold:
            return self.count[root]
        mid = (lo + hi) // 2
        return self.count_less(self.left[root], lo, mid, threshold) + self.count_less(
            self.right[root], mid + 1, hi, threshold
        )

    def kth_smallest(self, hi_root: int, lo_root: int, lo: int, hi: int, k: int) -> int:
        """The k-th smallest (1-indexed) value among positions counted by
        `hi_root` but not by `lo_root`. Descends BOTH trees at once,
        comparing the difference in left-subtree counts against k.
        """
        if lo == hi:
            return lo
        mid = (lo + hi) // 2
        left_hi = self.left[hi_root] if hi_root != NIL else NIL
        left_lo = self.left[lo_root] if lo_root != NIL else NIL
        left_count = self.count[left_hi] - self.count[left_lo]
        if k <= left_count:
            return self.kth_smallest(left_hi, left_lo, lo, mid, k)
        right_hi = self.right[hi_root] if hi_root != NIL else NIL
        right_lo = self.right[lo_root] if lo_root != NIL else NIL
        return self.kth_smallest(right_hi, right_lo, mid + 1, hi, k - left_count)


def main() -> None:
    # --- Part 1: branching multi-version array, range min --------------
    a = [0, 5, 9, 4, 7]  # 1-indexed; index 0 unused
    pa = PersistentArray()
    v0 = pa.build(1, 4, a)

    # Version 1: from v0, set position 2 to 2.
    v1 = pa.update(v0, 1, 4, 2, 2)
    # Version 2: branches from v0 too (not from v1!), set position 4 to 1.
    v2 = pa.update(v0, 1, 4, 4, 1)

    print(f"v0 = (5,9,4,7): range min [1,4] = {pa.query(v0, 1, 4, 1, 4)}")
    print(f"v1 = (5,2,4,7): range min [1,2] = {pa.query(v1, 1, 4, 1, 2)}")
    print(f"v2 = (5,9,4,1): range min [3,4] = {pa.query(v2, 1, 4, 3, 4)}")
    print(f"v0 is untouched: range min [1,4] = {pa.query(v0, 1, 4, 1, 4)} (still 4)")

    # --- Part 2: frequency-array chain, count<=x and k-th smallest ------
    vals = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
    n = len(vals)
    domain_hi = max(vals)
    ft = PersistentFrequencyTree()
    roots = [NIL] * (n + 1)  # roots[i] = version after inserting vals[0:i]
    for i in range(n):
        roots[i + 1] = ft.insert(roots[i], 1, domain_hi, vals[i])

    all_ok = True
    for l in range(1, n + 1):
        for r in range(l, n + 1):
            slice_ = vals[l - 1 : r]
            for x in range(1, domain_hi + 1):
                expected = sum(1 for v in slice_ if v <= x)
                got = ft.count_less(roots[r], 1, domain_hi, x + 1) - ft.count_less(
                    roots[l - 1], 1, domain_hi, x + 1
                )
                if got != expected:
                    all_ok = False
                    print(f"COUNT MISMATCH l={l} r={r} x={x}: got {got} want {expected}")
            for idx, want in enumerate(sorted(slice_)):
                k = idx + 1
                got = ft.kth_smallest(roots[r], roots[l - 1], 1, domain_hi, k)
                if got != want:
                    all_ok = False
                    print(f"KTH MISMATCH l={l} r={r} k={k}: got {got} want {want}")

    print(f"\nfrequency-tree cross-check over all (l,r,x) and (l,r,k): {'PASS' if all_ok else 'FAIL'}")

    print(f"\nvals = {vals}")
    print(
        "count of values <= 4 in [3,8] (1-indexed) = "
        f"{ft.count_less(roots[8], 1, domain_hi, 5) - ft.count_less(roots[2], 1, domain_hi, 5)}"
    )
    print(f"3rd smallest value in [3,8] (1-indexed) = {ft.kth_smallest(roots[8], roots[2], 1, domain_hi, 3)}")


if __name__ == "__main__":
    main()
