// Persistent segment tree: path copying turns an O(log N)-time update into
// an O(log N)-time, O(log N)-space *new version*, while every earlier
// version stays intact and queryable. Two flavors are demonstrated:
//
//  1. PersistentArray: a branching multi-version array with point update
//     and range-min query (the direct generalization of an ordinary
//     pointer-based segment tree).
//  2. PersistentFrequencyTree: a chain of versions S_0, S_1, ..., S_N where
//     S_i is the frequency array of A_1..A_i, giving range "count <= x"
//     and range k-th-smallest queries by comparing two versions.
//
// Both use an arena of nodes addressed by index, with index 0 reserved as
// the "empty subtree" sentinel -- the standard way to avoid materializing
// every empty leaf.

// PDFBEGIN
const NIL: usize = 0;

// ---------------------------------------------------------------------
// 1. Branching multi-version array (point update, range min)
// ---------------------------------------------------------------------

#[derive(Clone, Copy)]

struct MinNode {
    left: usize,
    right: usize,
    val: i64,
}

struct PersistentArray {
    lo: usize,
    hi: usize, // domain is 1-indexed positions [lo, hi]
    nodes: Vec<MinNode>,
}

const INF: i64 = i64::MAX;

impl PersistentArray {
    fn new(n: usize) -> Self {
        // nodes[0] is the NIL sentinel: an "empty" subtree contributes the
        // monoid identity (here, +infinity for min).
        PersistentArray { lo: 1, hi: n, nodes: vec![MinNode { left: NIL, right: NIL, val: INF }] }
    }

    fn val(&self, id: usize) -> i64 {
        self.nodes[id].val
    }

    /// Build the initial version (version 0) from a 1-indexed array.
    fn build(&mut self, lo: usize, hi: usize, a: &[i64]) -> usize {
        if lo == hi {
            let id = self.nodes.len();
            self.nodes.push(MinNode { left: NIL, right: NIL, val: a[lo] });
            return id;
        }
        let mid = (lo + hi) / 2;
        let l = self.build(lo, mid, a);
        let r = self.build(mid + 1, hi, a);
        let id = self.nodes.len();
        self.nodes.push(MinNode { left: l, right: r, val: self.val(l).min(self.val(r)) });
        id
    }

    /// Return the root of a NEW version with position `pos` set to `x`,
    /// built from an existing version's root `root`. Only the O(log N)
    /// nodes on the path to `pos` are copied; everything else is shared.
    fn update(&mut self, root: usize, lo: usize, hi: usize, pos: usize, x: i64) -> usize {
        if lo == hi {
            let id = self.nodes.len();
            self.nodes.push(MinNode { left: NIL, right: NIL, val: x });
            return id;
        }
        let mid = (lo + hi) / 2;
        let (old_l, old_r) = (self.nodes[root].left, self.nodes[root].right);
        let (new_l, new_r) = if pos <= mid {
            (self.update(old_l, lo, mid, pos, x), old_r)
        } else {
            (old_l, self.update(old_r, mid + 1, hi, pos, x))
        };
        let id = self.nodes.len();
        self.nodes.push(MinNode { left: new_l, right: new_r, val: self.val(new_l).min(self.val(new_r)) });
        id
    }

    /// Range min of [l, r] under the version rooted at `root`.
    fn query(&self, root: usize, lo: usize, hi: usize, l: usize, r: usize) -> i64 {
        if r < lo || hi < l {
            return INF;
        }
        if l <= lo && hi <= r {
            return self.val(root);
        }
        let mid = (lo + hi) / 2;
        let node = self.nodes[root];
        self.query(node.left, lo, mid, l, r).min(self.query(node.right, mid + 1, hi, l, r))
    }
}

// ---------------------------------------------------------------------
// 2. Chained frequency-array versions (range count <= x, range k-th smallest)
// ---------------------------------------------------------------------

#[derive(Clone, Copy)]

struct CountNode {
    left: usize,
    right: usize,
    count: i64,
}

struct PersistentFrequencyTree {
    lo: usize,
    hi: usize, // value domain [lo, hi]
    nodes: Vec<CountNode>,
}

impl PersistentFrequencyTree {
    fn new(lo: usize, hi: usize) -> Self {
        PersistentFrequencyTree { lo, hi, nodes: vec![CountNode { left: NIL, right: NIL, count: 0 }] }
    }

    fn count(&self, id: usize) -> i64 {
        self.nodes[id].count
    }

    /// New version = old version (root `prev`) with one more occurrence of
    /// `value` recorded.
    fn insert(&mut self, prev: usize, lo: usize, hi: usize, value: usize) -> usize {
        let base = if prev == NIL { 0 } else { self.nodes[prev].count };
        if lo == hi {
            let id = self.nodes.len();
            self.nodes.push(CountNode { left: NIL, right: NIL, count: base + 1 });
            return id;
        }
        let mid = (lo + hi) / 2;
        let (old_l, old_r) = if prev == NIL { (NIL, NIL) } else { (self.nodes[prev].left, self.nodes[prev].right) };
        let (new_l, new_r) = if value <= mid {
            (self.insert(old_l, lo, mid, value), old_r)
        } else {
            (old_l, self.insert(old_r, mid + 1, hi, value))
        };
        let id = self.nodes.len();
        self.nodes.push(CountNode { left: new_l, right: new_r, count: base + 1 });
        id
    }

    /// Count of values < threshold recorded in the version rooted at `root`.
    fn count_less(&self, root: usize, lo: usize, hi: usize, threshold: usize) -> i64 {
        if root == NIL || threshold <= lo {
            return 0;
        }
        if hi < threshold {
            return self.count(root);
        }
        let mid = (lo + hi) / 2;
        let node = self.nodes[root];
        self.count_less(node.left, lo, mid, threshold) + self.count_less(node.right, mid + 1, hi, threshold)
    }

    /// The k-th smallest (1-indexed) value among positions counted by
    /// `hi_root` but not by `lo_root` -- i.e. among a half-open range of
    /// insert-versions [lo_version, hi_version). Descends BOTH trees at
    /// once, comparing the difference in left-subtree counts against k.
    fn kth_smallest(&self, hi_root: usize, lo_root: usize, lo: usize, hi: usize, k: i64) -> usize {
        if lo == hi {
            return lo;
        }
        let mid = (lo + hi) / 2;
        let left_hi = if hi_root == NIL { NIL } else { self.nodes[hi_root].left };
        let left_lo = if lo_root == NIL { NIL } else { self.nodes[lo_root].left };
        let left_count = self.count(left_hi) - self.count(left_lo);
        if k <= left_count {
            self.kth_smallest(left_hi, left_lo, lo, mid, k)
        } else {
            let right_hi = if hi_root == NIL { NIL } else { self.nodes[hi_root].right };
            let right_lo = if lo_root == NIL { NIL } else { self.nodes[lo_root].right };
            self.kth_smallest(right_hi, right_lo, mid + 1, hi, k - left_count)
        }
    }
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn main() {
    // --- Part 1: branching multi-version array, range min --------------
    let a = vec![0, 5, 9, 4, 7]; // 1-indexed; index 0 unused
    let mut pa = PersistentArray::new(4);
    let v0 = pa.build(1, 4, &a);

    // Version 1: from v0, set position 2 to 2.
    let v1 = pa.update(v0, 1, 4, 2, 2);
    // Version 2: branches from v0 too (not from v1!), set position 4 to 1.
    let v2 = pa.update(v0, 1, 4, 4, 1);

    println!("v0 = (5,9,4,7): range min [1,4] = {}", pa.query(v0, 1, 4, 1, 4));
    println!("v1 = (5,2,4,7): range min [1,2] = {}", pa.query(v1, 1, 4, 1, 2));
    println!("v2 = (5,9,4,1): range min [3,4] = {}", pa.query(v2, 1, 4, 3, 4));
    println!("v0 is untouched: range min [1,4] = {} (still 4)", pa.query(v0, 1, 4, 1, 4));

    // --- Part 2: frequency-array chain, count<=x and k-th smallest ------
    let vals = vec![3usize, 1, 4, 1, 5, 9, 2, 6, 5, 3];
    let n = vals.len();
    let domain_hi = *vals.iter().max().unwrap();
    let mut ft = PersistentFrequencyTree::new(1, domain_hi);
    let mut roots = vec![NIL; n + 1]; // roots[i] = version after inserting vals[0..i]
    for i in 0..n {
        roots[i + 1] = ft.insert(roots[i], 1, domain_hi, vals[i]);
    }

    let mut all_ok = true;
    for l in 1..=n {
        for r in l..=n {
            let slice = &vals[l - 1..r];
            for x in 1..=domain_hi {
                let expected = slice.iter().filter(|&&v| v <= x).count() as i64;
                let got = ft.count_less(roots[r], 1, domain_hi, x + 1)
                    - ft.count_less(roots[l - 1], 1, domain_hi, x + 1);
                if got != expected {
                    all_ok = false;
                    println!("COUNT MISMATCH l={l} r={r} x={x}: got {got} want {expected}");
                }
            }
            let mut sorted_slice = slice.to_vec();
            sorted_slice.sort();
            for (idx, &want) in sorted_slice.iter().enumerate() {
                let k = (idx + 1) as i64;
                let got = ft.kth_smallest(roots[r], roots[l - 1], 1, domain_hi, k);
                if got != want {
                    all_ok = false;
                    println!("KTH MISMATCH l={l} r={r} k={k}: got {got} want {want}");
                }
            }
        }
    }
    println!("\nfrequency-tree cross-check over all (l,r,x) and (l,r,k): {}", if all_ok { "PASS" } else { "FAIL" });

    println!(
        "\nvals = {:?}\ncount of values <= 4 in [3,8] (1-indexed) = {}",
        vals,
        ft.count_less(roots[8], 1, domain_hi, 5) - ft.count_less(roots[2], 1, domain_hi, 5)
    );
    println!(
        "3rd smallest value in [3,8] (1-indexed) = {}",
        ft.kth_smallest(roots[8], roots[2], 1, domain_hi, 3)
    );
}
