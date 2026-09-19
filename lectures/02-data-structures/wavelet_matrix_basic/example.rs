// Wavelet Matrix: answers range k-th-smallest and range count-in-[a,b)
// queries on a static array of small non-negative integers, by re-encoding
// the array as one bit-plane per level plus a rank (prefix-count) table.

struct WaveletMatrix {
    height: u32,
    n: usize,
    // data[h][i] = number of 1-bits among the first i elements of the
    // level-h bit-plane (a "rank" prefix-sum table).
    data: Vec<Vec<u32>>,
}

impl WaveletMatrix {
    fn new(height: u32, mut a: Vec<u64>) -> Self {
        let n = a.len();
        let mut data = vec![Vec::new(); height as usize];

        // Process levels from the most significant bit down. At each level,
        // stably partition the current order into "bit = 0" then "bit = 1",
        // and record the bit-plane's prefix sum before overwriting `a`.
        for h in (0..height).rev() {
            let mut dir = vec![0u8; n];
            let (mut left, mut right) = (Vec::new(), Vec::new());
            for i in 0..n {
                let bit = ((a[i] >> h) & 1) as u8;
                dir[i] = bit;
                if bit == 0 {
                    left.push(a[i]);
                } else {
                    right.push(a[i]);
                }
            }
            let split = left.len();
            for (i, v) in left.into_iter().enumerate() {
                a[i] = v;
            }
            for (i, v) in right.into_iter().enumerate() {
                a[split + i] = v;
            }

            let mut prefix = vec![0u32; n + 1];
            for i in 0..n {
                prefix[i + 1] = prefix[i] + dir[i] as u32;
            }
            data[h as usize] = prefix;
        }

        WaveletMatrix { height, n, data }
    }

    /// Given a node at height h+1 covering [l, r), return the ranges its
    /// left and right children cover at height h: (left_l, left_r, right_l, right_r).
    fn subtree_range(&self, h: u32, l: usize, r: usize) -> (usize, usize, usize, usize) {
        let dat = &self.data[h as usize];
        let a1 = dat[l] as usize;
        let a0 = l - a1;
        let b1 = dat[r] as usize;
        let b0 = r - b1;
        let c0 = self.n - *dat.last().unwrap() as usize; // total zeros at this level
        (a0, b0, c0 + a1, c0 + b1)
    }

    /// 0-indexed k-th smallest value among A[l..r).
    fn kth_smallest(&self, mut l: usize, mut r: usize, mut k: usize) -> u64 {
        let mut value = 0u64;
        for h in (0..self.height).rev() {
            let (l0, r0, l1, r1) = self.subtree_range(h, l, r);
            let left_size = r0 - l0;
            if k < left_size {
                l = l0;
                r = r0;
            } else {
                k -= left_size;
                value |= 1u64 << h;
                l = l1;
                r = r1;
            }
        }
        value
    }

    /// Count of A[l..r) strictly less than x.
    fn count_less(&self, mut l: usize, mut r: usize, x: u64) -> usize {
        let limit = 1u64 << self.height;
        if x == 0 {
            return 0;
        }
        if x >= limit {
            return r - l;
        }
        let mut count = 0usize;
        for h in (0..self.height).rev() {
            let bit = (x >> h) & 1;
            let (l0, r0, l1, r1) = self.subtree_range(h, l, r);
            if bit == 1 {
                count += r0 - l0; // everything landing in the 0-branch is < x
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        count
    }

    /// Count of A[l..r) with a <= value < b.
    fn range_count(&self, l: usize, r: usize, a: u64, b: u64) -> usize {
        self.count_less(l, r, b) - self.count_less(l, r, a)
    }
}

fn main() {
    let a: Vec<u64> = vec![3, 1, 2, 1, 1, 0, 3, 2, 2];
    let n = a.len();
    let height = 2u32; // values fit in 0..4

    let wm = WaveletMatrix::new(height, a.clone());

    // Cross-check against brute force on every (l, r, k) and (l, r, a, b).
    let mut all_ok = true;
    for l in 0..n {
        for r in (l + 1)..=n {
            let mut sorted_slice: Vec<u64> = a[l..r].to_vec();
            sorted_slice.sort();
            for k in 0..sorted_slice.len() {
                let got = wm.kth_smallest(l, r, k);
                if got != sorted_slice[k] {
                    all_ok = false;
                    println!("MISMATCH kth_smallest({l},{r},{k}): got {got}, want {}", sorted_slice[k]);
                }
            }
            for lo in 0u64..=4 {
                for hi in lo..=4 {
                    let expected = a[l..r].iter().filter(|&&v| lo <= v && v < hi).count();
                    let got = wm.range_count(l, r, lo, hi);
                    if got != expected {
                        all_ok = false;
                        println!("MISMATCH range_count({l},{r},{lo},{hi}): got {got}, want {expected}");
                    }
                }
            }
        }
    }
    println!("brute-force cross-check over all queries: {}", if all_ok { "PASS" } else { "FAIL" });

    println!("\nA = {:?}", a);
    println!("kth_smallest(2, 7, 0) = {}  (median-ish of A[2..7))", wm.kth_smallest(2, 7, 0));
    println!("kth_smallest(2, 7, 4) = {}  (largest of A[2..7))", wm.kth_smallest(2, 7, 4));
    println!("range_count(0, 9, 1, 3) = {}  (values in [1,3) over all of A)", wm.range_count(0, 9, 1, 3));
}
