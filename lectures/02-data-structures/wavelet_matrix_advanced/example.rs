// Wavelet Matrix, advanced queries: XOR-shifted k-th smallest, range sums
// over an attached value array, and 2D rectangle counting via coordinate
// compression. Builds on the plain range k-th-smallest / range-count
// wavelet matrix from the companion "basics" lecture.

struct WaveletMatrixSum {
    height: u32,
    n: usize,
    // dat[h][i]  = number of 1-bits among the first i elements of the
    //              level-h bit-plane (same rank table as the basic matrix).
    dat: Vec<Vec<u32>>,
    // sum0[h][i] = sum of X over the first i elements of level h's
    //              *pre-partition* order, restricted to those with bit 0.
    //              This is exactly what a "count < threshold" walk needs,
    //              weighted by X instead of by 1.
    sum0: Vec<Vec<i64>>,
    // total_prefix[i] = sum of X_0 .. X_{i-1} in the ORIGINAL index order
    //                   (unpermuted). Used to answer threshold >= 2^height,
    //                   where every element qualifies.
    total_prefix: Vec<i64>,
}

impl WaveletMatrixSum {
    fn new(height: u32, mut a: Vec<u64>, mut x: Vec<i64>) -> Self {
        let n = a.len();
        let mut total_prefix = vec![0i64; n + 1];
        for i in 0..n {
            total_prefix[i + 1] = total_prefix[i] + x[i];
        }

        let mut dat = vec![Vec::new(); height as usize];
        let mut sum0 = vec![Vec::new(); height as usize];

        for h in (0..height).rev() {
            let mut dir = vec![0u8; n];
            let mut dat_h = vec![0u32; n + 1];
            let mut sum0_h = vec![0i64; n + 1];
            for i in 0..n {
                let bit = ((a[i] >> h) & 1) as u8;
                dir[i] = bit;
                dat_h[i + 1] = dat_h[i] + bit as u32;
                sum0_h[i + 1] = sum0_h[i] + if bit == 0 { x[i] } else { 0 };
            }
            dat[h as usize] = dat_h;
            sum0[h as usize] = sum0_h;

            let (mut left_a, mut right_a) = (Vec::new(), Vec::new());
            let (mut left_x, mut right_x) = (Vec::new(), Vec::new());
            for i in 0..n {
                if dir[i] == 0 {
                    left_a.push(a[i]);
                    left_x.push(x[i]);
                } else {
                    right_a.push(a[i]);
                    right_x.push(x[i]);
                }
            }
            let split = left_a.len();
            for (i, v) in left_a.into_iter().enumerate() {
                a[i] = v;
            }
            for (i, v) in right_a.into_iter().enumerate() {
                a[split + i] = v;
            }
            for (i, v) in left_x.into_iter().enumerate() {
                x[i] = v;
            }
            for (i, v) in right_x.into_iter().enumerate() {
                x[split + i] = v;
            }
        }

        WaveletMatrixSum { height, n, dat, sum0, total_prefix }
    }

    fn subtree_range(&self, h: u32, l: usize, r: usize) -> (usize, usize, usize, usize) {
        let dat = &self.dat[h as usize];
        let a1 = dat[l] as usize;
        let a0 = l - a1;
        let b1 = dat[r] as usize;
        let b0 = r - b1;
        let c0 = self.n - *dat.last().unwrap() as usize;
        (a0, b0, c0 + a1, c0 + b1)
    }

    /// Plain range k-th smallest (as in the basics lecture).
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

    /// k-th smallest of (A_i XOR mask) over i in [l, r). Bit h of A_i XOR
    /// mask equals bit h of A_i, flipped whenever mask's bit h is 1 -- so
    /// we just swap which child counts as "left" at levels where mask has
    /// a 1 bit.
    fn kth_smallest_xor(&self, mut l: usize, mut r: usize, mut k: usize, mask: u64) -> u64 {
        let mut value = 0u64;
        for h in (0..self.height).rev() {
            let (l0, r0, l1, r1) = self.subtree_range(h, l, r);
            let flip = (mask >> h) & 1 == 1;
            let (real_l0, real_r0, real_l1, real_r1) =
                if flip { (l1, r1, l0, r0) } else { (l0, r0, l1, r1) };
            let left_size = real_r0 - real_l0;
            if k < left_size {
                l = real_l0;
                r = real_r0;
            } else {
                k -= left_size;
                value |= 1u64 << h;
                l = real_l1;
                r = real_r1;
            }
        }
        value
    }

    /// Sum of X_i over i in [l, r) with A_i < threshold.
    fn sum_less(&self, mut l: usize, mut r: usize, threshold: u64) -> i64 {
        let limit = 1u64 << self.height;
        if threshold == 0 {
            return 0;
        }
        if threshold >= limit {
            // Every element satisfies A_i < 2^height, so this is a plain
            // range sum over the original (unpermuted) index order.
            return self.total_prefix[r] - self.total_prefix[l];
        }
        let mut total = 0i64;
        for h in (0..self.height).rev() {
            let bit = (threshold >> h) & 1;
            let (l0, r0, l1, r1) = self.subtree_range(h, l, r);
            if bit == 1 {
                total += self.sum0[h as usize][r] - self.sum0[h as usize][l];
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        total
    }

    /// Sum of X_i over i in [l, r) with a <= A_i < b.
    fn range_sum(&self, l: usize, r: usize, a: u64, b: u64) -> i64 {
        self.sum_less(l, r, b) - self.sum_less(l, r, a)
    }

    /// Count of i in [l, r) with a <= A_i < b (range_sum with X == all ones).
    fn range_count(&self, l: usize, r: usize, a: u64, b: u64) -> i64 {
        self.range_sum(l, r, a, b)
    }
}

/// 2D rectangle counting: given N points, answer "how many points have
/// L <= x < R and D <= y < U" by coordinate-compressing y, sorting by x,
/// and reusing the 1D wavelet matrix range-count query.
struct RectangleCounter {
    sorted_x: Vec<i64>,
    sorted_unique_y: Vec<i64>,
    wm: WaveletMatrixSum,
}

impl RectangleCounter {
    fn new(points: &[(i64, i64)]) -> Self {
        let mut pts = points.to_vec();
        pts.sort_by_key(|&(x, _)| x);
        let sorted_x: Vec<i64> = pts.iter().map(|&(x, _)| x).collect();

        let mut sorted_unique_y: Vec<i64> = pts.iter().map(|&(_, y)| y).collect();
        sorted_unique_y.sort();
        sorted_unique_y.dedup();

        let ranks: Vec<u64> = pts
            .iter()
            .map(|&(_, y)| sorted_unique_y.partition_point(|&v| v < y) as u64)
            .collect();

        let n = pts.len();
        let mut height = 1u32;
        while (1usize << height) < sorted_unique_y.len().max(1) {
            height += 1;
        }
        let ones = vec![1i64; n];
        let wm = WaveletMatrixSum::new(height, ranks, ones);

        RectangleCounter { sorted_x, sorted_unique_y, wm }
    }

    fn query(&self, l: i64, r: i64, d: i64, u: i64) -> i64 {
        let il = self.sorted_x.partition_point(|&v| v < l);
        let ir = self.sorted_x.partition_point(|&v| v < r);
        let dl = self.sorted_unique_y.partition_point(|&v| v < d) as u64;
        let du = self.sorted_unique_y.partition_point(|&v| v < u) as u64;
        self.wm.range_count(il, ir, dl, du)
    }
}

fn main() {
    // --- XOR-shifted k-th smallest -------------------------------------
    let a: Vec<u64> = vec![3, 1, 2, 1, 1, 0, 3, 2, 2];
    let n = a.len();
    let height = 2u32;
    let x_attached: Vec<i64> = a.iter().map(|&v| v as i64).collect();
    let wm = WaveletMatrixSum::new(height, a.clone(), x_attached);

    let mut xor_ok = true;
    for l in 0..n {
        for r in (l + 1)..=n {
            for mask in 0u64..4 {
                let mut xored: Vec<u64> = a[l..r].iter().map(|&v| v ^ mask).collect();
                xored.sort();
                for k in 0..xored.len() {
                    let got = wm.kth_smallest_xor(l, r, k, mask);
                    if got != xored[k] {
                        xor_ok = false;
                        println!("XOR MISMATCH l={l} r={r} k={k} mask={mask}: got {got} want {}", xored[k]);
                    }
                }
            }
        }
    }
    println!("XOR k-th-smallest cross-check: {}", if xor_ok { "PASS" } else { "FAIL" });
    println!(
        "sanity: kth_smallest(2,7,0) = {}, kth_smallest_xor(2,7,0,0) = {}",
        wm.kth_smallest(2, 7, 0),
        wm.kth_smallest_xor(2, 7, 0, 0)
    );

    // --- Range sum over an attached array -------------------------------
    let x_vals: Vec<i64> = vec![10, 20, 30, 40, 50, 60, 70, 80, 90];
    let wm_sum = WaveletMatrixSum::new(height, a.clone(), x_vals.clone());
    let mut sum_ok = true;
    for l in 0..n {
        for r in (l + 1)..=n {
            for lo in 0u64..=4 {
                for hi in lo..=4 {
                    let expected: i64 = (l..r).filter(|&i| lo <= a[i] && a[i] < hi).map(|i| x_vals[i]).sum();
                    let got = wm_sum.range_sum(l, r, lo, hi);
                    if got != expected {
                        sum_ok = false;
                        println!("SUM MISMATCH l={l} r={r} [{lo},{hi}): got {got} want {expected}");
                    }
                }
            }
        }
    }
    println!("range-sum cross-check: {}", if sum_ok { "PASS" } else { "FAIL" });

    // --- 2D rectangle counting -------------------------------------------
    let points: Vec<(i64, i64)> = vec![(1, 5), (2, 3), (2, 8), (4, 1), (5, 5), (5, 9), (7, 2), (8, 8)];
    let rc = RectangleCounter::new(&points);
    let mut rect_ok = true;
    let queries = [(0i64, 6, 0i64, 6), (2, 8, 2, 9), (1, 9, 0, 10), (4, 5, 0, 10)];
    for &(l, r, d, u) in &queries {
        let expected = points.iter().filter(|&&(x, y)| l <= x && x < r && d <= y && y < u).count() as i64;
        let got = rc.query(l, r, d, u);
        if got != expected {
            rect_ok = false;
            println!("RECT MISMATCH ({l},{r},{d},{u}): got {got} want {expected}");
        }
        println!("rect [{l},{r}) x [{d},{u}) -> {got} points");
    }
    println!("rectangle-count cross-check: {}", if rect_ok { "PASS" } else { "FAIL" });
}
