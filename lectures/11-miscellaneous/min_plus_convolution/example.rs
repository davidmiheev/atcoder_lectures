// Min-plus convolution of convex sequences: C_k = min_{i+j=k} (A_i + B_j).
// In general this takes O(NM) time (and is believed to need essentially
// that much). But if A and B are both convex ("convex below": nondecreasing
// differences), it collapses to O(N+M): reinterpret A_i as "cheapest way to
// pick i items from a bag priced by A's differences," and the convolution
// becomes "cheapest way to pick k items from the two bags combined" --
// i.e. just merge the two (already-sorted) difference sequences.

// PDFBEGIN
/// "Convex below": differences A_{i+1} - A_i are nondecreasing.
fn is_convex(a: &[i64]) -> bool {
    let diffs: Vec<i64> = a.windows(2).map(|w| w[1] - w[0]).collect();
    diffs.windows(2).all(|w| w[0] <= w[1])
}

fn differences(a: &[i64]) -> Vec<i64> {
    a.windows(2).map(|w| w[1] - w[0]).collect()
}

/// Merge two already-sorted (nondecreasing) slices in O(n+m).
fn merge_sorted(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] <= b[j] {
            result.push(a[i]);
            i += 1;
        } else {
            result.push(b[j]);
            j += 1;
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

/// Min-plus convolution of two convex sequences, O(N+M): merge their
/// (individually nondecreasing) difference sequences and re-accumulate.
fn min_plus_convolution_convex(a: &[i64], b: &[i64]) -> Vec<i64> {
    debug_assert!(is_convex(a) && is_convex(b));
    let merged_diffs = merge_sorted(&differences(a), &differences(b));
    let mut c = Vec::with_capacity(merged_diffs.len() + 1);
    c.push(a[0] + b[0]);
    for d in merged_diffs {
        c.push(c.last().unwrap() + d);
    }
    c
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

/// General min-plus convolution, O(N*M): the definition, directly.
fn min_plus_convolution_naive(a: &[i64], b: &[i64]) -> Vec<i64> {
    let (n, m) = (a.len() - 1, b.len() - 1);
    let mut c = vec![i64::MAX; n + m + 1];
    for i in 0..=n {
        for j in 0..=m {
            let s = a[i] + b[j];
            if s < c[i + j] {
                c[i + j] = s;
            }
        }
    }
    c
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
}

/// A random convex sequence of length `len`: nondecreasing random
/// differences, accumulated from a random starting value.
fn random_convex_sequence(len: usize, rng: &mut Rng) -> Vec<i64> {
    let mut diffs: Vec<i64> = (0..len.saturating_sub(1)).map(|_| rng.range(-5, 5)).collect();
    diffs.sort();
    let mut seq = vec![rng.range(-10, 10)];
    for d in diffs {
        seq.push(seq.last().unwrap() + d);
    }
    seq
}

fn main() {
    // --- The lecture's own worked examples ---------------------------------
    let a1 = vec![5i64, 1, 4, 2, 2];
    let b1 = vec![1i64, 3, 2];
    let expected1 = vec![6i64, 2, 4, 3, 3, 4, 4];
    let got1 = min_plus_convolution_naive(&a1, &b1);
    println!("worked example 1 (general): {:?} (expected {:?})", got1, expected1);
    assert_eq!(got1, expected1);

    let a2 = vec![2i64, 0, 0, 1, 2];
    let b2 = vec![0i64, 0, 1, 3];
    let expected2 = vec![2i64, 0, 0, 0, 1, 2, 3, 5];
    assert!(is_convex(&a2) && is_convex(&b2));
    let got2_fast = min_plus_convolution_convex(&a2, &b2);
    let got2_naive = min_plus_convolution_naive(&a2, &b2);
    println!("worked example 2 (convex):  {:?} (expected {:?})", got2_fast, expected2);
    assert_eq!(got2_fast, expected2);
    assert_eq!(got2_naive, expected2);

    // --- Random convex sequences: fast algorithm vs. naive definition -----
    let mut rng = Rng(2024);
    let mut convolution_ok = true;
    let mut result_convex_ok = true;
    for _ in 0..300 {
        let n = 1 + (rng.next() % 15) as usize;
        let m = 1 + (rng.next() % 15) as usize;
        let a = random_convex_sequence(n, &mut rng);
        let b = random_convex_sequence(m, &mut rng);
        let fast = min_plus_convolution_convex(&a, &b);
        let naive = min_plus_convolution_naive(&a, &b);
        if fast != naive {
            convolution_ok = false;
            println!("CONVOLUTION MISMATCH a={:?} b={:?}: fast={:?} naive={:?}", a, b, fast, naive);
        }
        if !is_convex(&fast) {
            result_convex_ok = false;
            println!("RESULT NOT CONVEX a={:?} b={:?}: c={:?}", a, b, fast);
        }
    }
    println!("fast convex-convolution matches the naive O(NM) definition: {}", if convolution_ok { "PASS" } else { "FAIL" });
    println!("the convolution of two convex sequences is itself convex: {}", if result_convex_ok { "PASS" } else { "FAIL" });

    // --- Problem 4: "choose exactly k elements, minimize the sum" ----------
    // The textbook O(N^2) DP is repeated min-plus convolution with (0, x).
    // Understanding that the dp table stays convex, and that convolving
    // with (0, x) just inserts x into the sorted difference multiset,
    // shows this dp is secretly just "sort A and take prefix sums."
    fn slow_dp(a: &[i64]) -> Vec<i64> {
        let mut dp = vec![0i64];
        for &x in a {
            let mut newdp = vec![i64::MAX; dp.len() + 1];
            for (i, &v) in dp.iter().enumerate() {
                newdp[i] = newdp[i].min(v); // don't use x
                newdp[i + 1] = newdp[i + 1].min(v + x); // use x
            }
            dp = newdp;
        }
        dp
    }
    fn fast_via_sort(a: &[i64]) -> Vec<i64> {
        let mut sorted = a.to_vec();
        sorted.sort();
        let mut dp = vec![0i64];
        for x in sorted {
            dp.push(dp.last().unwrap() + x);
        }
        dp
    }
    let mut dp_ok = true;
    let mut dp_stays_convex = true;
    for _ in 0..100 {
        let n = 1 + (rng.next() % 10) as usize;
        let a: Vec<i64> = (0..n).map(|_| rng.range(-5, 10)).collect();
        let slow = slow_dp(&a);
        let fast = fast_via_sort(&a);
        if slow != fast {
            dp_ok = false;
            println!("DP MISMATCH a={:?}: slow_dp={:?} fast_via_sort={:?}", a, slow, fast);
        }
        // Confirm the dp table is convex after every prefix of A, as the
        // theory predicts (each step is a convex-convex convolution).
        let mut dp = vec![0i64];
        for &x in &a {
            let mut newdp = vec![i64::MAX; dp.len() + 1];
            for (i, &v) in dp.iter().enumerate() {
                newdp[i] = newdp[i].min(v);
                newdp[i + 1] = newdp[i + 1].min(v + x);
            }
            dp = newdp;
            if !is_convex(&dp) {
                dp_stays_convex = false;
                println!("DP NOT CONVEX a={:?} partial dp={:?}", a, dp);
            }
        }
    }
    println!("naive O(N^2) knapsack-min dp matches sort-and-prefix-sum: {}", if dp_ok { "PASS" } else { "FAIL" });
    println!("the dp table stays convex at every step (as the theorem predicts): {}", if dp_stays_convex { "PASS" } else { "FAIL" });
}
