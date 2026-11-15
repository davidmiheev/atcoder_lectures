// The Cycle Lemma: for an integer sequence A of length N with every
// A_i <= 1 and nonnegative sum K, exactly K of its N cyclic shifts have
// every prefix sum (after the empty prefix) strictly positive. Several
// classical counting formulas -- Catalan, Fuss-Catalan, Narayana, and the
// Catalan convolution identity -- follow from it via a "prepend +1" trick.

// PDFBEGIN
fn binomial(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    let mut row = vec![0u128; n + 1];
    row[0] = 1;
    for i in 1..=n {
        for j in (1..=i).rev() {
            row[j] += row[j - 1];
        }
    }
    row[k]
}

fn catalan(n: usize) -> u128 {
    binomial(2 * n, n) / (n as u128 + 1)
}

fn prefix_always_nonneg(a: &[i64]) -> bool {
    let mut s = 0i64;
    for &x in a {
        s += x;
        if s < 0 {
            return false;
        }
    }
    true
}

fn prefix_positive_after_first(a: &[i64]) -> bool {
    let mut s = 0i64;
    for &x in a {
        s += x;
        if s <= 0 {
            return false;
        }
    }
    true
}

/// The Cycle Lemma count: how many of A's N cyclic shifts have every
/// prefix sum (excluding S_0) strictly positive.
fn cycle_lemma_count(a: &[i64]) -> usize {
    let n = a.len();
    (0..n)
        .filter(|&j| {
            let shifted: Vec<i64> = a[j..].iter().chain(a[..j].iter()).cloned().collect();
            prefix_positive_after_first(&shifted)
        })
        .count()
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

/// All ways to choose `k` of `n` positions (as sorted index vectors).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![vec![]];
    }
    if k > n {
        return vec![];
    }
    let mut result = Vec::new();
    for first in 0..=(n - k) {
        for mut rest in combinations(n - first - 1, k - 1) {
            for v in rest.iter_mut() {
                *v += first + 1;
            }
            let mut combo = vec![first];
            combo.append(&mut rest);
            result.push(combo);
        }
    }
    result
}

/// Every length-`total` sequence with `plus_val` at `num_plus` positions
/// (chosen every possible way) and `minus_val` everywhere else.
fn all_sequences(total: usize, num_plus: usize, plus_val: i64, minus_val: i64) -> Vec<Vec<i64>> {
    combinations(total, num_plus)
        .into_iter()
        .map(|combo| {
            let mut seq = vec![minus_val; total];
            for idx in combo {
                seq[idx] = plus_val;
            }
            seq
        })
        .collect()
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

fn main() {
    // --- Core theorem, on the lecture's own worked example ---------------
    let example = vec![1i64, 0, -2, 1, -1, 1, 0, 0, 1, 1];
    let k: i64 = example.iter().sum();
    let got = cycle_lemma_count(&example);
    println!("lecture's own example: sum K = {k}, cyclic shifts with all-positive prefix sums = {got}");
    assert_eq!(k, 2);
    assert_eq!(got, 2);

    // --- Core theorem, on random sequences with A_i <= 1, sum >= 0 -------
    let mut rng = Rng(99);
    let mut core_ok = true;
    for _ in 0..500 {
        let n = 1 + (rng.next() % 12) as usize;
        let mut a: Vec<i64> = (0..n).map(|_| rng.range(-3, 1)).collect();
        let sum: i64 = a.iter().sum();
        if sum < 0 {
            // bump one entry up (still <= 1) until the sum is nonnegative
            let deficit = -sum;
            let idx = (rng.next() as usize) % n;
            a[idx] += deficit;
            if a[idx] > 1 {
                continue; // skip the rare case this pushes a value above 1
            }
        }
        let k = a.iter().sum::<i64>();
        let got = cycle_lemma_count(&a) as i64;
        if got != k {
            core_ok = false;
            println!("CYCLE LEMMA MISMATCH a={:?}: K={k} got={got}", a);
        }
    }
    println!("Cycle Lemma holds on random sequences (A_i <= 1, sum >= 0): {}", if core_ok { "PASS" } else { "FAIL" });

    // --- Catalan numbers: both formulas from Theorem 6 --------------------
    let mut catalan_ok = true;
    for n in 0..=6usize {
        // Condition 1: correct bracket sequences of length 2n.
        let count1 = if n == 0 {
            1
        } else {
            all_sequences(2 * n, n, 1, -1).iter().filter(|a| prefix_always_nonneg(a)).count()
        };
        // Condition 2: length 2n+1, n+1 plus-ones, prefix positive after first.
        let count2 = all_sequences(2 * n + 1, n + 1, 1, -1)
            .iter()
            .filter(|a| prefix_positive_after_first(a))
            .count();
        // Condition 4: length 2n, n ones and n minus-ones... no: condition 4
        // is length 2n+1 sequences (n+1 plus, n minus) that START with +1.
        let count4 = all_sequences(2 * n + 1, n + 1, 1, -1).iter().filter(|a| a[0] == 1).count();

        let formula = catalan(n) as usize;
        let total3 = binomial(2 * n + 1, n) as usize; // condition 3: all such sequences

        let cycle_lemma_relation = count2 * (2 * n + 1) == total3;
        let second_relation = count2 * (n + 1) == count4;
        if count1 != formula || count2 != formula || count4 != binomial(2 * n, n) as usize || !cycle_lemma_relation || !second_relation {
            catalan_ok = false;
            println!(
                "CATALAN MISMATCH n={n}: count1={count1} count2={count2} count4={count4} formula={formula} cycle_lemma_relation={cycle_lemma_relation} second_relation={second_relation}"
            );
        }
    }
    println!("Catalan numbers: both formulas C_N = C(2N+1,N)/(2N+1) = C(2N,N)/(N+1) verified via Cycle Lemma: {}", if catalan_ok { "PASS" } else { "FAIL" });

    // --- Fuss-Catalan numbers (Theorem 8) ----------------------------------
    let mut fuss_ok = true;
    for n in 1..=4usize {
        for m in 1..=3usize {
            let total = n * m + n;
            let num_plus = n * m;
            let count = all_sequences(total, num_plus, 1, -(m as i64))
                .iter()
                .filter(|a| prefix_always_nonneg(a))
                .count();
            let formula = (binomial(total, n) / (n as u128 * m as u128 + 1)) as usize;
            if count != formula {
                fuss_ok = false;
                println!("FUSS-CATALAN MISMATCH n={n} m={m}: count={count} formula={formula}");
            }
        }
    }
    println!("Fuss-Catalan numbers (Theorem 8) verified for small N,M: {}", if fuss_ok { "PASS" } else { "FAIL" });

    // --- Narayana numbers (Theorem 9) --------------------------------------
    let mut narayana_ok = true;
    for n in 1..=6usize {
        for m in 1..=n {
            let count = all_sequences(2 * n, n, 1, -1)
                .iter()
                .filter(|a| {
                    prefix_always_nonneg(a) && a.windows(2).filter(|w| w[0] == 1 && w[1] == -1).count() == m
                })
                .count();
            let formula = (binomial(n, m) * binomial(n, m - 1) / n as u128) as usize;
            if count != formula {
                narayana_ok = false;
                println!("NARAYANA MISMATCH n={n} m={m}: count={count} formula={formula}");
            }
        }
    }
    println!("Narayana numbers (Theorem 9) verified for n<=6: {}", if narayana_ok { "PASS" } else { "FAIL" });

    // --- Catalan convolution formula (Theorem 7) ---------------------------
    fn compositions(total: usize, parts: usize) -> Vec<Vec<usize>> {
        if parts == 1 {
            return vec![vec![total]];
        }
        let mut result = Vec::new();
        for first in 0..=total {
            for mut rest in compositions(total - first, parts - 1) {
                let mut c = vec![first];
                c.append(&mut rest);
                result.push(c);
            }
        }
        result
    }
    let mut conv_ok = true;
    for n in 0..=5usize {
        for m in 1..=4usize {
            let sum: u128 = compositions(n, m).iter().map(|parts| parts.iter().map(|&x| catalan(x)).product::<u128>()).sum();
            let formula = (m as u128 * binomial(2 * n + m, n)) / (2 * n as u128 + m as u128);
            if sum != formula {
                conv_ok = false;
                println!("CONVOLUTION MISMATCH n={n} m={m}: sum={sum} formula={formula}");
            }
        }
    }
    println!("Catalan convolution formula (Theorem 7) verified for small N,M: {}", if conv_ok { "PASS" } else { "FAIL" });
}
