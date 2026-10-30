// The hook length formula: for a partition lambda of n, the number of
// standard Young tableaux of shape lambda is
//
//     f(lambda) = n! / prod over cells of (hook length of the cell).
//
// The core below computes hook lengths, the count (exactly, and modulo a
// prime), the same count from the determinant-derived product formula
// (Frobenius), and a uniform random tableau generator built on the hook
// walk of the first proof. The test harness after PDFEND cross-checks all
// of them against each other and against definition-level brute force.

// PDFBEGIN
const P: u64 = 998_244_353;

/// A partition: row lengths lambda[0] >= lambda[1] >= ... > 0 (row 0 is the top row).
type Shape = Vec<usize>;

/// Column lengths (the conjugate partition).
fn conjugate(lambda: &[usize]) -> Vec<usize> {
    let width = lambda.first().copied().unwrap_or(0);
    (0..width).map(|c| lambda.iter().filter(|&&len| len > c).count()).collect()
}

/// hook[r][c] = arm + leg + 1, where the arm is the number of cells to the
/// right of (r, c) in its row and the leg the number of cells below it in
/// its column.
fn hook_lengths(lambda: &[usize]) -> Vec<Vec<usize>> {
    let conj = conjugate(lambda);
    (0..lambda.len())
        .map(|r| (0..lambda[r]).map(|c| (lambda[r] - c - 1) + (conj[c] - r - 1) + 1).collect())
        .collect()
}

fn size(lambda: &[usize]) -> usize {
    lambda.iter().sum()
}

fn mod_pow(mut base: u64, mut exp: u64) -> u64 {
    let mut result = 1u64;
    base %= P;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % P;
        }
        base = base * base % P;
        exp >>= 1;
    }
    result
}

fn mod_inv(a: u64) -> u64 {
    mod_pow(a, P - 2)
}

/// Exact f(lambda) by the hook length formula; fine for n up to about 30.
fn count_syt_exact(lambda: &[usize]) -> u128 {
    let n = size(lambda);
    let numerator: u128 = (1..=n as u128).product();
    let denominator: u128 = hook_lengths(lambda).iter().flatten().map(|&h| h as u128).product();
    assert_eq!(numerator % denominator, 0, "hook product must divide n!");
    numerator / denominator
}

/// f(lambda) modulo the prime P, for any n < P.
fn count_syt_mod(lambda: &[usize]) -> u64 {
    let n = size(lambda) as u64;
    let numerator = (1..=n).fold(1u64, |acc, k| acc * k % P);
    let denominator = hook_lengths(lambda).iter().flatten().fold(1u64, |acc, &h| acc * h as u64 % P);
    numerator * mod_inv(denominator) % P
}

/// f(lambda) from Frobenius' product formula
///   n! * prod_{i<j}(b_i - b_j) / prod_i b_i!,   b_i = lambda_i + m - 1 - i,
/// which the determinant proof produces directly (rows indexed from 0).
fn count_syt_frobenius_mod(lambda: &[usize]) -> u64 {
    let m = lambda.len();
    let n = size(lambda) as u64;
    let b: Vec<u64> = (0..m).map(|i| (lambda[i] + m - 1 - i) as u64).collect();
    let mut numerator = (1..=n).fold(1u64, |acc, k| acc * k % P);
    let mut denominator = 1u64;
    for i in 0..m {
        for j in i + 1..m {
            numerator = numerator * (b[i] - b[j]) % P;
        }
        denominator = (1..=b[i]).fold(denominator, |acc, k| acc * k % P);
    }
    numerator * mod_inv(denominator) % P
}

/// A uniformly random standard Young tableau of shape `lambda`, built by
/// repeatedly running a hook walk from a uniformly random cell: the
/// corner where the walk stops receives the largest remaining entry.
/// `rand_below(k)` must return a uniform integer in 0..k.
fn sample_tableau(lambda: &[usize], rand_below: &mut dyn FnMut(usize) -> usize) -> Vec<Vec<usize>> {
    let mut cur: Shape = lambda.to_vec();
    let mut tableau: Vec<Vec<usize>> = lambda.iter().map(|&len| vec![0; len]).collect();
    for value in (1..=size(lambda)).rev() {
        // Pick a uniformly random cell of the current shape.
        let mut pick = rand_below(size(&cur));
        let mut r = 0;
        while pick >= cur[r] {
            pick -= cur[r];
            r += 1;
        }
        let mut c = pick;
        // Hook walk: jump to a uniformly random other cell of the hook until a corner.
        loop {
            let arm = cur[r] - c - 1;
            let leg = cur[r + 1..].iter().filter(|&&len| len > c).count();
            if arm + leg == 0 {
                break;
            }
            let k = rand_below(arm + leg);
            if k < arm {
                c += 1 + k;
            } else {
                r += 1 + (k - arm);
            }
        }
        tableau[r][c] = value;
        cur[r] -= 1;
        if cur[r] == 0 {
            cur.pop();
        }
    }
    tableau
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

use std::collections::HashMap;

/// All partitions of n in non-increasing lexicographic order.
fn partitions(n: usize) -> Vec<Shape> {
    fn go(remaining: usize, max_part: usize, current: &mut Shape, out: &mut Vec<Shape>) {
        if remaining == 0 {
            out.push(current.clone());
            return;
        }
        for part in (1..=remaining.min(max_part)).rev() {
            current.push(part);
            go(remaining - part, part, current, out);
            current.pop();
        }
    }
    let mut out = Vec::new();
    go(n, n, &mut Vec::new(), &mut out);
    out
}

/// f(lambda) = sum of f(lambda - corner) over removable corners.
fn count_syt_recurrence(lambda: &[usize], memo: &mut HashMap<Shape, u128>) -> u128 {
    if lambda.is_empty() {
        return 1;
    }
    if let Some(&v) = memo.get(lambda) {
        return v;
    }
    let mut total = 0u128;
    for r in 0..lambda.len() {
        let is_corner = r + 1 == lambda.len() || lambda[r + 1] < lambda[r];
        if is_corner {
            let mut smaller = lambda.to_vec();
            smaller[r] -= 1;
            if smaller[r] == 0 {
                smaller.pop();
            }
            total += count_syt_recurrence(&smaller, memo);
        }
    }
    memo.insert(lambda.to_vec(), total);
    total
}

/// Definition-level brute force: fill the cells in row-major order with a
/// permutation of 1..n so that every row and every column increases.
fn count_syt_brute_force(lambda: &[usize]) -> u128 {
    fn go(lambda: &[usize], grid: &mut Vec<Vec<usize>>, used: &mut Vec<bool>, r: usize, c: usize) -> u128 {
        if r == lambda.len() {
            return 1;
        }
        let (nr, nc) = if c + 1 == lambda[r] { (r + 1, 0) } else { (r, c + 1) };
        let mut total = 0;
        for v in 1..used.len() {
            if used[v] {
                continue;
            }
            if c > 0 && grid[r][c - 1] > v {
                continue;
            }
            if r > 0 && grid[r - 1][c] > v {
                continue;
            }
            used[v] = true;
            grid[r][c] = v;
            total += go(lambda, grid, used, nr, nc);
            used[v] = false;
        }
        total
    }
    let n = size(lambda);
    let mut grid: Vec<Vec<usize>> = lambda.iter().map(|&len| vec![0; len]).collect();
    go(lambda, &mut grid, &mut vec![false; n + 1], 0, 0)
}

/// n! * det[1/(b_j - a_i)!] modulo P, with a_i = m-1-i and b_j = lambda_j + m-1-j,
/// by Gaussian elimination: the determinant of the first proof, evaluated numerically.
fn count_syt_determinant_mod(lambda: &[usize]) -> u64 {
    let m = lambda.len();
    let n = size(lambda);
    let max_arg = lambda[0] + m;
    let mut fact = vec![1u64; max_arg + 1];
    for k in 1..=max_arg {
        fact[k] = fact[k - 1] * k as u64 % P;
    }
    let inv_fact: Vec<u64> = fact.iter().map(|&f| mod_inv(f)).collect();
    let a: Vec<i64> = (0..m).map(|i| (m - 1 - i) as i64).collect();
    let b: Vec<i64> = (0..m).map(|j| (lambda[j] + m - 1 - j) as i64).collect();
    let mut mat: Vec<Vec<u64>> = (0..m)
        .map(|i| (0..m).map(|j| if b[j] >= a[i] { inv_fact[(b[j] - a[i]) as usize] } else { 0 }).collect())
        .collect();
    let mut det = 1u64;
    for col in 0..m {
        let pivot = (col..m).find(|&r| mat[r][col] != 0);
        let Some(pivot) = pivot else { return 0 };
        if pivot != col {
            mat.swap(pivot, col);
            det = (P - det) % P;
        }
        det = det * mat[col][col] % P;
        let inv = mod_inv(mat[col][col]);
        for r in col + 1..m {
            let factor = mat[r][col] * inv % P;
            for k in col..m {
                mat[r][k] = (mat[r][k] + P - factor * mat[col][k] % P) % P;
            }
        }
    }
    det * (1..=n as u64).fold(1u64, |acc, k| acc * k % P) % P
}

/// P(u, c) of the hook walk, as floating point, for every cell u and the
/// corner c = corner_cell: probability that a walk started at u ends at c.
fn walk_probabilities(lambda: &[usize], corner_cell: (usize, usize)) -> Vec<Vec<f64>> {
    let conj = conjugate(lambda);
    let hooks = hook_lengths(lambda);
    let mut prob: Vec<Vec<f64>> = lambda.iter().map(|&len| vec![0.0; len]).collect();
    for r in (0..lambda.len()).rev() {
        for c in (0..lambda[r]).rev() {
            if hooks[r][c] == 1 {
                prob[r][c] = if (r, c) == corner_cell { 1.0 } else { 0.0 };
                continue;
            }
            let mut sum = 0.0;
            for c2 in c + 1..lambda[r] {
                sum += prob[r][c2];
            }
            for r2 in r + 1..conj[c] {
                sum += prob[r2][c];
            }
            prob[r][c] = sum / (hooks[r][c] - 1) as f64;
        }
    }
    prob
}

struct Lcg(u64);
impl Lcg {
    fn below(&mut self, bound: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) % bound as u64) as usize
    }
}

fn main() {
    // --- A worked shape ---------------------------------------------------
    let lambda: Shape = vec![4, 2, 1];
    println!("hook lengths of (4,2,1): {:?}", hook_lengths(&lambda));
    println!("f(4,2,1) = {} (7! / 144)", count_syt_exact(&lambda));

    // --- Four independent counts agree on every partition of n <= 12 --------
    let mut memo: HashMap<Shape, u128> = HashMap::new();
    let mut all_ok = true;
    let mut shapes_checked = 0;
    for n in 0..=12 {
        for lam in partitions(n) {
            shapes_checked += 1;
            let rec = count_syt_recurrence(&lam, &mut memo);
            let hook = count_syt_exact(&lam);
            let hook_mod = count_syt_mod(&lam);
            let ok = rec == hook && hook_mod == (hook % P as u128) as u64;
            if !ok {
                all_ok = false;
                println!("MISMATCH (recurrence/hook) {:?}: {} vs {}", lam, rec, hook);
            }
            if !lam.is_empty() {
                let frob = count_syt_frobenius_mod(&lam);
                let det = count_syt_determinant_mod(&lam);
                if frob != hook_mod || det != hook_mod {
                    all_ok = false;
                    println!("MISMATCH (Frobenius/determinant) {:?}: {} {} vs {}", lam, frob, det, hook_mod);
                }
            }
        }
    }
    println!("recurrence, hook formula, Frobenius product and determinant agree on all {shapes_checked} partitions of n <= 12: {}", if all_ok { "PASS" } else { "FAIL" });

    // --- Definition-level brute force for n <= 8 ---------------------------
    let mut brute_ok = true;
    for n in 0..=8 {
        for lam in partitions(n) {
            if count_syt_brute_force(&lam) != count_syt_exact(&lam) {
                brute_ok = false;
                println!("BRUTE FORCE MISMATCH {:?}", lam);
            }
        }
    }
    println!("brute-force enumeration of tableaux matches the formula for all partitions of n <= 8: {}", if brute_ok { "PASS" } else { "FAIL" });

    // --- sum over lambda of f^2 = n! (a consequence of the RSK correspondence)
    let mut sq_ok = true;
    for n in 1..=12u128 {
        let total: u128 = partitions(n as usize).iter().map(|lam| count_syt_exact(lam).pow(2)).sum();
        if total != (1..=n).product::<u128>() {
            sq_ok = false;
        }
    }
    println!("sum of f^2 over partitions of n equals n! for n <= 12: {}", if sq_ok { "PASS" } else { "FAIL" });

    // --- Two-row rectangles give the Catalan numbers -------------------------
    let mut catalan_ok = true;
    let mut catalan = 1u128;
    for n in 0..=10u128 {
        if n > 0 {
            catalan = catalan * 2 * (2 * n - 1) / (n + 1);
        }
        let lam: Shape = if n == 0 { vec![] } else { vec![n as usize, n as usize] };
        if count_syt_exact(&lam) != catalan {
            catalan_ok = false;
        }
    }
    println!("f(n,n) equals the n-th Catalan number for n <= 10: {}", if catalan_ok { "PASS" } else { "FAIL" });

    // --- A large shape, modulo P ---------------------------------------------
    let big: Shape = vec![500, 500];
    let catalan_500 = {
        let fact = |k: u64| (1..=k).fold(1u64, |acc, x| acc * x % P);
        fact(1000) * mod_inv(fact(501) * fact(500) % P) % P
    };
    println!("f(500,500) mod p equals Catalan(500) mod p: {}", if count_syt_mod(&big) == catalan_500 && count_syt_frobenius_mod(&big) == catalan_500 { "PASS" } else { "FAIL" });

    // --- The row-hook lemma: hooks of row i together with the gaps b_i - b_j are exactly {1..b_i}
    let mut lemma_ok = true;
    for n in 1..=12 {
        for lam in partitions(n) {
            let m = lam.len();
            let b: Vec<usize> = (0..m).map(|i| lam[i] + m - 1 - i).collect();
            let hooks = hook_lengths(&lam);
            for i in 0..m {
                let mut values: Vec<usize> = hooks[i].clone();
                values.extend((i + 1..m).map(|j| b[i] - b[j]));
                values.sort();
                if values != (1..=b[i]).collect::<Vec<_>>() {
                    lemma_ok = false;
                }
            }
        }
    }
    println!("row-hook lemma ({{hooks of row i}} + {{b_i - b_j}} = {{1..b_i}}) holds for all partitions of n <= 12: {}", if lemma_ok { "PASS" } else { "FAIL" });

    // --- Hook walk reaches corner c with total probability H(lambda)/H(lambda - c)
    let mut walk_ok = true;
    for n in 1..=10 {
        for lam in partitions(n) {
            let h_lam: f64 = hook_lengths(&lam).iter().flatten().map(|&h| h as f64).product();
            for r in 0..lam.len() {
                if r + 1 == lam.len() || lam[r + 1] < lam[r] {
                    let corner = (r, lam[r] - 1);
                    let mut smaller = lam.clone();
                    smaller[r] -= 1;
                    if smaller[r] == 0 {
                        smaller.pop();
                    }
                    let h_small: f64 = hook_lengths(&smaller).iter().flatten().map(|&h| h as f64).product();
                    let total: f64 = walk_probabilities(&lam, corner).iter().flatten().sum();
                    if (total - h_lam / h_small).abs() > 1e-9 * (1.0 + total) {
                        walk_ok = false;
                    }
                }
            }
        }
    }
    println!("sum_u P(u, c) = H(lambda)/H(lambda - c) for every corner of every partition of n <= 10: {}", if walk_ok { "PASS" } else { "FAIL" });

    // --- The hook-walk sampler is uniform --------------------------------------
    let shape: Shape = vec![3, 2, 1];
    let f = count_syt_exact(&shape) as usize;
    let trials = 16_000usize;
    let mut rng = Lcg(2026);
    let mut counts: HashMap<Vec<Vec<usize>>, usize> = HashMap::new();
    for _ in 0..trials {
        let t = sample_tableau(&shape, &mut |k| rng.below(k));
        *counts.entry(t).or_insert(0) += 1;
    }
    let expected = trials as f64 / f as f64;
    let sigma = (expected * (1.0 - 1.0 / f as f64)).sqrt();
    let worst = counts.values().map(|&c| (c as f64 - expected).abs()).fold(0.0, f64::max);
    println!("hook-walk sampler on (3,2,1): {} distinct tableaux seen (f = {f}), every count within 4 sigma of uniform: {}", counts.len(), if counts.len() == f && worst < 4.0 * sigma { "PASS" } else { "FAIL" });
}
