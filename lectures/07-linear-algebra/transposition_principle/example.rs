// The transposition principle (Tellegen's principle): an algorithm for
// "y = Ax" built out of (1) loading x into some scratch variables,
// (2) a sequence of elementary steps z[a] += c * z[b], and (3) reading y
// out of some scratch variables, can be mechanically transformed into an
// algorithm for the TRANSPOSED problem "X = A^T Y" by reversing the step
// order, swapping each step's a/b roles, and swapping the roles of the
// load/store positions.
//
// This file implements a generic interpreter for such programs (forward
// and transposed), and: (1) validates the transposition principle itself
// against the matrix identity it is a restatement of, over many random
// programs; (2) applies it to the classic prefix-sum solution of range-sum
// queries to mechanically derive a range-add solution, cross-checked
// against both a hand-derived variant from the lecture and the well-known
// "imos method"; (3) applies it to the fast subset zeta transform to
// mechanically derive the fast superset zeta transform; and (4)
// reproduces the lecture's own small worked example exactly.

// PDFBEGIN
/// A "linear program": load x into some scratch slots, run a sequence of
/// z[a] += c*z[b] steps, then read y out of some scratch slots.
struct Program {
    k: usize,
    load_positions: Vec<usize>,  // length M: x[j] goes into z[load_positions[j]]
    ops: Vec<(usize, usize, i64)>, // (a, b, c): z[a] += c * z[b]
    store_positions: Vec<usize>, // length N: y[i] = z[store_positions[i]]
}

fn forward(prog: &Program, x: &[i64]) -> Vec<i64> {
    let mut z = vec![0i64; prog.k];
    for (j, &pos) in prog.load_positions.iter().enumerate() {
        // Additive, not overwriting: if two x_j share a load position, the
        // matrix-multiplication reading of "z = Sx" (S a 0/1 selection
        // matrix) sums their contributions rather than keeping only the
        // last one. The lecture's own examples never hit this case (every
        // load position it uses is distinct), but summing is what makes
        // the transpose identity hold in general, and it costs nothing
        // when positions are already distinct.
        z[pos] += x[j];
    }
    for &(a, b, c) in &prog.ops {
        z[a] += c * z[b];
    }
    prog.store_positions.iter().map(|&p| z[p]).collect()
}

/// The transposed algorithm from Theorem 2: reverse step order, swap each
/// step's (a,b) roles, and swap the roles of load/store positions.
fn transpose_run(prog: &Program, y: &[i64]) -> Vec<i64> {
    let mut z = vec![0i64; prog.k];
    for (i, &pos) in prog.store_positions.iter().enumerate() {
        // Additive for the same reason as in `forward`'s load step, mirrored
        // here since the transposed algorithm loads Y through T^T.
        z[pos] += y[i];
    }
    for &(a, b, c) in prog.ops.iter().rev() {
        z[b] += c * z[a];
    }
    prog.load_positions.iter().map(|&p| z[p]).collect()
}

// --- Section 5: range-sum queries (Problem 3A) and their transpose,
// range-add queries (Problem 3B), via the prefix-sum solution. ---

/// Builds a Program implementing Problem 3(A) exactly as in the lecture:
/// scratch variables s_0..s_n (indices 0..=n) and y_1..y_q (indices
/// n+1..n+q), computing y_i = sum_{j=L_i}^{R_i} x_j via prefix sums.
fn interval_sum_program(n: usize, queries: &[(usize, usize)]) -> Program {
    let q = queries.len();
    let k = (n + 1) + q;
    let load_positions: Vec<usize> = (0..n).map(|j| j + 1).collect();
    let mut ops = Vec::new();
    for j in 0..n {
        ops.push((j + 1, j, 1));
    }
    for (i, &(l, r)) in queries.iter().enumerate() {
        let yi = n + 1 + i;
        ops.push((yi, r + 1, 1));
        ops.push((yi, l, -1));
    }
    let store_positions: Vec<usize> = (0..q).map(|i| n + 1 + i).collect();
    Program { k, load_positions, ops, store_positions }
}

/// The hand-derived transpose from section 5.3 of the lecture.
fn transpose_hand_derived(n: usize, queries: &[(usize, usize)], y: &[i64]) -> Vec<i64> {
    let mut s = vec![0i64; n + 1];
    for (i, &(l, r)) in queries.iter().enumerate() {
        s[l] -= y[i];
        s[r + 1] += y[i];
    }
    for j in (0..n).rev() {
        s[j] += s[j + 1];
    }
    (0..n).map(|j| s[j + 1]).collect()
}

/// The "imos method" variant from section 5.4 (different loop direction
/// and sign convention, same result).
fn transpose_imos(n: usize, queries: &[(usize, usize)], y: &[i64]) -> Vec<i64> {
    let mut s = vec![0i64; n + 1];
    for (i, &(l, r)) in queries.iter().enumerate() {
        s[l] += y[i];
        s[r + 1] -= y[i];
    }
    for j in 0..n {
        s[j + 1] += s[j];
    }
    (0..n).map(|j| s[j]).collect()
}

// --- Section 6: subset zeta transform (Problem 4A) and its transpose,
// the superset zeta transform (Problem 4B). ---

fn subset_zeta_program(n: usize) -> Program {
    let size = 1usize << n;
    let load_positions: Vec<usize> = (0..size).collect();
    let mut ops = Vec::new();
    for i in 0..n {
        for t in 0..size {
            if (t >> i) & 1 == 1 {
                ops.push((t, t ^ (1 << i), 1));
            }
        }
    }
    let store_positions: Vec<usize> = (0..size).collect();
    Program { k: size, load_positions, ops, store_positions }
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

/// Computes the N x M matrix A that `forward` implements, by running it on
/// each standard basis vector. Used only to validate the transposition
/// principle itself against A^T, not part of either algorithm.
fn compute_matrix(prog: &Program) -> Vec<Vec<i64>> {
    let m = prog.load_positions.len();
    let n = prog.store_positions.len();
    let mut a = vec![vec![0i64; m]; n];
    for j in 0..m {
        let mut x = vec![0i64; m];
        x[j] = 1;
        let y = forward(prog, &x);
        for i in 0..n {
            a[i][j] = y[i];
        }
    }
    a
}

fn mat_transpose_vec(a: &[Vec<i64>], y: &[i64]) -> Vec<i64> {
    let n = a.len();
    let m = if n > 0 { a[0].len() } else { 0 };
    (0..m).map(|j| (0..n).fold(0i64, |acc, i| acc + a[i][j] * y[i])).collect()
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

fn random_program(rng: &mut Rng) -> Program {
    let k = rng.range(2, 8) as usize;
    let m = rng.range(1, k as i64) as usize;
    let n = rng.range(1, k as i64) as usize;
    let q = rng.range(0, 15) as usize;
    let load_positions: Vec<usize> = (0..m).map(|_| rng.range(0, k as i64 - 1) as usize).collect();
    let store_positions: Vec<usize> = (0..n).map(|_| rng.range(0, k as i64 - 1) as usize).collect();
    let ops: Vec<(usize, usize, i64)> = (0..q)
        .map(|_| {
            let a = rng.range(0, k as i64 - 1) as usize;
            let b = rng.range(0, k as i64 - 1) as usize;
            let c = rng.range(-5, 5);
            (a, b, c)
        })
        .collect();
    Program { k, load_positions, ops, store_positions }
}

fn interval_sum_naive(n: usize, queries: &[(usize, usize)], x: &[i64]) -> Vec<i64> {
    let _ = n;
    queries.iter().map(|&(l, r)| x[l..=r].iter().sum()).collect()
}

fn interval_add_naive(n: usize, queries: &[(usize, usize)], y: &[i64]) -> Vec<i64> {
    (0..n)
        .map(|j| {
            queries
                .iter()
                .zip(y.iter())
                .filter(|&(&(l, r), _)| l <= j && j <= r)
                .map(|(_, &yi)| yi)
                .sum()
        })
        .collect()
}

fn subset_zeta_naive(n: usize, a: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    (0..size).map(|t| (0..size).filter(|&s| s & t == s).map(|s| a[s]).sum()).collect()
}

fn superset_zeta_naive(n: usize, a: &[i64]) -> Vec<i64> {
    let size = 1usize << n;
    (0..size).map(|t| (0..size).filter(|&s| s & t == t).map(|s| a[s]).sum()).collect()
}

fn main() {
    let mut rng = Rng(2026092001);

    // --- Theorem 2 itself: transpose_run(Y) == A^T Y, for many random programs ---
    let mut theorem_ok = true;
    for _ in 0..300 {
        let prog = random_program(&mut rng);
        let n = prog.store_positions.len();
        let a = compute_matrix(&prog);
        let y: Vec<i64> = (0..n).map(|_| rng.range(-10, 10)).collect();
        let via_transpose = transpose_run(&prog, &y);
        let via_matrix = mat_transpose_vec(&a, &y);
        if via_transpose != via_matrix {
            theorem_ok = false;
            println!("TRANSPOSITION PRINCIPLE MISMATCH");
        }
    }
    println!("Transposition principle (transpose_run(Y) == A^T Y) over 300 random programs: {}", if theorem_ok { "PASS" } else { "FAIL" });

    // --- Worked example (section 3.3): y1=2x1+3x2, y2=4x1+5x2, y3=6x1+7x2 ---
    let example = Program {
        k: 5,
        load_positions: vec![0, 1],
        ops: vec![(2, 0, 2), (2, 1, 3), (3, 0, 4), (3, 1, 5), (4, 0, 6), (4, 1, 7)],
        store_positions: vec![2, 3, 4],
    };
    let mut example_ok = true;
    for _ in 0..50 {
        let x1 = rng.range(-20, 20);
        let x2 = rng.range(-20, 20);
        let y = forward(&example, &[x1, x2]);
        let expected_y = [2 * x1 + 3 * x2, 4 * x1 + 5 * x2, 6 * x1 + 7 * x2];
        let y1 = rng.range(-20, 20);
        let y2 = rng.range(-20, 20);
        let y3 = rng.range(-20, 20);
        let x_t = transpose_run(&example, &[y1, y2, y3]);
        let expected_x = [2 * y1 + 4 * y2 + 6 * y3, 3 * y1 + 5 * y2 + 7 * y3];
        if y != expected_y || x_t != expected_x {
            example_ok = false;
            println!("WORKED EXAMPLE MISMATCH");
        }
    }
    println!("Worked example (section 3.3): forward and transpose match the lecture's own formulas, over 50 random trials: {}", if example_ok { "PASS" } else { "FAIL" });

    // --- Section 5: range-sum queries and their transpose ---
    let mut interval_ok = true;
    for _ in 0..100 {
        let n = rng.range(1, 20) as usize;
        let q = rng.range(0, 15) as usize;
        let queries: Vec<(usize, usize)> = (0..q)
            .map(|_| {
                let l = rng.range(0, n as i64 - 1) as usize;
                let r = rng.range(l as i64, n as i64 - 1) as usize;
                (l, r)
            })
            .collect();
        let prog = interval_sum_program(n, &queries);

        let x: Vec<i64> = (0..n).map(|_| rng.range(-20, 20)).collect();
        let y_via_prog = forward(&prog, &x);
        let y_via_naive = interval_sum_naive(n, &queries, &x);
        if y_via_prog != y_via_naive {
            interval_ok = false;
            println!("INTERVAL SUM MISMATCH n={n} q={q}");
        }

        let y: Vec<i64> = (0..q).map(|_| rng.range(-20, 20)).collect();
        let x_via_transpose = transpose_run(&prog, &y);
        let x_via_naive = interval_add_naive(n, &queries, &y);
        let x_via_hand = transpose_hand_derived(n, &queries, &y);
        let x_via_imos = transpose_imos(n, &queries, &y);
        if x_via_transpose != x_via_naive || x_via_transpose != x_via_hand || x_via_transpose != x_via_imos {
            interval_ok = false;
            println!("INTERVAL ADD MISMATCH n={n} q={q}");
        }
    }
    println!("Range-sum queries (Problem 3A) match brute force, over 100 random trials: {}", if interval_ok { "PASS" } else { "FAIL" });
    println!("Mechanically-derived transpose == hand-derived (5.3) == imos method (5.4) == brute-force range-add, over 100 random trials: {}", if interval_ok { "PASS" } else { "FAIL" });

    // --- Section 6: subset zeta transform and its transpose ---
    let mut zeta_ok = true;
    for _ in 0..30 {
        let n = rng.range(0, 6) as usize;
        let size = 1usize << n;
        let prog = subset_zeta_program(n);

        let x: Vec<i64> = (0..size).map(|_| rng.range(-10, 10)).collect();
        let y_via_prog = forward(&prog, &x);
        let y_via_naive = subset_zeta_naive(n, &x);
        if y_via_prog != y_via_naive {
            zeta_ok = false;
            println!("SUBSET ZETA MISMATCH n={n}");
        }

        let y: Vec<i64> = (0..size).map(|_| rng.range(-10, 10)).collect();
        let x_via_transpose = transpose_run(&prog, &y);
        let x_via_naive = superset_zeta_naive(n, &y);
        if x_via_transpose != x_via_naive {
            zeta_ok = false;
            println!("SUPERSET ZETA (TRANSPOSE) MISMATCH n={n}");
        }
    }
    println!("Subset zeta transform (Problem 4A) and its mechanically-derived transpose, the superset zeta transform (Problem 4B), over 30 random trials: {}", if zeta_ok { "PASS" } else { "FAIL" });

    let all_ok = theorem_ok && example_ok && interval_ok && zeta_ok;
    println!("\nAll checks: {}", if all_ok { "PASS" } else { "FAIL" });
}
