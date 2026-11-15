// Parking functions: a length-n sequence a_1..a_n with values in [1,n] is a
// parking function if, sending car i to look for the first free spot at or
// after a_i (spots numbered 1..n), every car finds a spot. There are
// exactly (n+1)^(n-1) of them, proved here via a circular extension trick.

// PDFBEGIN
/// Definition 1, directly: simulate the parking procedure.
fn is_parking_function_sim(a: &[usize], n: usize) -> bool {
    let mut occupied = vec![false; n + 1]; // 1-indexed spots 1..=n
    for &ai in a {
        let mut parked = false;
        for j in ai..=n {
            if !occupied[j] {
                occupied[j] = true;
                parked = true;
                break;
            }
        }
        if !parked {
            return false;
        }
    }
    true
}

/// Proposition 3, condition 2: for every k, at least k of the a_i are <= k.
fn is_parking_function_condition(a: &[usize], n: usize) -> bool {
    (1..=n).all(|k| a.iter().filter(|&&ai| ai <= k).count() >= k)
}

/// The circular extension from the counting proof: n cars, spots 0..=n
/// (n+1 of them), each a_i in 0..=n, searching circularly from a_i. This
/// never fails; returns the one spot left empty.
fn circular_park_empty_slot(a: &[usize], n: usize) -> usize {
    let m = n + 1;
    let mut occupied = vec![false; m];
    for &ai in a {
        let mut j = ai;
        loop {
            if !occupied[j] {
                occupied[j] = true;
                break;
            }
            j = (j + 1) % m;
        }
    }
    (0..m).find(|&s| !occupied[s]).unwrap()
}

/// Proposition 5's generalized condition, for length-m sequences with
/// values in [1,n]: for every k, at least (m-n+k) of the a_i are <= k.
fn satisfies_generalized_condition(a: &[usize], n: usize, m: usize) -> bool {
    (1..=n).all(|k| {
        let need = m as i64 - n as i64 + k as i64;
        let have = a.iter().filter(|&&ai| ai <= k).count() as i64;
        have >= need
    })
}

fn binom_pow(base: u64, exp: u32) -> u64 {
    base.pow(exp)
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn main() {
    // --- The source's own worked examples (Problem 2) --------------------
    let cases: [(&[usize], bool); 4] = [
        (&[4, 1, 5, 1, 3], true),
        (&[5, 3, 1, 4, 3], false),
        (&[4, 1, 6, 5, 4, 1], false),
        (&[3, 6, 1, 5, 3, 1], true),
    ];
    let mut worked_ok = true;
    for (a, expected) in cases {
        let n = a.len();
        let sim = is_parking_function_sim(a, n);
        let cond = is_parking_function_condition(a, n);
        if sim != expected || cond != expected {
            worked_ok = false;
            println!("WORKED-EXAMPLE MISMATCH a={:?}: sim={sim} cond={cond} expected={expected}", a);
        }
    }
    println!("lecture's own worked examples (Problem 2) match: {}", if worked_ok { "PASS" } else { "FAIL" });

    // --- Equivalence of the simulation and the threshold condition -------
    // Brute force over all a in {1,...,n}^n for small n.
    let mut equiv_ok = true;
    let mut count_ok = true;
    for n in 1..=6usize {
        let mut count = 0u64;
        let mut a = vec![1usize; n];
        loop {
            let sim = is_parking_function_sim(&a, n);
            let cond = is_parking_function_condition(&a, n);
            if sim != cond {
                equiv_ok = false;
                println!("EQUIV MISMATCH n={n} a={:?}: sim={sim} cond={cond}", a);
            }
            if sim {
                count += 1;
            }
            // next tuple in {1..=n}^n, odometer style
            let mut i = n;
            loop {
                if i == 0 {
                    break;
                }
                i -= 1;
                a[i] += 1;
                if a[i] <= n {
                    break;
                }
                a[i] = 1;
                if i == 0 {
                    a = vec![n + 1; n]; // sentinel: overflow past the last tuple
                    break;
                }
            }
            if a[0] > n {
                break;
            }
        }
        let expected = binom_pow(n as u64 + 1, n as u32 - 1);
        if count != expected {
            count_ok = false;
            println!("COUNT MISMATCH n={n}: got {count} expected {expected}");
        }
    }
    println!("simulation matches threshold condition, for n=1..=6: {}", if equiv_ok { "PASS" } else { "FAIL" });
    println!("parking function count matches (n+1)^(n-1), for n=1..=6: {}", if count_ok { "PASS" } else { "FAIL" });

    // --- Circular extension: bijection + symmetry across slots -----------
    let mut bijection_ok = true;
    for n in 1..=6usize {
        let mut a = vec![1usize; n]; // values in 1..=n (never 0)
        loop {
            let leaves_zero_empty = circular_park_empty_slot(&a, n) == 0;
            let is_pf = is_parking_function_sim(&a, n);
            if leaves_zero_empty != is_pf {
                bijection_ok = false;
                println!("BIJECTION MISMATCH n={n} a={:?}: circular={leaves_zero_empty} linear={is_pf}", a);
            }
            let mut i = n;
            let mut done = false;
            loop {
                if i == 0 {
                    done = true;
                    break;
                }
                i -= 1;
                a[i] += 1;
                if a[i] <= n {
                    break;
                }
                a[i] = 1;
            }
            if done {
                break;
            }
        }
    }
    println!("circular parking leaves spot 0 empty <=> a is an ordinary parking function, for n=1..=6: {}", if bijection_ok { "PASS" } else { "FAIL" });

    let mut symmetry_ok = true;
    for n in 1..=5usize {
        let m = n + 1;
        let mut tally = vec![0u64; m];
        let mut a = vec![0usize; n]; // values in 0..=n
        loop {
            let empty = circular_park_empty_slot(&a, n);
            tally[empty] += 1;
            let mut i = n;
            let mut done = false;
            loop {
                if i == 0 {
                    done = true;
                    break;
                }
                i -= 1;
                a[i] += 1;
                if a[i] <= n {
                    break;
                }
                a[i] = 0;
            }
            if done {
                break;
            }
        }
        let expected = binom_pow(n as u64 + 1, n as u32 - 1);
        if tally.iter().any(|&c| c != expected) {
            symmetry_ok = false;
            println!("SYMMETRY MISMATCH n={n}: tally={:?} expected each={expected}", tally);
        }
    }
    println!("every spot is left empty equally often ((n+1)^(n-1) times), for n=1..=5: {}", if symmetry_ok { "PASS" } else { "FAIL" });

    // --- Proposition 5: the generalized (length-m, values in [1,n]) count -
    let mut prop5_ok = true;
    for n in 1..=6usize {
        for m in 1..=n {
            let mut count = 0u64;
            let mut a = vec![1usize; m];
            loop {
                if satisfies_generalized_condition(&a, n, m) {
                    count += 1;
                }
                let mut i = m;
                let mut done = false;
                loop {
                    if i == 0 {
                        done = true;
                        break;
                    }
                    i -= 1;
                    a[i] += 1;
                    if a[i] <= n {
                        break;
                    }
                    a[i] = 1;
                }
                if done {
                    break;
                }
            }
            let expected = (n as u64 + 1 - m as u64) * binom_pow(n as u64 + 1, m as u32 - 1);
            if count != expected {
                prop5_ok = false;
                println!("PROP5 MISMATCH n={n} m={m}: got {count} expected {expected}");
            }
        }
    }
    println!("Proposition 5 count (n+1-m)(n+1)^(m-1) matches brute force, for n,m<=6: {}", if prop5_ok { "PASS" } else { "FAIL" });
}
