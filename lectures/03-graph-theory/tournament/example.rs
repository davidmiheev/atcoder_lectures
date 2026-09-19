// Tournaments: a simple directed graph where every pair of distinct
// vertices has exactly one directed edge between them (the result of a
// round-robin contest). Demonstrates: transitivity testing, strongly-
// connected-component decomposition purely from in-degrees, a Hamiltonian
// path (Redei's theorem), a Hamiltonian cycle for strongly connected
// tournaments (Camion's theorem), and a 2-king (Landau's theorem).

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn bool(&mut self) -> bool {
        self.next() & 1 == 1
    }
}

/// beats[i][j] = true iff i -> j (i beat j). Exactly one of beats[i][j],
/// beats[j][i] is true for every i != j.
fn random_tournament(n: usize, rng: &mut Rng) -> Vec<Vec<bool>> {
    let mut beats = vec![vec![false; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if rng.bool() {
                beats[i][j] = true;
            } else {
                beats[j][i] = true;
            }
        }
    }
    beats
}

fn is_transitive(beats: &[Vec<bool>], n: usize) -> bool {
    // Equivalent to "contains no 3-cycle".
    for a in 0..n {
        for b in 0..n {
            if !beats[a][b] {
                continue;
            }
            for c in 0..n {
                if beats[b][c] && beats[c][a] {
                    return false;
                }
            }
        }
    }
    true
}

fn in_degrees(beats: &[Vec<bool>], n: usize) -> Vec<usize> {
    (0..n).map(|v| (0..n).filter(|&u| beats[u][v]).count()).collect()
}

/// Reference SCC decomposition via mutual reachability (O(n^3) transitive
/// closure), used only to cross-check the in-degree-based method below.
fn scc_by_reachability(beats: &[Vec<bool>], n: usize) -> Vec<Vec<usize>> {
    let mut reach = vec![vec![false; n]; n];
    for i in 0..n {
        reach[i][i] = true;
        for j in 0..n {
            if beats[i][j] {
                reach[i][j] = true;
            }
        }
    }
    for k in 0..n {
        for i in 0..n {
            if reach[i][k] {
                for j in 0..n {
                    if reach[k][j] {
                        reach[i][j] = true;
                    }
                }
            }
        }
    }
    let mut assigned = vec![false; n];
    let mut comps = Vec::new();
    for v in 0..n {
        if assigned[v] {
            continue;
        }
        let comp: Vec<usize> = (0..n).filter(|&u| reach[v][u] && reach[u][v]).collect();
        for &u in &comp {
            assigned[u] = true;
        }
        comps.push(comp);
    }
    comps
}

/// SCC decomposition using ONLY in-degrees: sort vertices ascending by
/// in-degree, then cut wherever the running sum of in-degrees among the
/// first k equals C(k,2) -- the signature of "no edges enter from outside".
fn scc_by_indegree(beats: &[Vec<bool>], n: usize) -> Vec<Vec<usize>> {
    let indeg = in_degrees(beats, n);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&v| indeg[v]);

    let mut comps = Vec::new();
    let mut start = 0usize;
    let mut running = 0i64;
    for k in 1..=n {
        running += indeg[order[k - 1]] as i64;
        let target = (k as i64) * (k as i64 - 1) / 2;
        if running == target {
            comps.push(order[start..k].to_vec());
            start = k;
        }
    }
    comps
}

/// Redei's theorem: build a Hamiltonian path by inserting vertices one at
/// a time into the correct position of the path built so far.
fn hamiltonian_path(beats: &[Vec<bool>], n: usize) -> Vec<usize> {
    let mut path: Vec<usize> = Vec::with_capacity(n);
    for v in 0..n {
        if path.is_empty() {
            path.push(v);
            continue;
        }
        if beats[v][path[0]] {
            path.insert(0, v);
        } else if beats[*path.last().unwrap()][v] {
            path.push(v);
        } else {
            // find i with path[i] -> v -> path[i+1]
            let mut idx = 0;
            for i in 0..path.len() - 1 {
                if beats[path[i]][v] && beats[v][path[i + 1]] {
                    idx = i + 1;
                    break;
                }
            }
            path.insert(idx, v);
        }
    }
    path
}

fn is_valid_hamiltonian_path(beats: &[Vec<bool>], path: &[usize], n: usize) -> bool {
    if path.len() != n {
        return false;
    }
    let mut seen = vec![false; n];
    for &v in path {
        if seen[v] {
            return false;
        }
        seen[v] = true;
    }
    path.windows(2).all(|w| beats[w[0]][w[1]])
}

/// Camion's theorem: for a strongly connected tournament with n >= 3,
/// extend a Hamiltonian path's growing prefix cycle one step at a time
/// until it covers everything.
fn hamiltonian_cycle(beats: &[Vec<bool>], n: usize, path: &[usize]) -> Vec<usize> {
    let mut order = path.to_vec();
    // Find a starting cycle length k with order[k-1] -> order[0].
    let mut k = 0usize;
    for j in (1..n).rev() {
        if beats[order[j]][order[0]] {
            k = j + 1;
            break;
        }
    }
    while k < n {
        let next = order[k]; // this is v_{k+1}
        if beats[next][order[0]] {
            // Case A: extends directly.
            k += 1;
            continue;
        }
        // Search for minimal p0 in 0..k with next -> order[p0].
        let mut found_p: Option<usize> = None;
        for p0 in 0..k {
            if beats[next][order[p0]] {
                found_p = Some(p0);
                break;
            }
        }
        if let Some(p0) = found_p {
            // p0 >= 1 always here (p0 == 0 would have been caught above).
            let mut new_order = Vec::with_capacity(n);
            new_order.extend_from_slice(&order[p0..k]);
            new_order.extend_from_slice(&order[0..p0]);
            new_order.push(next);
            new_order.extend_from_slice(&order[k + 1..n]);
            order = new_order;
            k += 1;
        } else {
            // Case C: every current cycle vertex beats `next`. Use strong
            // connectivity to find the first later vertex with an edge
            // back into the cycle.
            let mut found: Option<(usize, usize)> = None;
            'outer: for q0 in (k + 1)..n {
                for p0 in 0..k {
                    if beats[order[q0]][order[p0]] {
                        found = Some((q0, p0));
                        break 'outer;
                    }
                }
            }
            let (q0, p0) = found.expect("strongly connected tournament guarantees this exists");
            if p0 == 0 {
                // order[q0] -> order[0] directly: no reordering needed.
                k = q0 + 1;
            } else {
                let mut new_order = Vec::with_capacity(n);
                new_order.extend_from_slice(&order[p0..k]);
                new_order.extend_from_slice(&order[0..p0]);
                new_order.extend_from_slice(&order[k..=q0]);
                new_order.extend_from_slice(&order[q0 + 1..n]);
                order = new_order;
                k = q0 + 1;
            }
        }
    }
    order
}

fn is_valid_hamiltonian_cycle(beats: &[Vec<bool>], cycle: &[usize], n: usize) -> bool {
    if !is_valid_hamiltonian_path(beats, cycle, n) {
        return false;
    }
    beats[cycle[n - 1]][cycle[0]]
}

/// Landau's theorem: the maximum out-degree vertex reaches every other
/// vertex within distance 2.
fn find_2_king(beats: &[Vec<bool>], n: usize) -> usize {
    (0..n).max_by_key(|&v| (0..n).filter(|&u| beats[v][u]).count()).unwrap()
}

fn is_2_king(beats: &[Vec<bool>], v: usize, n: usize) -> bool {
    (0..n).all(|x| {
        x == v || beats[v][x] || (0..n).any(|y| beats[v][y] && beats[y][x])
    })
}

fn main() {
    let mut rng = Rng(12345);

    // --- Transitivity, SCC decomposition, in-degree cross-check --------
    let mut scc_ok = true;
    for trial in 0..200 {
        let n = 2 + (trial % 6);
        let beats = random_tournament(n, &mut rng);
        let mut ref_comps = scc_by_reachability(&beats, n);
        let mut deg_comps = scc_by_indegree(&beats, n);
        for c in ref_comps.iter_mut() {
            c.sort();
        }
        for c in deg_comps.iter_mut() {
            c.sort();
        }
        let mut ref_set: Vec<Vec<usize>> = ref_comps;
        let mut deg_set: Vec<Vec<usize>> = deg_comps;
        ref_set.sort();
        deg_set.sort();
        if ref_set != deg_set {
            scc_ok = false;
            println!("SCC MISMATCH n={n}: {:?} vs {:?}", ref_set, deg_set);
        }
    }
    println!("in-degree-only SCC decomposition matches reachability-based SCC: {}", if scc_ok { "PASS" } else { "FAIL" });

    // --- Hamiltonian path, on many random tournaments -------------------
    let mut path_ok = true;
    for trial in 0..200 {
        let n = 1 + (trial % 8);
        let beats = random_tournament(n, &mut rng);
        let path = hamiltonian_path(&beats, n);
        if !is_valid_hamiltonian_path(&beats, &path, n) {
            path_ok = false;
            println!("HAM PATH FAIL n={n} beats={:?} path={:?}", beats, path);
        }
    }
    println!("Redei Hamiltonian path valid on every trial: {}", if path_ok { "PASS" } else { "FAIL" });

    // --- Hamiltonian cycle, on strongly connected tournaments ------------
    let mut cycle_ok = true;
    let mut cycle_trials = 0;
    let mut trial = 0;
    while cycle_trials < 200 {
        trial += 1;
        let n = 3 + (trial % 6);
        let beats = random_tournament(n, &mut rng);
        if scc_by_reachability(&beats, n).len() != 1 {
            continue; // not strongly connected; skip
        }
        cycle_trials += 1;
        let path = hamiltonian_path(&beats, n);
        let cycle = hamiltonian_cycle(&beats, n, &path);
        if !is_valid_hamiltonian_cycle(&beats, &cycle, n) {
            cycle_ok = false;
            println!("HAM CYCLE FAIL n={n} beats={:?} cycle={:?}", beats, cycle);
        }
    }
    println!("Camion Hamiltonian cycle valid on every strongly-connected trial: {}", if cycle_ok { "PASS" } else { "FAIL" });

    // --- 2-king ------------------------------------------------------------
    let mut king_ok = true;
    for trial in 0..200 {
        let n = 1 + (trial % 8);
        let beats = random_tournament(n, &mut rng);
        let king = find_2_king(&beats, n);
        if !is_2_king(&beats, king, n) {
            king_ok = false;
            println!("2-KING FAIL n={n} beats={:?} king={king}", beats);
        }
    }
    println!("max-out-degree vertex is always a 2-king: {}", if king_ok { "PASS" } else { "FAIL" });

    // --- A concrete worked example: rock-paper-scissors plus two more ----
    // 0=Rock 1=Scissors 2=Paper (Rock beats Scissors, Scissors beats Paper,
    // Paper beats Rock): a 3-cycle, so NOT transitive.
    let n = 3;
    let mut beats = vec![vec![false; n]; n];
    beats[0][1] = true; // Rock beats Scissors
    beats[1][2] = true; // Scissors beats Paper
    beats[2][0] = true; // Paper beats Rock
    println!("\nrock-paper-scissors tournament transitive? {}", is_transitive(&beats, n));
    println!("its SCCs (in-degree method): {:?}", scc_by_indegree(&beats, n));
}
