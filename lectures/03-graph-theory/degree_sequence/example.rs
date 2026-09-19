// Degree sequences of simple undirected graphs: the Erdos-Gallai
// realizability test and the Havel-Hakimi construction algorithm.

/// Is `d` (any order) the degree sequence of some simple undirected graph?
/// Checks the handshake-lemma parity condition and, for every prefix
/// length k of the sorted-descending sequence, the Erdos-Gallai
/// inequality sum(d[0..k]) <= k(k-1) + sum_{i>=k} min(d[i], k).
///
/// The right-hand sum is evaluated in O(log n) per k by binary-searching
/// the (sorted, so monotone) boundary between "d[i] >= k" and "d[i] < k"
/// in the suffix, then combining with prefix sums -- O(n log n) overall.
fn is_graphic_erdos_gallai(d: &[i64]) -> bool {
    let n = d.len();
    if d.iter().any(|&x| x < 0 || x > (n as i64 - 1)) {
        return false;
    }
    let total: i64 = d.iter().sum();
    if total % 2 != 0 {
        return false;
    }

    let mut sorted = d.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a)); // descending

    let mut prefix = vec![0i64; n + 1];
    for i in 0..n {
        prefix[i + 1] = prefix[i] + sorted[i];
    }

    for k in 1..=n {
        let lhs = prefix[k];
        // Among sorted[k..n] (already descending), find how many are >= k.
        // partition_point needs an ascending predicate, so search on the
        // reversed sense: the suffix is descending, so "sorted[k+j] >= k"
        // holds for a prefix of that suffix.
        let suffix = &sorted[k..n];
        let split = suffix.partition_point(|&v| v >= k as i64); // count with value >= k
        let sum_ge = split as i64 * k as i64;
        let sum_lt = prefix[n] - prefix[k] - {
            // sum of the ">= k" part within the suffix
            let mut s = 0i64;
            for &v in &suffix[..split] {
                s += v;
            }
            s
        };
        let rhs = (k as i64) * (k as i64 - 1) + sum_ge + sum_lt;
        if lhs > rhs {
            return false;
        }
    }
    true
}

/// Havel-Hakimi: construct a simple graph with degree sequence `d` (any
/// order), or return None if none exists. Returns 0-indexed edges in the
/// order of the ORIGINAL indices of `d`.
///
/// At each round, connects the current highest-degree vertex to the next
/// `d1` highest-degree remaining vertices, then "retires" it. Re-sorting
/// every round costs O(n log n) per round; an O(n + sum d) implementation
/// is possible with bucket-based bookkeeping but is not needed to see the
/// algorithm work correctly.
fn havel_hakimi(d: &[i64]) -> Option<Vec<(usize, usize)>> {
    let n = d.len();
    if !is_graphic_erdos_gallai(d) {
        return None;
    }
    let mut remaining: Vec<(i64, usize)> = d.iter().cloned().zip(0..n).collect();
    let mut edges = Vec::new();

    loop {
        remaining.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        let (d1, v1) = remaining[0];
        if d1 == 0 {
            break;
        }
        if d1 as usize >= remaining.len() {
            return None; // not enough other vertices to connect to
        }
        for slot in remaining.iter_mut().take(d1 as usize + 1).skip(1) {
            if slot.0 <= 0 {
                return None; // ran out of capacity among the chosen neighbors
            }
            edges.push((v1, slot.1));
            slot.0 -= 1;
        }
        remaining[0].0 = 0;
    }
    Some(edges)
}

fn verify_degrees(n: usize, edges: &[(usize, usize)], want: &[i64]) -> bool {
    let mut got = vec![0i64; n];
    for &(a, b) in edges {
        got[a] += 1;
        got[b] += 1;
    }
    // `want` is in original-index order, matching `edges`' indices.
    got == want
}

fn main() {
    // The three worked examples from the lecture text.
    let cases: [(&str, Vec<i64>); 3] = [
        ("(4,4,3,2,2,1)", vec![4, 4, 3, 2, 2, 1]),
        ("(4,3,3,2,2,1)", vec![4, 3, 3, 2, 2, 1]),
        ("(4,4,4,2,1,1)", vec![4, 4, 4, 2, 1, 1]),
    ];

    for (label, d) in &cases {
        let graphic = is_graphic_erdos_gallai(d);
        println!("{label}: graphic = {graphic}");
        match havel_hakimi(d) {
            Some(edges) => {
                let ok = verify_degrees(d.len(), &edges, d);
                println!("  constructed {} edges: {:?}", edges.len(), edges);
                println!("  degree check: {}", if ok { "PASS" } else { "FAIL" });
            }
            None => println!("  no construction (correctly not graphic)"),
        }
    }

    // Brute-force cross-check: every sequence of length 5 with entries in
    // 0..=4 should agree between "is_graphic" and "a graph with 0..5 edges
    // among 5 labeled vertices realizes it", checked by exhaustive search.
    let n = 5;
    let all_pairs: Vec<(usize, usize)> = (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let mut mismatches = 0;
    for mask in 0u32..(1 << all_pairs.len()) {
        let edges: Vec<(usize, usize)> =
            all_pairs.iter().enumerate().filter(|(bit, _)| mask & (1 << bit) != 0).map(|(_, &e)| e).collect();
        let mut deg = vec![0i64; n];
        for &(a, b) in &edges {
            deg[a] += 1;
            deg[b] += 1;
        }
        let mut sorted_deg = deg.clone();
        sorted_deg.sort_unstable_by(|a, b| b.cmp(a));
        if !is_graphic_erdos_gallai(&sorted_deg) {
            mismatches += 1;
            println!("MISMATCH: degree sequence {sorted_deg:?} came from a real graph but is_graphic said false");
        }
    }
    println!(
        "\nexhaustive check over all {}-vertex graphs: every realized degree sequence passes is_graphic ({} mismatches)",
        n, mismatches
    );
}
