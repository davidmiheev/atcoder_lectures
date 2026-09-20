// Sperner's theorem: the largest antichain in the Boolean lattice B_n (all
// subsets of {0,...,n-1}, ordered by inclusion) has exactly C(n, floor(n/2))
// elements. Demonstrates both proofs from the lecture -- the symmetric
// chain decomposition (constructive) and the LYM inequality (probabilistic)
// -- plus an independent cross-check via Dilworth's theorem (max antichain
// = min chain cover, computed by bipartite matching).

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

fn popcount(x: u32) -> u32 {
    x.count_ones()
}

// ---------------------------------------------------------------------
// Proof 2: symmetric chain decomposition
// ---------------------------------------------------------------------

/// Partition B_n into disjoint symmetric chains, via the induction in the
/// lecture's proof: extend each chain of B_{n-1} by adding bit (n-1) at
/// the top, and peel off a length-1 leftover chain from the middle when
/// the parent chain had length > 1.
fn symmetric_chain_decomposition(n: usize) -> Vec<Vec<u32>> {
    if n == 0 {
        return vec![vec![0u32]];
    }
    let prev = symmetric_chain_decomposition(n - 1);
    let bit = 1u32 << (n - 1);
    let mut result = Vec::new();
    for chain in prev {
        let m = chain.len();
        if m == 1 {
            result.push(vec![chain[0], chain[0] | bit]);
        } else {
            let mut chain1 = chain.clone();
            chain1.push(chain[m - 1] | bit);
            let chain2: Vec<u32> = chain[0..m - 1].iter().map(|&s| s | bit).collect();
            result.push(chain1);
            result.push(chain2);
        }
    }
    result
}

fn verify_symmetric_chain_decomposition(n: usize, chains: &[Vec<u32>]) -> bool {
    let total = 1usize << n;
    let mut seen = vec![false; total];
    for chain in chains {
        // strictly increasing by exactly one element each step
        for w in chain.windows(2) {
            if w[0] & w[1] != w[0] || popcount(w[1]) != popcount(w[0]) + 1 {
                return false;
            }
        }
        // symmetric: min size + max size = n
        let lo = popcount(chain[0]);
        let hi = popcount(*chain.last().unwrap());
        if lo + hi != n as u32 {
            return false;
        }
        for &s in chain {
            let idx = s as usize;
            if seen[idx] {
                return false; // covered twice
            }
            seen[idx] = true;
        }
    }
    seen.iter().all(|&b| b) // every subset covered
}

/// Corollary: a maximum matching between B_{n,k} and B_{n,k+1} (each edge a
/// proper-subset pair), read directly off the chains that pass through
/// both levels.
fn level_matching(chains: &[Vec<u32>], k: usize) -> Vec<(u32, u32)> {
    let mut edges = Vec::new();
    for chain in chains {
        for w in chain.windows(2) {
            if popcount(w[0]) as usize == k {
                edges.push((w[0], w[1]));
            }
        }
    }
    edges
}

// ---------------------------------------------------------------------
// Proof 1: LYM inequality (exact rational arithmetic)
// ---------------------------------------------------------------------

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

#[derive(Clone, Copy)]
struct Frac {
    num: i128,
    den: i128,
}
impl Frac {
    fn new(num: i128, den: i128) -> Self {
        let g = gcd(num, den).max(1);
        Frac { num: num / g, den: den / g }
    }
    fn add(self, other: Frac) -> Frac {
        Frac::new(self.num * other.den + other.num * self.den, self.den * other.den)
    }
    fn le_one(self) -> bool {
        self.num <= self.den
    }
}

/// LYM sum: for an antichain A, sum_{S in A} 1 / C(n, |S|).
fn lym_sum(n: usize, antichain: &[u32]) -> Frac {
    let mut acc = Frac::new(0, 1);
    for &s in antichain {
        let k = popcount(s) as usize;
        acc = acc.add(Frac::new(1, binomial(n, k) as i128));
    }
    acc
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
}

/// A random MAXIMAL (not necessarily maximum) antichain, built greedily:
/// shuffle all subsets, add each one that's incomparable to everything
/// already chosen.
fn random_maximal_antichain(n: usize, rng: &mut Rng) -> Vec<u32> {
    let total = 1u32 << n;
    let mut order: Vec<u32> = (0..total).collect();
    for i in (1..order.len()).rev() {
        let j = (rng.next() as usize) % (i + 1);
        order.swap(i, j);
    }
    let mut antichain: Vec<u32> = Vec::new();
    for s in order {
        let comparable = antichain.iter().any(|&t| (s & t == s) || (s & t == t));
        if !comparable {
            antichain.push(s);
        }
    }
    antichain
}

// ---------------------------------------------------------------------
// Independent cross-check: Dilworth's theorem via bipartite matching.
// max antichain size = min chain cover = (#elements) - (max matching in
// the "u properly-contained-in v" bipartite graph).
// ---------------------------------------------------------------------

fn max_bipartite_matching(adj: &[Vec<usize>], n_right: usize) -> usize {
    let n_left = adj.len();
    let mut match_right = vec![usize::MAX; n_right];
    let mut result = 0;
    for u in 0..n_left {
        let mut visited = vec![false; n_right];
        if try_kuhn(u, adj, &mut visited, &mut match_right) {
            result += 1;
        }
    }
    result
}

fn try_kuhn(u: usize, adj: &[Vec<usize>], visited: &mut [bool], match_right: &mut [usize]) -> bool {
    for &v in &adj[u] {
        if visited[v] {
            continue;
        }
        visited[v] = true;
        if match_right[v] == usize::MAX || try_kuhn(match_right[v], adj, visited, match_right) {
            match_right[v] = u;
            return true;
        }
    }
    false
}

fn max_antichain_via_dilworth(n: usize) -> usize {
    let total = 1usize << n;
    let mut adj = vec![Vec::new(); total];
    for u in 0..total {
        for v in 0..total {
            if u != v && (u & v) == u {
                adj[u].push(v); // u properly (or not) contained in v; equality excluded below
            }
        }
    }
    // remove u==v duplicates (u&v==u includes u==v; drop those)
    for u in 0..total {
        adj[u].retain(|&v| v != u);
    }
    let matching = max_bipartite_matching(&adj, total);
    total - matching
}

fn main() {
    // --- Symmetric chain decomposition: validity + count -----------------
    let mut decomp_ok = true;
    for n in 0..=10usize {
        let chains = symmetric_chain_decomposition(n);
        let valid = verify_symmetric_chain_decomposition(n, &chains);
        let expected = binomial(n, n / 2) as usize;
        if !valid || chains.len() != expected {
            decomp_ok = false;
            println!("DECOMP FAIL n={n}: valid={valid} chains={} expected={expected}", chains.len());
        }
    }
    println!("symmetric chain decomposition valid, with C(n, floor(n/2)) chains, for n=0..=10: {}", if decomp_ok { "PASS" } else { "FAIL" });

    // --- Level matching corollary -----------------------------------------
    let mut matching_ok = true;
    for n in 1..=8usize {
        let chains = symmetric_chain_decomposition(n);
        for k in 0..n {
            let edges = level_matching(&chains, k);
            let expected = binomial(n, k).min(binomial(n, k + 1)) as usize;
            let mut ss: Vec<u32> = edges.iter().map(|&(s, _)| s).collect();
            let mut ts: Vec<u32> = edges.iter().map(|&(_, t)| t).collect();
            ss.sort();
            ts.sort();
            ss.dedup();
            ts.dedup();
            let injective = ss.len() == edges.len() && ts.len() == edges.len();
            let proper = edges.iter().all(|&(s, t)| s & t == s && s != t);
            if edges.len() != expected || !injective || !proper {
                matching_ok = false;
                println!("MATCHING FAIL n={n} k={k}: got {} expected {expected} injective={injective} proper={proper}", edges.len());
            }
        }
    }
    println!("B_n,k <-> B_n,k+1 matching has size min(C(n,k),C(n,k+1)) and is injective: {}", if matching_ok { "PASS" } else { "FAIL" });

    // --- LYM inequality on random maximal antichains ----------------------
    let mut rng = Rng(7);
    let mut lym_ok = true;
    for _ in 0..300 {
        let n = 1 + (rng.next() % 8) as usize;
        let antichain = random_maximal_antichain(n, &mut rng);
        let sum = lym_sum(n, &antichain);
        if !sum.le_one() {
            lym_ok = false;
            println!("LYM FAIL n={n}: sum={}/{} > 1", sum.num, sum.den);
        }
    }
    println!("LYM inequality (sum of 1/C(n,|S|) <= 1) over random maximal antichains: {}", if lym_ok { "PASS" } else { "FAIL" });

    // Equality case: the full middle layer.
    let n = 6;
    let middle: Vec<u32> = (0..(1u32 << n)).filter(|&s| popcount(s) as usize == n / 2).collect();
    let sum = lym_sum(n, &middle);
    println!("LYM sum for the full middle layer B_{{6,3}}: {}/{} (should be exactly 1)", sum.num, sum.den);

    // --- Independent cross-check via Dilworth's theorem --------------------
    let mut dilworth_ok = true;
    for n in 0..=5usize {
        let via_dilworth = max_antichain_via_dilworth(n);
        let via_formula = binomial(n, n / 2) as usize;
        if via_dilworth != via_formula {
            dilworth_ok = false;
            println!("DILWORTH FAIL n={n}: matching-based={via_dilworth} formula={via_formula}");
        }
    }
    println!("max antichain size via Dilworth's theorem (bipartite matching) matches C(n,floor(n/2)) for n=0..=5: {}", if dilworth_ok { "PASS" } else { "FAIL" });
}
