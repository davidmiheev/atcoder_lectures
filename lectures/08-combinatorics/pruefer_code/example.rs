// Pruefer code: a bijection between labeled trees on {0,...,n-1} and
// sequences in {0,...,n-1}^(n-2).  Encoding repeatedly deletes the smallest
// leaf and records its neighbor; decoding rebuilds the tree from the fact
// that a vertex occurring c times in the code has degree c+1.  The harness
// below cross-checks the heap and linear versions and Cayley's formula.

// PDFBEGIN
use std::cmp::Reverse;
use std::collections::BinaryHeap;

type Edge = (usize, usize);

/// Pruefer code of a labeled tree, O(n log n) with a min-heap.
fn encode_heap(n: usize, edges: &[Edge]) -> Vec<usize> {
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut deg: Vec<usize> = adj.iter().map(|x| x.len()).collect();
    let mut leaves: BinaryHeap<Reverse<usize>> =
        (0..n).filter(|&v| deg[v] == 1).map(Reverse).collect();
    let mut code = Vec::with_capacity(n.saturating_sub(2));
    for _ in 0..n.saturating_sub(2) {
        let Reverse(u) = leaves.pop().unwrap(); // smallest current leaf
        let v = *adj[u].iter().find(|&&w| deg[w] > 0).unwrap(); // its live neighbor
        code.push(v);
        deg[u] = 0;
        deg[v] -= 1;
        if deg[v] == 1 {
            leaves.push(Reverse(v));
        }
    }
    code
}

/// Inverse of `encode_heap`, O(n log n).
fn decode_heap(n: usize, code: &[usize]) -> Vec<Edge> {
    let mut left = vec![0usize; n]; // occurrences in the unread part of the code
    for &v in code {
        left[v] += 1;
    }
    let mut free: BinaryHeap<Reverse<usize>> =
        (0..n).filter(|&v| left[v] == 0).map(Reverse).collect();
    let mut edges = Vec::with_capacity(n - 1);
    for &v in code {
        let Reverse(u) = free.pop().unwrap();
        edges.push((u, v));
        left[v] -= 1;
        if left[v] == 0 {
            free.push(Reverse(v));
        }
    }
    let Reverse(a) = free.pop().unwrap();
    let Reverse(b) = free.pop().unwrap();
    edges.push((a, b));
    edges
}

/// Pruefer code in O(n): root the tree at n-1 and sweep a pointer.
fn encode_linear(n: usize, edges: &[Edge]) -> Vec<usize> {
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut parent = vec![usize::MAX; n];
    let mut stack = vec![n - 1];
    while let Some(v) = stack.pop() {
        // orient every edge toward the root n-1
        for &w in &adj[v] {
            if w != parent[v] {
                parent[w] = v;
                stack.push(w);
            }
        }
    }
    let mut deg: Vec<usize> = adj.iter().map(|x| x.len()).collect();
    let mut ptr = 0;
    while deg[ptr] != 1 {
        ptr += 1;
    }
    let mut leaf = ptr; // smallest current leaf
    let mut code = Vec::with_capacity(n.saturating_sub(2));
    for _ in 0..n.saturating_sub(2) {
        let p = parent[leaf];
        code.push(p);
        deg[p] -= 1;
        if deg[p] == 1 && p < ptr {
            leaf = p; // p became a leaf smaller than everything unseen
        } else {
            ptr += 1;
            while deg[ptr] != 1 {
                ptr += 1;
            }
            leaf = ptr;
        }
    }
    code
}

/// Inverse of `encode_linear`, O(n).
fn decode_linear(n: usize, code: &[usize]) -> Vec<Edge> {
    let mut deg = vec![1usize; n];
    for &v in code {
        deg[v] += 1;
    }
    let mut ptr = 0;
    while deg[ptr] != 1 {
        ptr += 1;
    }
    let mut leaf = ptr;
    let mut edges = Vec::with_capacity(n - 1);
    for &v in code {
        edges.push((leaf, v));
        deg[v] -= 1;
        if deg[v] == 1 && v < ptr {
            leaf = v;
        } else {
            ptr += 1;
            while deg[ptr] != 1 {
                ptr += 1;
            }
            leaf = ptr;
        }
    }
    edges.push((leaf, n - 1));
    edges
}
// PDFEND -- everything below is the test harness and is not printed in the PDF

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, m: usize) -> usize {
        (self.next() % m as u64) as usize
    }
}

fn canon(edges: &[Edge]) -> Vec<Edge> {
    let mut v: Vec<Edge> = edges.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
    v.sort();
    v
}

fn find(p: &mut Vec<usize>, x: usize) -> usize {
    let mut r = x;
    while p[r] != r {
        r = p[r];
    }
    let mut c = x;
    while p[c] != r {
        let nx = p[c];
        p[c] = r;
        c = nx;
    }
    r
}

fn is_tree(n: usize, edges: &[Edge]) -> bool {
    if edges.len() != n - 1 {
        return false;
    }
    let mut p: Vec<usize> = (0..n).collect();
    for &(a, b) in edges {
        let (ra, rb) = (find(&mut p, a), find(&mut p, b));
        if ra == rb {
            return false;
        }
        p[ra] = rb;
    }
    true
}

/// Random labeled tree NOT produced by Pruefer decoding, so the round-trip
/// tests are not circular: attach each vertex to an earlier one, relabel.
fn random_tree(n: usize, rng: &mut Rng) -> Vec<Edge> {
    let mut perm: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        perm.swap(i, rng.below(i + 1));
    }
    (1..n).map(|i| (perm[i], perm[rng.below(i)])).collect()
}

/// Number of labeled trees on n vertices, by brute force over edge subsets.
fn count_trees_bruteforce(n: usize) -> usize {
    let mut pairs = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            pairs.push((a, b));
        }
    }
    let m = pairs.len();
    let mut count = 0;
    for mask in 0u32..(1 << m) {
        if mask.count_ones() as usize != n - 1 {
            continue;
        }
        let edges: Vec<Edge> = (0..m).filter(|&i| mask >> i & 1 == 1).map(|i| pairs[i]).collect();
        if is_tree(n, &edges) {
            count += 1;
        }
    }
    count
}

fn main() {
    // --- the lecture's worked examples -----------------------------------
    let ex = [(0, 4), (1, 4), (4, 7), (2, 7), (7, 3), (3, 6), (5, 6), (6, 8)];
    let code = encode_heap(9, &ex);
    assert_eq!(code, vec![4, 4, 7, 7, 6, 3, 6]);
    assert_eq!(code, encode_linear(9, &ex));
    assert_eq!(canon(&decode_heap(9, &code)), canon(&ex));
    let t2 = decode_heap(8, &[3, 3, 0, 5, 0, 3]);
    assert_eq!(canon(&t2), vec![(0, 3), (0, 4), (0, 5), (1, 3), (2, 3), (3, 7), (5, 6)]);
    println!("worked examples: PASS");

    // --- Cayley's formula by brute force, n = 2..=6 ----------------------
    for n in 2..=6usize {
        let c = count_trees_bruteforce(n);
        assert_eq!(c, n.pow(n as u32 - 2), "n={}", n);
    }
    println!("brute-force tree counts equal n^(n-2) for n<=6: PASS");

    // --- every sequence decodes to a tree whose code is the sequence -----
    for n in 2..=7usize {
        let len = n - 2;
        let total = n.pow(len as u32);
        for mut x in 0..total {
            let mut code = vec![0usize; len];
            for c in code.iter_mut() {
                *c = x % n;
                x /= n;
            }
            let t = decode_heap(n, &code);
            assert!(is_tree(n, &t));
            assert_eq!(encode_heap(n, &t), code);
            assert_eq!(canon(&t), canon(&decode_linear(n, &code)));
        }
    }
    println!("decode/encode round trip on all sequences, n<=7: PASS");

    // --- random larger trees: heap vs linear, degree lemma ---------------
    let mut rng = Rng(0x9E3779B97F4A7C15);
    for _ in 0..500 {
        let n = 2 + rng.below(300);
        let t = random_tree(n, &mut rng);
        let c = encode_heap(n, &t);
        assert_eq!(c, encode_linear(n, &t));
        assert_eq!(canon(&decode_linear(n, &c)), canon(&t));
        assert_eq!(canon(&decode_heap(n, &c)), canon(&t));
        let mut deg = vec![0usize; n];
        for &(a, b) in &t {
            deg[a] += 1;
            deg[b] += 1;
        }
        let mut cnt = vec![0usize; n];
        for &v in &c {
            cnt[v] += 1;
        }
        assert!((0..n).all(|v| deg[v] == 1 + cnt[v]));
    }
    println!("random trees: heap == linear, deg(v) = 1 + #occurrences: PASS");
    println!("ALL PASS");
}
