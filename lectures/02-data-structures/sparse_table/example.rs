// Sparse table for idempotent range queries (here: range minimum).
//
// Build: O(n log n). Query: O(1).

struct SparseTable {
    // table[k][i] = min of a[i .. i + 2^k)
    table: Vec<Vec<i64>>,
    log2floor: Vec<u32>,
}

impl SparseTable {
    fn new(a: &[i64]) -> Self {
        let n = a.len();
        let mut log2floor = vec![0u32; n + 1];
        for i in 2..=n {
            log2floor[i] = log2floor[i / 2] + 1;
        }
        let max_k = if n == 0 { 1 } else { log2floor[n] as usize + 1 };

        let mut table = vec![a.to_vec()];
        for k in 1..max_k {
            let half = 1usize << (k - 1);
            let len = 1usize << k;
            if len > n {
                break;
            }
            let prev = &table[k - 1];
            let mut row = Vec::with_capacity(n - len + 1);
            for i in 0..=n - len {
                row.push(prev[i].min(prev[i + half]));
            }
            table.push(row);
        }
        SparseTable { table, log2floor }
    }

    /// Minimum of a[l..r), 0 <= l < r <= n.
    fn query(&self, l: usize, r: usize) -> i64 {
        let k = self.log2floor[r - l] as usize;
        let window = 1usize << k;
        self.table[k][l].min(self.table[k][r - window])
    }
}

fn main() {
    let a: Vec<i64> = vec![5, 2, 4, 7, 1, 3, 6];
    let st = SparseTable::new(&a);

    let queries = [(1usize, 6usize), (0, 7), (3, 4), (4, 6)];
    for (l, r) in queries {
        println!("min(a[{}..{})) = {}", l, r, st.query(l, r));
    }
}
