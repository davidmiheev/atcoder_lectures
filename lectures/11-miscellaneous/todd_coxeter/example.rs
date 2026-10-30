// The Todd-Coxeter algorithm for monoids. Given a finite alphabet and a
// finite list of relations u = v between words, it builds the right Cayley
// graph of the monoid they define, provided that monoid is finite. The graph
// is a transition table T[v][x]: the element reached from element v by
// multiplying on the right with the generator x, with the identity numbered 0.
//
// The algorithm keeps a growing set of provisional vertices, creates a vertex
// whenever a needed edge is undefined, and merges (with a union-find
// structure) vertices that must stand for the same element. If the monoid is
// finite it always stops; if it is infinite it never does, so the code takes
// a cap on the number of vertices and reports when that cap is hit.

// PDFBEGIN
use std::collections::{HashMap, VecDeque};

const UNDEFINED: usize = usize::MAX;

type Word = Vec<usize>;
type Relation = (Word, Word);

struct Enumerator {
    k: usize,
    relations: Vec<Relation>,
    table: Vec<Vec<usize>>, // table[v][x] = target vertex or UNDEFINED
    parent: Vec<usize>,     // union-find forest over vertices
    done: Vec<bool>,        // class fully processed?
    created: usize,         // total vertices ever created
}

impl Enumerator {
    fn new(num_letters: usize, relations: Vec<Relation>) -> Self {
        Enumerator { k: num_letters, relations, table: vec![], parent: vec![], done: vec![], created: 0 }
    }

    // ---- union-find (the smaller index always becomes the representative) ----
    fn find(&mut self, mut v: usize) -> usize {
        while self.parent[v] != v {
            self.parent[v] = self.parent[self.parent[v]];
            v = self.parent[v];
        }
        v
    }

    fn add_vertex(&mut self) -> usize {
        let v = self.table.len();
        self.table.push(vec![UNDEFINED; self.k]);
        self.parent.push(v);
        self.done.push(false);
        self.created += 1;
        v
    }

    /// Follow the edge labelled x from v, creating a vertex if it is undefined.
    fn step(&mut self, v: usize, x: usize) -> usize {
        let v = self.find(v);
        if self.table[v][x] == UNDEFINED {
            let w = self.add_vertex();
            self.table[v][x] = w;
        }
        let target = self.table[v][x];
        self.find(target)
    }

    fn walk(&mut self, mut v: usize, word: &[usize]) -> usize {
        for &x in word {
            v = self.step(v, x);
        }
        self.find(v)
    }

    /// Declare vertices p and q to be the same element, and propagate.
    fn merge(&mut self, p: usize, q: usize) {
        let mut pending = VecDeque::new();
        pending.push_back((p, q));
        while let Some((p, q)) = pending.pop_front() {
            let (mut p, mut q) = (self.find(p), self.find(q));
            if p == q {
                continue;
            }
            if q < p {
                std::mem::swap(&mut p, &mut q); // keep the smaller index as representative
            }
            self.parent[q] = p;
            self.done[p] = self.done[p] || self.done[q];
            for x in 0..self.k {
                let (a, b) = (self.table[p][x], self.table[q][x]);
                if a == UNDEFINED {
                    self.table[p][x] = b;
                } else if b != UNDEFINED {
                    pending.push_back((a, b)); // same element times x: targets coincide
                }
            }
        }
    }

    /// Returns the transition table (identity = 0, numbered in breadth-first
    /// order), or None if more than max_vertices vertices were needed.
    fn run(&mut self, max_vertices: usize) -> Option<Vec<Vec<usize>>> {
        self.add_vertex();
        let mut p = 0;
        while p < self.table.len() {
            if self.table.len() > max_vertices {
                return None;
            }
            let root = self.find(p);
            if !self.done[root] {
                for i in 0..self.relations.len() {
                    let (u, v) = self.relations[i].clone();
                    let (a, b) = (self.walk(p, &u), self.walk(p, &v));
                    self.merge(a, b);
                }
                for x in 0..self.k {
                    self.step(p, x);
                }
                let root = self.find(p);
                self.done[root] = true;
            }
            p += 1;
        }
        Some(self.canonical_table())
    }

    /// Renumber the surviving classes in breadth-first order from the identity.
    fn canonical_table(&mut self) -> Vec<Vec<usize>> {
        let root = self.find(0);
        let mut label: HashMap<usize, usize> = HashMap::new();
        label.insert(root, 0);
        let mut queue = VecDeque::new();
        queue.push_back(root);
        let mut order = vec![];
        while let Some(v) = queue.pop_front() {
            order.push(v);
            for x in 0..self.k {
                let target = self.table[v][x];
                let w = self.find(target);
                if !label.contains_key(&w) {
                    let next = label.len();
                    label.insert(w, next);
                    queue.push_back(w);
                }
            }
        }
        order
            .iter()
            .map(|&v| {
                (0..self.k)
                    .map(|x| {
                        let target = self.table[v][x];
                        label[&self.find(target)]
                    })
                    .collect()
            })
            .collect()
    }
}

fn todd_coxeter(num_letters: usize, relations: Vec<Relation>, max_vertices: usize) -> Option<Vec<Vec<usize>>> {
    Enumerator::new(num_letters, relations).run(max_vertices)
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn follow(table: &[Vec<usize>], word: &[usize], start: usize) -> usize {
    word.iter().fold(start, |v, &x| table[v][x])
}

fn parse(text: &str, letters: &str) -> Word {
    text.chars().filter(|&ch| ch != 'e').map(|ch| letters.find(ch).unwrap()).collect() // 'e' denotes the empty word
}

type Perm = Vec<usize>;

/// All permutations generated (as a monoid) by the given ones, by breadth-first search.
fn closure_of_permutations(generators: &[Perm]) -> Vec<Perm> {
    let size = generators[0].len();
    let identity: Perm = (0..size).collect();
    let mut seen: HashMap<Perm, usize> = HashMap::new();
    seen.insert(identity.clone(), 0);
    let mut all = vec![identity.clone()];
    let mut queue = VecDeque::new();
    queue.push_back(identity);
    while let Some(g) = queue.pop_front() {
        for h in generators {
            let gh: Perm = (0..size).map(|i| g[h[i]]).collect(); // (g h)(i) = g(h(i))
            if !seen.contains_key(&gh) {
                seen.insert(gh.clone(), seen.len());
                all.push(gh.clone());
                queue.push_back(gh);
            }
        }
    }
    all
}

fn check_group_presentation(name: &str, letters: &str, relations: &[(&str, &str)], generators: &[Perm], expected: usize) -> bool {
    let rels: Vec<Relation> = relations.iter().map(|(u, v)| (parse(u, letters), parse(v, letters))).collect();
    let mut enumerator = Enumerator::new(letters.len(), rels);
    let table = match enumerator.run(100_000) {
        Some(t) if t.len() == expected => t,
        other => {
            println!("  {name}: FAIL (got {:?}, expected {expected})", other.map(|t| t.len()));
            return false;
        }
    };
    // The map word -> permutation must be a bijection from the enumerated elements
    // onto the permutation group, which proves the presentation defines that group.
    let mut word_of: HashMap<usize, Word> = HashMap::new();
    word_of.insert(0, vec![]);
    let mut queue = VecDeque::new();
    queue.push_back(0usize);
    while let Some(v) = queue.pop_front() {
        for x in 0..letters.len() {
            let w = table[v][x];
            if !word_of.contains_key(&w) {
                let mut word = word_of[&v].clone();
                word.push(x);
                word_of.insert(w, word);
                queue.push_back(w);
            }
        }
    }
    let size = generators[0].len();
    let mut images: Vec<Perm> = word_of
        .values()
        .map(|word| {
            let mut perm: Perm = (0..size).collect();
            for &x in word {
                perm = (0..size).map(|i| perm[generators[x][i]]).collect();
            }
            perm
        })
        .collect();
    images.sort();
    images.dedup();
    let mut group = closure_of_permutations(generators);
    group.sort();
    let ok = images.len() == expected && images == group;
    println!("  {name}: {} elements, {} vertices created, isomorphic to the permutation group: {}", table.len(), enumerator.created, if ok { "PASS" } else { "FAIL" });
    ok
}

struct Lcg(u64);
impl Lcg {
    fn below(&mut self, bound: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) % bound as u64) as usize
    }
}

fn main() {
    let mut all_ok = true;

    // --- Worked example: a small monoid that is not a group ------------------
    let rels = vec![(parse("aa", "ab"), parse("a", "ab")), (parse("bb", "ab"), parse("b", "ab")), (parse("ba", "ab"), parse("ab", "ab"))];
    let mut enumerator = Enumerator::new(2, rels);
    let table = enumerator.run(1000).unwrap();
    println!("<a,b | aa=a, bb=b, ba=ab>: transition table (rows = elements, columns = a, b):");
    for (v, row) in table.iter().enumerate() {
        println!("  {v}: {:?}", row);
    }
    all_ok &= table == vec![vec![1, 2], vec![1, 3], vec![3, 2], vec![3, 3]];
    println!("  4 elements ({} vertices created along the way): {}", enumerator.created, if table.len() == 4 { "PASS" } else { "FAIL" });

    // --- Groups given by presentations, checked against permutation groups --------
    println!("\ngroups from presentations:");
    let s3: Vec<Perm> = vec![vec![1, 0, 2], vec![0, 2, 1]];
    all_ok &= check_group_presentation("S3 = <a,b | aaa=e, bb=e, ba=aab>", "ab", &[("aaa", "e"), ("bb", "e"), ("ba", "aab")], &[vec![1, 2, 0], vec![1, 0, 2]], 6);
    all_ok &= check_group_presentation("S3 = <a,b | aa=e, bb=e, aba=bab>", "ab", &[("aa", "e"), ("bb", "e"), ("aba", "bab")], &s3, 6);
    for n in 3..=8usize {
        let rot: Perm = (0..n).map(|i| (i + 1) % n).collect();
        let refl: Perm = (0..n).map(|i| (n - i) % n).collect();
        let a_n = "a".repeat(n);
        let a_n1 = format!("{}b", "a".repeat(n - 1));
        all_ok &= check_group_presentation(&format!("dihedral D{n}"), "ab", &[(&a_n, "e"), ("bb", "e"), ("ba", &a_n1)], &[rot, refl], 2 * n);
    }
    for (n, expected) in [(4usize, 24usize), (5, 120)] {
        let letters: String = "abcdefgh".chars().take(n - 1).collect();
        let chars: Vec<char> = letters.chars().collect();
        let mut rels_text: Vec<(String, String)> = chars.iter().map(|c| (format!("{c}{c}"), "e".to_string())).collect();
        for i in 0..n - 1 {
            for j in i + 1..n - 1 {
                let pair = format!("{}{}", chars[i], chars[j]);
                rels_text.push((pair.repeat(if j == i + 1 { 3 } else { 2 }), "e".to_string()));
            }
        }
        let gens: Vec<Perm> = (0..n - 1)
            .map(|i| {
                let mut perm: Perm = (0..n).collect();
                perm.swap(i, i + 1);
                perm
            })
            .collect();
        let rels_ref: Vec<(&str, &str)> = rels_text.iter().map(|(u, v)| (u.as_str(), v.as_str())).collect();
        all_ok &= check_group_presentation(&format!("symmetric group S{n} (Coxeter presentation)"), &letters, &rels_ref, &gens, expected);
    }

    // --- Infinite monoids never finish, so a cap is needed ---------------------
    let mut bicyclic = Enumerator::new(2, vec![(parse("ab", "ab"), parse("e", "ab"))]);
    let capped = bicyclic.run(300).is_none();
    println!("\n<a,b | ab=e> is infinite; the run is stopped by the vertex cap: {}", if capped { "PASS" } else { "FAIL" });
    all_ok &= capped;

    // --- Many provisional vertices, one final element ----------------------------
    for m in [10usize, 50, 200] {
        let mut e = Enumerator::new(1, vec![(vec![0; m], vec![]), (vec![0], vec![])]);
        let t = e.run(10_000).unwrap();
        println!("<a | a^{m}=e, a=e>: {} element, {} vertices created", t.len(), e.created);
        all_ok &= t == vec![vec![0]];
    }

    // --- Random presentations: soundness checks that need no known answer ----------
    let mut rng = Lcg(12345);
    let mut finite = 0;
    let tested = 200;
    let mut problems = 0;
    for _ in 0..tested {
        let k = 2 + rng.below(2);
        let mut rels: Vec<Relation> = vec![];
        for _ in 0..1 + rng.below(3) {
            let len_u = rng.below(4);
            let u: Word = (0..len_u).map(|_| rng.below(k)).collect();
            let len_v = rng.below(4);
            let v: Word = (0..len_v).map(|_| rng.below(k)).collect();
            rels.push((u, v));
        }
        let Some(table) = todd_coxeter(k, rels.clone(), 400) else { continue };
        finite += 1;
        let n = table.len();
        // (1) every relation holds at every element
        for p in 0..n {
            for (u, v) in &rels {
                if follow(&table, u, p) != follow(&table, v, p) {
                    problems += 1;
                }
            }
        }
        // (2) the multiplication read off the graph is associative
        let mut rep: HashMap<usize, Word> = HashMap::new();
        rep.insert(0, vec![]);
        let mut queue = VecDeque::new();
        queue.push_back(0usize);
        while let Some(v) = queue.pop_front() {
            for x in 0..k {
                let w = table[v][x];
                if !rep.contains_key(&w) {
                    let mut word = rep[&v].clone();
                    word.push(x);
                    rep.insert(w, word);
                    queue.push_back(w);
                }
            }
        }
        let mult: Vec<Vec<usize>> = (0..n).map(|g| (0..n).map(|h| follow(&table, &rep[&h], g)).collect()).collect();
        for g in 0..n {
            for h in 0..n {
                for i in 0..n {
                    if mult[mult[g][h]][i] != mult[g][mult[h][i]] {
                        problems += 1;
                    }
                }
            }
        }
        // (3) words connected by rewriting with the relations end at the same element
        for length in 0..6usize {
            for code in 0..k.pow(length as u32) {
                let word: Word = (0..length).map(|j| (code / k.pow(j as u32)) % k).collect();
                for (u, v) in &rels {
                    if u.len() > length {
                        continue;
                    }
                    for pos in 0..=length - u.len() {
                        if word[pos..pos + u.len()] == u[..] {
                            let mut other = word[..pos].to_vec();
                            other.extend_from_slice(v);
                            other.extend_from_slice(&word[pos + u.len()..]);
                            if follow(&table, &word, 0) != follow(&table, &other, 0) {
                                problems += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!("\n{finite} of {tested} random presentations define finite monoids (cap 400 vertices); relations hold, multiplication is associative, rewriting is respected: {}", if problems == 0 { "PASS" } else { "FAIL" });
    all_ok &= problems == 0;

    println!("\nAll checks: {}", if all_ok { "PASS" } else { "FAIL" });
}
