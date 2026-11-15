// Prime sieves: Eratosthenes, the linear sieve (with smallest prime factors),
// bulk multiplicative functions, and a segmented sieve for a high window.

// PDFBEGIN
/// Primes up to n, striking multiples of p starting at p*p.
fn sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    if n < 2 {
        return Vec::new();
    }
    let mut is_prime = vec![true; n + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    let mut p = 2usize;
    while p * p <= n {
        if is_prime[p] {
            let mut m = p * p;
            while m <= n {
                is_prime[m] = false;
                m += p;
            }
        }
        p += 1;
    }
    (2..=n).filter(|&i| is_prime[i]).collect()
}

/// Linear sieve. Every composite is written exactly once, by its smallest
/// prime factor. Returns (primes, spf) with spf[0] = spf[1] = 0.
fn linear_sieve(n: usize) -> (Vec<usize>, Vec<usize>) {
    let mut spf = vec![0usize; n + 1];
    let mut primes: Vec<usize> = Vec::new();

    for i in 2..=n {
        if spf[i] == 0 {
            spf[i] = i;
            primes.push(i);
        }
        for &p in &primes {
            // The guard p <= spf[i] is what makes each composite unique.
            if p > spf[i] || p * i > n {
                break;
            }
            spf[p * i] = p;
        }
    }
    (primes, spf)
}

/// Factor n into (prime, exponent) pairs via repeated spf lookups.
fn factorize(mut n: usize, spf: &[usize]) -> Vec<(usize, u32)> {
    let mut factors: Vec<(usize, u32)> = Vec::new();
    while n > 1 {
        let p = spf[n];
        let mut e = 0u32;
        while n % p == 0 {
            n /= p;
            e += 1;
        }
        factors.push((p, e));
    }
    factors
}

/// Euler's totient for 1..=n, using phi(p^e) = p^(e-1) (p - 1).
fn totient_table(n: usize, spf: &[usize]) -> Vec<u64> {
    let mut phi = vec![0u64; n + 1];
    if n >= 1 {
        phi[1] = 1;
    }
    for i in 2..=n {
        let p = spf[i];
        let (mut m, mut e) = (i, 0u32);
        while m % p == 0 {
            m /= p;
            e += 1;
        }
        let mut pe = (p - 1) as u64;
        for _ in 1..e {
            pe *= p as u64;
        }
        phi[i] = phi[m] * pe;
    }
    phi
}

/// Divisor count for 1..=n, using d(p^e) = e + 1.
fn divisor_count_table(n: usize, spf: &[usize]) -> Vec<u64> {
    let mut d = vec![0u64; n + 1];
    if n >= 1 {
        d[1] = 1;
    }
    for i in 2..=n {
        let p = spf[i];
        let (mut m, mut e) = (i, 0u64);
        while m % p == 0 {
            m /= p;
            e += 1;
        }
        d[i] = d[m] * (e + 1);
    }
    d
}

/// Primes in [lo, hi], for hi far beyond what a full sieve could hold.
fn segmented_sieve(lo: u64, hi: u64) -> Vec<u64> {
    if hi < 2 {
        return Vec::new();
    }
    let lo = lo.max(2);
    if lo > hi {
        return Vec::new();
    }

    let mut root = (hi as f64).sqrt() as u64 + 2;
    while root * root > hi {
        root -= 1;
    }
    let base = sieve_of_eratosthenes(root as usize);

    let len = (hi - lo + 1) as usize;
    let mut is_prime = vec![true; len];

    for &p in &base {
        let p = p as u64;
        // First multiple of p inside the window, but never below p*p.
        let mut start = ((lo + p - 1) / p) * p;
        if start < p * p {
            start = p * p;
        }
        let mut m = start;
        while m <= hi {
            is_prime[(m - lo) as usize] = false;
            m += p;
        }
    }

    (0..len)
        .filter(|&i| is_prime[i])
        .map(|i| lo + i as u64)
        .collect()
}

// PDFEND -- everything below is the test harness and is not printed in the PDF

fn main() {
    let n = 50usize;

    let primes = sieve_of_eratosthenes(n);
    println!("primes up to {}: {:?}", n, primes);

    let (lin_primes, spf) = linear_sieve(n);
    println!("linear sieve agrees: {}", lin_primes == primes);
    println!("spf[2..=12] = {:?}", &spf[2..=12]);

    for m in [360usize, 97, 1024] {
        let f = factorize(m, &linear_sieve(m).1);
        let pretty: Vec<String> = f.iter().map(|(p, e)| format!("{}^{}", p, e)).collect();
        println!("{} = {}", m, pretty.join(" * "));
    }

    let phi = totient_table(n, &spf);
    println!("phi(1..=12) = {:?}", &phi[1..=12]);

    let d = divisor_count_table(n, &spf);
    println!("d(1..=12)   = {:?}", &d[1..=12]);

    let (lo, hi) = (1_000_000_000_000u64, 1_000_000_000_100u64);
    let window = segmented_sieve(lo, hi);
    println!("primes in [{}, {}]: {:?}", lo, hi, window);
}
