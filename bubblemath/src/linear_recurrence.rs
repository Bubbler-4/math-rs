use crate::number_theory::pow_mod;
use num::BigUint;

/// Berlekamp-Massey algorithm to find the minimal recurrence relation of a given sequence under a field.
/// 
/// Assumes that `modulo` is a prime under `2^31`, and every value in the given sequence is less than `modulo`.
/// The implementation is not very clean because it is almost directly ported from
/// [koosaga's C++ implementation for competitive programming](https://gist.github.com/koosaga/d4afc4434dbaa348d5bef0d60ac36aa4).
pub fn berlekamp_massey(x: &[u64], modulo: u64) -> Vec<u64> {
    let mut ls: Vec<u64> = vec![];
    let mut cur: Vec<u64> = vec![];
    let mut lf = 0usize;
    let mut ld = 0u64;
    for i in 0..x.len() {
        let mut t = 0u64;
        for j in 0..cur.len() {
            t = (t + x[i-j-1] * cur[j]) % modulo;
        }
        if t == x[i] { continue; }
        if cur.is_empty() {
            cur.resize(i+1, 0);
            lf = i;
            ld = (t + modulo - x[i]) % modulo;
            continue;
        }
        let k = (t + modulo - x[i]) % modulo * pow_mod(ld, modulo - 2, modulo) % modulo;
        let mut c = vec![0; i - lf - 1];
        c.push(k);
        for &j in &ls { c.push((modulo - j) * k % modulo); }
        if c.len() < cur.len() { c.resize(cur.len(), 0); }
        for j in 0..cur.len() {
            c[j] = (c[j] + cur[j]) % modulo;
        }
        if i - lf + ls.len() >= cur.len() {
            ls = cur; lf = i; ld = (t + modulo - x[i]) % modulo;
        }
        cur = c;
    }
    for i in &mut cur { *i %= modulo; }
    cur.reverse();
    cur
}

/// Faster-than-naive polynomial multiplication under a field.
/// 
/// Delegates the main algorithm to bigint multiplication implemented inside the `num` library.
/// The expected complexity (between `O(n^1.46)` and `O(n^1.58)`) is worse than the optimal
/// `O(n log n)` by FFT or related algorithms, but this works well enough in practice.
/// (NTT quickly becomes impractical when the modulus and the polynomial degree grow enough so that
/// `modulus^2 * degree > 2^64`.)
pub fn poly_mul(p1: &[u64], p2: &[u64], modulo: u64) -> Vec<u64> {
    let mut p1_adjusted: Vec<u32> = Vec::with_capacity(p1.len() * 3);
    for &p1_i in p1 {
        p1_adjusted.push(p1_i as u32);
        p1_adjusted.push(0);
        p1_adjusted.push(0);
    }
    let mut p2_adjusted: Vec<u32> = Vec::with_capacity(p2.len() * 3);
    for &p2_i in p2 {
        p2_adjusted.push(p2_i as u32);
        p2_adjusted.push(0);
        p2_adjusted.push(0);
    }
    let p1_biguint = BigUint::new(p1_adjusted);
    let p2_biguint = BigUint::new(p2_adjusted);
    let res = (&p1_biguint * &p2_biguint).to_u32_digits();
    let mut digits: Vec<u64> = res.chunks(3).map(|chunk| {
        match chunk.len() {
            1 => chunk[0] as u64 % modulo,
            2 => (chunk[0] as u64 + chunk[1] as u64 * 2u64.pow(32)) % modulo,
            3 => (chunk[0] as u64 + (chunk[1] as u64 + chunk[2] as u64 * 2u64.pow(32)) % modulo * 2u64.pow(32)) % modulo,
            _ => unreachable!()
        }
    }).collect();
    digits.resize(p1.len() + p2.len() - 1, 0);
    digits
}

fn one_coeff(p: &mut [u64], q: &mut [u64], mut n: u64, modulo: u64) -> u64 {
    while n >= 1 {
        let d = p.len();
        let q_minus: Vec<u64> = q.iter().enumerate().map(|(i, &q_i)| if i % 2 == 0 { q_i } else { (modulo - q_i) % modulo }).collect();
        let u = poly_mul(p, &q_minus, modulo);
        let n_bit = (n % 2) as usize;
        for i in 0..d {
            p[i] = u[2 * i + n_bit];
        }
        let a = poly_mul(q, &q_minus, modulo);
        for i in 0..=d {
            q[i] = a[2 * i];
        }
        n /= 2;
    }
    p[0] * pow_mod(q[0], modulo - 2, modulo) % modulo
}

/// Fast algorithm to get the n-th term of a linear recurrence relation under a field.
/// 
/// Assumes `modulo` is a prime under `2^31`, and every input value in `recurrence` and `initial` is less than `modulo`.
/// Uses the algorithm given in [the paper published in Jan 2021](https://hal.inria.fr/hal-02917827v2/document).
pub fn nth_term(recurrence: &[u64], initial: &[u64], n: u64, modulo: u64) -> u64 {
    let d = recurrence.len();
    let mut q: Vec<u64> = Vec::with_capacity(d+1);
    q.push(1);
    for &ci in recurrence.iter().rev() { q.push((modulo - ci) % modulo); }
    let mut p = poly_mul(initial, &q, modulo);
    p.resize(d, 0);
    one_coeff(&mut p, &mut q, n, modulo)
}