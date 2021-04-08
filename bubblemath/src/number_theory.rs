use num::bigint::Sign;
use num::integer::*;
use num::iter::*;
use itertools::Itertools;

/// Modulo multiplication; `x*y mod m`.
fn mul_mod(x: u64, y: u64, m: u64) -> u64 {
    (x as u128 * y as u128 % m as u128) as u64
}

/// Modulo power; `n^p mod m`.
pub fn pow_mod(n: u64, mut p: u64, m: u64) -> u64 {
    let mut ret = 1u64;
    let mut square = n;
    while p > 0 {
        if p % 2 != 0 { ret = mul_mod(ret, square, m); }
        square = mul_mod(square, square, m);
        p /= 2;
    }
    ret
}

/// Jacobi symbol `(a/n)`. The value is one of `-1, 0, 1`, which is represented as respective `Sign`.
/// 
/// Implements the algorithm from https://en.wikipedia.org/wiki/Jacobi_symbol.
pub fn jacobi_symbol(mut a: u64, mut n: u64) -> Sign {
    if gcd(a, n) != 1 { return Sign::NoSign; }
    let mut cur_sign = Sign::Plus;
    n = n >> n.trailing_zeros();
    a = a % n;
    loop {
        let twos_sign = if (n % 8 == 3 || n % 8 == 5) && a.trailing_zeros() % 2 != 0 { Sign::Minus } else { Sign::Plus };
        a = a >> a.trailing_zeros();
        cur_sign = cur_sign * twos_sign;
        if a == 1 { break; }
        let flip_sign = if a % 4 == 3 && n % 4 == 3 { Sign::Minus } else { Sign::Plus };
        let new_n = a;
        a = n % a;
        n = new_n;
        cur_sign = cur_sign * flip_sign;
    }
    cur_sign
}

/// Strong Fermat probable prime test with base 2 for odd `n >= 3`.
/// 
/// `n = d * 2^s + 1` (`d` is odd) is a strong Fermat probable prime if
/// `2^d == 1 mod n` or `2^(d * 2^r) == -1 mod n` for some `0 <= r < s`.
/// Refer to https://en.wikipedia.org/wiki/Strong_pseudoprime for details.
pub fn is_strong_probable_prime(n: u64) -> bool {
    let s = (n - 1).trailing_zeros();
    let d = n >> s;
    let mut lhs = pow_mod(2, d, n);
    if lhs == 1 { return true; }
    for _ in 0..s {
        if lhs == n - 1 { return true; }
        if lhs == 1 { return false; }
        lhs = mul_mod(lhs, lhs, n);
    }
    false
}

/// Strong Lucas probable prime test.
/// 
/// Refer to https://en.wikipedia.org/wiki/Lucas_pseudoprime#Strong_Lucas_pseudoprimes for details.
pub fn is_lucas_probable_prime(n: u64) -> bool {
    if sqrt(n).pow(2) == n { return false; }
    let d = range_step_from(5, 2)
        .find(|&x| jacobi_symbol(if x % 4 == 1 { x } else { n - x % n }, n) == Sign::Minus)
        .unwrap();
    let d = if d % 4 == 1 { d % n } else { n - d % n };
    let delta = n + 1;
    let small_s = delta.trailing_zeros();
    let small_d = delta >> small_s;
    let mut u = 1;
    let mut v = 1;
    let bit_width = 64 - small_d.leading_zeros();
    let halve_mod_n = |x| if x % 2 == 0 { x / 2 } else { (x + n) / 2 } % n;
    for shift in (0..bit_width - 1).rev() {
        let bit = (small_d >> shift) % 2;
        let mut u_next = mul_mod(u, v, n);
        let mut v_next = halve_mod_n(mul_mod(v, v, n) + mul_mod(mul_mod(d, u, n), u, n));
        if bit == 1 {
            let u_next2 = halve_mod_n(u_next + v_next);
            let v_next2 = halve_mod_n(mul_mod(d, u_next, n) + v_next);
            u_next = u_next2;
            v_next = v_next2;
        }
        u = u_next;
        v = v_next;
    }
    if u == 0 { return true; }
    for _ in 0..small_s {
        if v == 0 { return true; }
        let u_next = mul_mod(u, v, n);
        let v_next = halve_mod_n(mul_mod(v, v, n) + mul_mod(mul_mod(d, u, n), u, n));
        u = u_next;
        v = v_next;
    }
    false
}

/// Baillie-PSW primality test.
/// 
/// Refer to https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test for details.
/// It is known that no counterexamples exist in the range of u64.
pub fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    if n == 2 { return true; }
    if n % 2 == 0 { return false; }
    is_strong_probable_prime(n) && is_lucas_probable_prime(n)
}

fn find_factor_general(n: u64, init: u64, offset: u64) -> Option<u64> {
    let mut x = init;
    let mut y = init;
    let g = |x| (mul_mod(x, x, n) + offset) % n;
    loop {
        x = g(x);
        y = g(g(y));
        let gcd_xy = gcd(x + n - y, n);
        if gcd_xy != 1 {
            return if n == gcd_xy { None } else { Some(gcd_xy) };
        }
    }
}

fn find_factor(n: u64) -> u64 {
    [(2, 1), (2, n-1), (3, 1), (3, 2)].iter().filter_map(|&(init, offset)| find_factor_general(n, init, offset))
    .next().expect(&format!("Failed to find a non-trivial factor of {}!", n))
}

fn factorize_odd(n: u64) -> Vec<u64> {
    if n <= 1 { return vec![]; }
    if is_prime(n) { return vec![n]; }
    let fac = find_factor(n);
    let factors1 = factorize(fac);
    let factors2 = factorize(n / fac);
    factors1.iter().cloned().merge(factors2).collect()
}

/// Prime factorization.
/// 
/// Uses Pollard's rho algorithm internally. Refer to https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm for details.
/// "The rare special case" of failing with the default setting (initial value 2, increment by 1) is handled by trying
/// two more settings of (2, -1) and (3, 1).
pub fn factorize(n: u64) -> Vec<u64> {
    let twos = n.trailing_zeros() as usize;
    let mut factor_2s = vec![2; twos];
    let mut odd_fac = factorize_odd(n >> twos);
    factor_2s.append(&mut odd_fac);
    factor_2s
}

/// Tests if the given number is palindrome in the given base.
pub fn is_palindrome(n: u64, base: u64) -> bool {
    let mut cur_n = n;
    let mut mirror = 0;
    while cur_n > 0 {
        let (next_n, cur_digit) = div_rem(cur_n, base);
        mirror = mirror * base + cur_digit;
        cur_n = next_n;
    }
    mirror == n
}

/// Calculates the smallest prime greater than n.
pub fn next_prime(n: u64) -> u64 {
    (n+1..).find(|&x| is_prime(x)).unwrap()
}

/// Calculates the nth prime.
pub fn nth_prime(n: u64) -> u64 {
    let mut x = 1;
    for _ in 0..n { x = next_prime(x); }
    x
}

use std::collections::BinaryHeap;
use std::iter::Iterator;

#[derive(PartialEq, Eq, Default, Debug)]
struct PythagoreanTriplet(u64, u64, u64);
impl PartialOrd for PythagoreanTriplet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for PythagoreanTriplet {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // for use in max heap, compare c, a, then b in reverse
        (other.2, other.0, other.1).cmp(&(self.2, self.0, self.1))
    }
}

/// An iterator that infinitely generates all primitive Pythagorean triplets.
/// 
/// See https://en.wikipedia.org/wiki/Tree_of_primitive_Pythagorean_triples for details.
/// This can also be used to generate any Pythagorean-like triplets in the form of a^2 + b^2 = c^2 + k,
/// by giving a suitable initial triple.
/// This iterator internally uses a BinaryHeap of a 3-tuple (a, b, c) with custom sorting order,
/// namely the inverted order of (c, a, b).
/// Therefore, in order to get a useful result, it is recommended to cut the iterator
/// with a condition on c (e.g. using take_while), and then do more computation on it.
/// Also note that the heap will contain 2n+1 items after n items are produced.
pub struct PythagoreanTripletGenerator {
    heap: BinaryHeap<PythagoreanTriplet>,
}
impl Iterator for PythagoreanTripletGenerator {
    type Item = (u64, u64, u64);
    fn next(&mut self) -> Option<(u64, u64, u64)> {
        let ret = self.heap.pop().unwrap();
        let PythagoreanTriplet(a, b, c) = ret;
        self.heap.push(PythagoreanTriplet(a + 2 * c - 2 * b, 2 * a + 2 * c - b, 2 * a + 3 * c - 2 * b));
        self.heap.push(PythagoreanTriplet(a + 2 * c + 2 * b, 2 * a + 2 * c + b, 2 * a + 3 * c + 2 * b));
        self.heap.push(PythagoreanTriplet(2 * b + 2 * c - a, b + 2 * c - 2 * a, 2 * b + 3 * c - 2 * a));
        Some((a, b, c))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (usize::MAX, None)
    }
}

/// Creates a Pythagorean-like triplet generator with a custom starting triple.
/// 
/// It is guaranteed that all triples emitted share the value of a^2 + b^2 - c^2.
/// Also, if gcd(a, b, c) = 1, such a property is preserved across all triples.
pub fn pythagorean_like_triples(a: u64, b: u64, c: u64) -> PythagoreanTripletGenerator {
    let mut heap: BinaryHeap<_> = BinaryHeap::new();
    heap.push(PythagoreanTriplet(a, b, c));
    PythagoreanTripletGenerator { heap }
}

/// Creates a Pythagorean triplet generator with the starting triple of (3, 4, 5).
pub fn pythagorean_triples() -> PythagoreanTripletGenerator {
    pythagorean_like_triples(3, 4, 5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn primality_test() {
        let a001262 = range_step_from(3, 2).filter(|&x| is_strong_probable_prime(x) && !is_lucas_probable_prime(x)).take(10).collect_vec();
        assert_eq!(a001262, vec![2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633]);
        let a217255 = range_step_from(3, 2).filter(|&x| !is_strong_probable_prime(x) && is_lucas_probable_prime(x)).take(10).collect_vec();
        assert_eq!(a217255, vec![5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519]);
    }

    #[test]
    fn factorization() {
        // corner cases test
        assert_eq!(factorize(3486784401), vec![3; 20]);
        assert_eq!(factorize(95367431640625), vec![5; 20]);
        assert_eq!(factorize(5371), vec![41, 131]);
        for x in (2..2000).filter(|&x| is_prime(x)) {
            assert_eq!(factorize(x * x), vec![x, x]);
            assert_eq!(factorize(x * x * x), vec![x, x, x]);
        }
        // speed test
        assert_eq!(factorize(600851475143), vec![71, 839, 1471, 6857]);
        assert_eq!(factorize(123098432789), vec![132137, 931597]);
        assert_eq!(factorize(122233569105639659), vec![123789241, 987432899]);
    }
}
