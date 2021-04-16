use num::bigint::Sign;
use num::integer::*;
use num::iter::*;
use itertools::Itertools;

/// Modulo multiplication; `x*y mod m`.
pub fn mul_mod(x: u64, y: u64, m: u64) -> u64 {
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
/// Implements the algorithm from [Wikipedia: Jacobi symbol](https://en.wikipedia.org/wiki/Jacobi_symbol).
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
/// Refer to [Wikipedia: Strong pseudoprime](https://en.wikipedia.org/wiki/Strong_pseudoprime) for details.
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
/// Refer to [Wikipedia: Lucas pseudoprime](https://en.wikipedia.org/wiki/Lucas_pseudoprime#Strong_Lucas_pseudoprimes) for details.
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
/// Refer to [Wikipedia: Baillie-PSW primality test](https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test) for details.
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
/// Uses Pollard's rho algorithm internally. Refer to [Wikipedia: Pollard's rho algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm) for details.
/// "The rare special case" of failing with the default setting (initial value 2, increment by 1) is handled by trying
/// some more hardcoded settings of (2, -1), (3, 1), and (3, 2). Currently **panics** if all settings fail.
/// This will be patched when a failing input is found.
pub fn factorize(n: u64) -> Vec<u64> {
    let twos = n.trailing_zeros() as usize;
    let mut factor_2s = vec![2; twos];
    let mut odd_fac = factorize_odd(n >> twos);
    factor_2s.append(&mut odd_fac);
    factor_2s
}

/// Divisor sigma function, which includes the count and sum of divisors.
/// 
/// Refer to [Wikipedia: Divisor function](https://en.wikipedia.org/wiki/Divisor_function) for details.
/// `divisor_sigma(n, k)` is the sum of k-th power of all divisors of n.
/// The count of divisors of n can be found by giving `k = 0`, and the sum by giving `k = 1`.
/// 
/// Implemented using the prime factorization function and the formula given on the Wikipedia page.
pub fn divisor_sigma(n: u64, k: u32) -> u64 {
    let fac = factorize(n);
    let p_e_form = fac.iter().dedup_with_count();
    if k == 0 {
        p_e_form.map(|(e, _)| e as u64 + 1).product()
    } else {
        p_e_form.map(|(e, p)| (p.pow(k * (e as u32 + 1)) - 1) / (p.pow(k) - 1)).product()
    }
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

use num::integer::div_ceil;

/// An O(poly(log m)) algorithm to find the smallest x that satisfies `left <= a*x <= right (mod m)`.
pub fn mod_inv_range(a: u64, left: u64, right: u64, modulo: u64) -> u64 {
    if left == 0 { return 0; }
    let (a, left, right) = if 2 * a > modulo { (modulo-a, modulo-right, modulo-left) } else { (a, left, right) };
    let cc_1 = div_ceil(left, a);
    if a * cc_1 <= right { return cc_1; }
    let cc_2 = mod_inv_range(a - modulo % a, left % a, right % a, a);
    div_ceil(left as u128 + modulo as u128 * cc_2 as u128, a as u128) as u64
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

/// **Deprecated.** Please use PythagoreanTripletGenerator2 instead.
/// 
/// An iterator that infinitely generates all primitive Pythagorean triplets.
/// 
/// See [Wikipedia: Tree of primitive Pythagorean triples](https://en.wikipedia.org/wiki/Tree_of_primitive_Pythagorean_triples) for details.
/// This can also be used to generate any Pythagorean-like triplets in the form of `a^2 + b^2 = c^2 + k`,
/// by giving a suitable initial triple.
/// This iterator internally uses a BinaryHeap of a 3-tuple (a, b, c) with custom sorting order,
/// namely the inverted order of (c, a, b).
/// Therefore, in order to get a useful result, it is recommended to cut the iterator
/// with a condition on c (e.g. using take_while), and then do more computation on it.
/// Also note that the heap will contain 2n+1 items after n items are produced.
#[deprecated(
    since = "0.1.1",
    note = "Please use PythagoreanTripletGenerator2 instead"
)]
pub struct PythagoreanTripletGenerator {
    heap: BinaryHeap<PythagoreanTriplet>,
}

#[allow(deprecated)]
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

use std::collections::VecDeque;

/// An iterator that generates all primitive Pythagorean triplets which satisfy a given condition.
/// 
/// See [Wikipedia: Tree of primitive Pythagorean triples](https://en.wikipedia.org/wiki/Tree_of_primitive_Pythagorean_triples) for details.
/// This can also be used to generate any Pythagorean-like triplets in the form of `a^2 + b^2 = c^2 + k`,
/// by giving a suitable initial triple.
/// This iterator internally uses a VecDeque of a 3-tuple (a, b, c).
/// Note that the deque will contain at most 2n+1 items after n items are produced.
pub struct PythagoreanTripletGenerator2<F>
    where F: Fn(u64, u64, u64) -> bool
{
    deque: VecDeque<(u64, u64, u64)>,
    take_condition: F,
}
impl<F> Iterator for PythagoreanTripletGenerator2<F>
    where F: Fn(u64, u64, u64) -> bool
{
    type Item = (u64, u64, u64);
    fn next(&mut self) -> Option<(u64, u64, u64)> {
        if let Some(ret) = self.deque.pop_front() {
            let (a, b, c) = ret;
            let (a2, b2, c2) = (a + 2 * c - 2 * b, 2 * a + 2 * c - b, 2 * a + 3 * c - 2 * b);
            if (self.take_condition)(a2, b2, c2) {
                self.deque.push_back((a2, b2, c2));
            }
            let (a2, b2, c2) = (a + 2 * c + 2 * b, 2 * a + 2 * c + b, 2 * a + 3 * c + 2 * b);
            if (self.take_condition)(a2, b2, c2) {
                self.deque.push_back((a2, b2, c2));
            }
            let (a2, b2, c2) = (2 * b + 2 * c - a, b + 2 * c - 2 * a, 2 * b + 3 * c - 2 * a);
            if (self.take_condition)(a2, b2, c2) {
                self.deque.push_back((a2, b2, c2));
            }
            Some((a, b, c))
        }
        else { None }
    }
}

/// Creates a Pythagorean-like triplet generator with a custom starting triple.
/// 
/// It is guaranteed that all triples emitted share the value of a^2 + b^2 - c^2.
/// Also, if gcd(a, b, c) = 1, such a property is preserved across all triples.
#[deprecated(
    since = "0.1.1",
    note = "Please use pythagorean_like_triples2 instead"
)]
#[allow(deprecated)]
pub fn pythagorean_like_triples(a: u64, b: u64, c: u64) -> PythagoreanTripletGenerator {
    let mut heap: BinaryHeap<_> = BinaryHeap::new();
    heap.push(PythagoreanTriplet(a, b, c));
    PythagoreanTripletGenerator { heap }
}

/// Creates a Pythagorean triplet generator with the starting triple of (3, 4, 5).
#[deprecated(
    since = "0.1.1",
    note = "Please use pythagorean_triples2 instead"
)]
#[allow(deprecated)]
pub fn pythagorean_triples() -> PythagoreanTripletGenerator {
    pythagorean_like_triples(3, 4, 5)
}

/// Creates a Pythagorean-like triplet generator with a custom starting triple and a filtering function.
/// 
/// It is guaranteed that all triples emitted share the value of a^2 + b^2 - c^2.
/// Also, if gcd(a, b, c) = 1, such a property is preserved across all triples.
pub fn pythagorean_like_triples2<F>(a: u64, b: u64, c: u64, f: F) -> PythagoreanTripletGenerator2<F>
    where F: Fn(u64, u64, u64) -> bool
{
    let mut deque: VecDeque<_> = VecDeque::new();
    deque.push_back((a, b, c));
    PythagoreanTripletGenerator2 { deque, take_condition: f }
}

/// Creates a Pythagorean triplet generator with the starting triple of (3, 4, 5) and a filtering function.
pub fn pythagorean_triples2<F>(f: F) -> PythagoreanTripletGenerator2<F>
    where F: Fn(u64, u64, u64) -> bool
{
    pythagorean_like_triples2(3, 4, 5, f)
}

fn extended_gcd(a: u64, b: u64) -> (u64, u64, u64) {
    // returns (gcd, x, y) so that ax + by = gcd.
    let (mut old_r, mut r) = (a as i64, b as i64);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);
    while r != 0 {
        let q = old_r.div_euclid(r);
        let (temp_r, temp_s, temp_t) = (r, s, t);
        r = old_r - q * r;
        s = old_s - q * s;
        t = old_t - q * t;
        old_r = temp_r;
        old_s = temp_s;
        old_t = temp_t;
    }
    let gcd = old_r.abs();
    let a_coeff = old_s.rem_euclid(s.abs());
    let b_coeff = old_t.rem_euclid(t.abs());
    (gcd as u64, a_coeff as u64, b_coeff as u64)
}

/// Chinese remainder theorem for two moduli.
/// 
/// Returns the smallest integer x which satisfies `x == a1 (mod m1)` and `x == a2 (mod m2)`.
/// 
/// **Panics** if m1 and m2 are not coprime, even if the system is solvable.
pub fn chinese_remainder(a1: u64, m1: u64, a2: u64, m2: u64) -> u64 {
    let (gcd, n1, n2) = extended_gcd(m1, m2);
    if gcd != 1 {
        panic!("m1 and m2 are not coprime");
    }
    (a1 * n2 % m1 * m2 + a2 * n1 % m2 * m1) % (m1 * m2)
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
