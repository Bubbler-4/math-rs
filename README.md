# Bubblemath

![Crates.io](https://img.shields.io/crates/v/bubblemath)

A collection of mathematical algorithms in pure Rust. Documentation is available on [Docs.rs](https://docs.rs/bubblemath/0.1.0/bubblemath/).

## Notes

* I'm collecting the algorithms used while solving Project Euler problems.
* I won't publish my PE solutions directly, although it's being worked on side-by-side with the library.
* Please open an issue (or a pull request if applicable) if:
    * you know of a better implementation of an algorithm included in my library,
    * you want to optimize my library code further,
    * or you have any other suggestions :)

## Contents

* Number theory
    * Primality test (Baillie-PSW)
    * Integer factorization (Pollard's rho)
    * Next prime, n-th prime
    * Jacobi symbol
    * Iterator over Pythagorean triplets
    * Chinese remainder theorem for two moduli
* Linear recurrence
    * n-th term calculation (faster than matrix exponentiation by squaring)
    * Berlekamp-Massey algorithm (identifying linear recurrence modulo prime)

## Misc

![](https://projecteuler.net/profile/Bubbler.png)

Rust solution progress:

* 1-20 (revisiting early problems)
* 224, 258, 654