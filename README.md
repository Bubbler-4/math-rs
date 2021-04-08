# Bubblemath

A collection of mathematical algorithms in Rust

## Notes

* I'm collecting the algorithms used while solving Project Euler problems.
* I won't publish my PE solutions directly, although it's being worked on side-by-side with the library.
* I will publish this library to crates.io when I feel it's ready.
* Please open an issue (or a pull request if applicable) if:
    * you know of a better implementation of an algorithm included in my library,
    * you want to optimize my library code further,
    * you want the library published ASAP, or
    * you have any other suggestions :)

## Contents

* Number theory
    * Primality test (Baillie-PSW)
    * Integer factorization (Pollard's rho)
    * Next prime, n-th prime
    * Jacobi symbol
    * Iterator over Pythagorean triplets

## Misc

![](https://projecteuler.net/profile/Bubbler.png)

Rust solution progress:

* Level `5%`
    * 1-11 (revisiting early problems)
* Level `45%`
    * 654 (used Berlekamp-Massey and n-th term acceleration, will include in the library later)
* Level `75%`
    * 224 (solvable without any specific mathy algorithm)