//! A collection of mathematical algorithms written in pure Rust.
//! 
//! ## Contents
//! 
//! * Number theory
//!     * Primality test (Baillie-PSW)
//!     * Integer factorization (Pollard's rho)
//!     * Next prime, n-th prime
//!     * Jacobi symbol
//!     * Iterator over Pythagorean triplets
//!     * Chinese remainder theorem for two moduli
//! * Linear recurrence
//!     * n-th term calculation (faster than matrix exponentiation by squaring)
//!     * Berlekamp-Massey algorithm (identifying linear recurrence modulo prime)

pub mod number_theory;
pub mod linear_recurrence;