# Prime Walk: Deterministic Four-Square Prime Quaternions for $\mathrm{SO}(3)$ Quadrature

A high-performance C++17 computational framework for constructing balanced integer four-square prime decompositions ($p = c^2 + d^2 + e^2 + f^2$) and mapping them via exact $\mathrm{Beta}(1/2, 1/2)$ Haar geometry onto $S^3_+ / \mathrm{SO}(3)$ for seed-free Quasi-Monte Carlo (QMC) numerical integration.

---

## Overview

1. **Parity-Aligned Four-Square Solver (`prime_1billion.cpp`)**
   * Steers $(c, d)$ near target coordinates with strict parity enforcement ($c, d$ even for $p \equiv 1 \pmod 4$; odd for $p \equiv 3 \pmod 4$) so that $R = p - c^2 - d^2 \equiv 1 \pmod 4$.
   * When $R$ passes deterministic Miller–Rabin primality testing, solves $R = e^2 + f^2$ via Tonelli–Shanks modular square roots and Cornacchia's algorithm.
   * Verifies universal algebraic bounds ($\max |\mathrm{Det}| / p \le 0.5$) and cyclic quadrilateral area bounds ($\mathcal{A}/p \le 0.25$).

2. **Exact Haar-Measure Mapping & $\mathrm{SO}(3)$ Quadrature (`so3_prime_quadrature.cpp`)**
   * Maps the 2D Plastic-Constant ($R_2$) low-discrepancy sequence through the inverse $\mathrm{Beta}(1/2, 1/2)$ CDF ($F^{-1}(u) = \sin^2(\frac{\pi}{2}u)$) so normalized prime quaternions $(c,d,e,f)/\sqrt{p} \in S^3_+$ follow the exact $\mathrm{Dirichlet}(1/2, 1/2, 1/2, 1/2)$ Haar volume element.
   * Applies the 24-element $S_4$ coordinate permutation orbit and 8-element octant sign group to cover $\mathrm{SO}(3)$ without boundary or ordering bias.

---

## Key Empirical Results

### 1. Arithmetic & Cryptographic Prime Verification

| Test Scope / Family | Bit Size | Primes Tested | Success Rate | Max $\text{maxDet}$ |
| :--- | :--- | :--- | :--- | :--- |
| **Segmented Sieve Stream ($\le 10^9$)** | $\le 30$-bit | 50,847,534 (50,847 sampled) | **100%** | `0.500000` |
| **Trillion & Quadrillion Primes** | 40–50-bit | $10^{12}+39$, $10^{15}+37$ | **100%** ($k \le 38$) | `0.494812` |
| **Sophie Germain & Safe Primes** | 61–258-bit | 90 primes | **100%** ($k \le 449$) | `0.499757` |
| **Solinas, Curve & ZK-STARK Fields** | 31–521-bit | 20 primes (incl. M31, Curve25519, BLS12-381, P-521) | **100%** ($k \le 257$) | `0.499999` |
| **RFC 3526 Diffie–Hellman Safe Primes** | 1024–2048-bit | Group 2 & Group 14 | **100%** ($k \le 593$) | `0.427450` |

### 2. $S^3_+$ Uniformity & Unbiased $\mathrm{SO}(3)$ Integration ($N = 1024$, $p > 10^7$)
* **Acceptance Rate:** `2048 / 2048` (**100%**) within diamond search radius $\le 250$.
* **$S^3_+$ Riesz $s=2$ Energy:** `5.013241` (matching exact Haar Monte Carlo `5.023201`; 140$\times$ reduction vs. unsteered balanced slices `702.0227`).

| Integrand on $\mathrm{SO}(3)$ | Exact Mean | Sobol' Error | Halton Error | **Haar Prime Walk** | Monte Carlo RMSE |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $f_1(R) = \mathrm{Tr}(R)$ | `0.0` | `0.001873` | `0.004995` | **`0.006945`** | `0.034058` |
| $f_2(R) = \mathrm{Tr}(R^2)$ | `1.0` | `0.002271` | `0.004641` | **`0.032117`** | `0.048501` |
| $f_3(R) = \mathbf{v}^T R^T M R \mathbf{v}$ | `2.0` | `0.000349` | `0.000059` | **`0.003458`** | `0.009157` |

---

## Build & Run

Requires a C++17 compiler with OpenMP support (`g++` or `clang++`).

```bash
# 1. Compile and run the 1-Billion Prime Walk stream test
g++ -O3 -march=native -fopenmp -std=c++17 prime_1billion.cpp -o prime_1billion
./prime_1billion --limit 1000000000 --sample 1000

# 2. Compile and run the S^3_+ Riesz energy & SO(3) QMC quadrature benchmark
g++ -O3 -march=native -fopenmp -std=c++17 so3_prime_quadrature.cpp -o so3_prime_quadrature
./so3_prime_quadrature
```

---

## Repository Contents

* `prime_1billion.cpp` — Segmented prime sieve and OpenMP balanced four-square verifier up to $10^9$.
* `so3_prime_quadrature.cpp` — Unbiased $\mathrm{Beta}(1/2, 1/2)$ Haar Prime Walk generator and $\mathrm{SO}(3)$ integration benchmark against Sobol', Halton, and Monte Carlo.
* `Prime_Walk_SO3_Quadrature_IEEE_Paper.tex` / `.pdf` — Formal two-column IEEE manuscript detailing proofs, algorithms, and empirical benchmarks.
