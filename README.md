# Prime Walk: Deterministic Four-Square Prime Quaternions for SO(3)

A C++17 computational framework for constructing four-square decompositions of primes and exploring normalized quaternion geometry on $S^3$ and $SO(3)$.

## Overview

Prime Walk combines parity-guided coordinate search with integer four-square decomposition:

- Uses parity constraints to ensure that the remaining two-square term $R \equiv 1 \pmod 4$.
- Tests $R$ for primality and decomposes it with Tonelli–Shanks and Cornacchia's method.
- Normalizes the resulting $(c,d,e,f)$ tuple to a unit quaternion.
- Includes experiments for geometry, Riesz energy, and SO(3) quadrature.

The repository contains two implementation eras: an original 1-billion-prime cyclic test and later SO(3) / Haar-mapping versions.

## Implementations

| Path | Purpose |
| :--- | :--- |
| [`main.cpp`](main.cpp) | Main Haar-mapped SO(3) experiment for generating positive-octant unit quaternions and comparing prime-walk, Sobol, Halton, and Monte Carlo quadrature. |
| [`V0/main.cpp`](V0/main.cpp) | Original segmented-sieve and cyclic four-square test for primes up to $10^9$. |
| [`V1/A/main.cpp`](V1/A/main.cpp) | Earlier Riesz $s=2$ energy implementation using a Haar-mapped $S^3_+$ prime walk and $S_4$ permutations. |

The remaining implementation and report files are listed in [Repository Contents](#repository-contents).

## Build and Run

Requires a C++17 compiler. OpenMP is used by the supplied implementations.

### Main Haar / SO(3) benchmark

```bash
g++ -O3 -march=native -fopenmp -std=c++17 main.cpp -o prime_walk
./prime_walk
```

The executable streams primes above $10^7$, constructs four-square decompositions, and reports SO(3) integration errors and convergence exponents.

### 1-billion-prime cyclic test

```bash
g++ -O3 -march=native -fopenmp -std=c++17 V0/main.cpp -o prime_walk_v0
./prime_walk_v0 --limit 1000000000 --sample 1000
```

For a full check of every prime below 100 million:

```bash
./prime_walk_v0 --limit 100000000 --sample 1
```

> Under the default `--sample 1000`, only every sampled prime is decomposition-checked; streamed primes outside the sample are not individually four-square checked.

### V1 Riesz-energy benchmark

```bash
g++ -O3 -march=native -fopenmp -std=c++17 V1/A/main.cpp -o prime_walk_v1a
./prime_walk_v1a
```

## Repository Contents

- [`main.cpp`](main.cpp) — Haar-mapped $S^3_+$ generator and SO(3) quadrature benchmarks.
- [`V0/`](V0/) — original large-scale prime walk test and output/report artifacts.
- [`V1/A/`](V1/A/) — $S_4$-orbit Haar prime-walk and Riesz-energy benchmark.
- [`V1/B/`](V1/B/) — additional V1 source and results.
- [`output.txt`](output.txt) — representative output for the main implementation.
- [`Prime_Walk_SO3_Quadrature_IEEE_Paper.tex`](Prime_Walk_SO3_Quadrature_IEEE_Paper.tex) — LaTeX manuscript source.
- [`Prime_Walk_SO3_Quadrature_IEEE_Paper.pdf`](Prime_Walk_SO3_Quadrature_IEEE_Paper.pdf) — compiled IEEE-style manuscript.
- [`Prime_Walk_QMC_and_SO3_Quadrature_Report.pdf`](Prime_Walk_QMC_and_SO3_Quadrature_Report.pdf) — supplemental report.

## Method Summary

For an odd prime $p$, the implementation searches for positive integers $(c,d)$ with parity determined by $p \bmod 4$, leaving:

$$R = p-c^2-d^2 \equiv 1 \pmod 4.$$

If $R$ is prime, Cornacchia's method decomposes $R=e^2+f^2$, yielding:

$$p=c^2+d^2+e^2+f^2.$$

Dividing the tuple by $\sqrt{p}$ gives a unit quaternion $q=(c,d,e,f)/\sqrt{p}$ that is used by the SO(3) experiment.

## Geometry Checks

The original implementation checks algebraic geometry derived from the four-square components, including normalized determinant expressions and cyclic-quadrilateral area/Ptolemy tests. These results are experimental unless proved and stated in the accompanying manuscript.

## Reproducing the LaTeX Paper

```bash
pdflatex Prime_Walk_SO3_Quadrature_IEEE_Paper.tex
pdflatex Prime_Walk_SO3_Quadrature_IEEE_Paper.tex
```

## License

This project is licensed under the [MIT License](LICENSE).

## Language Composition

- C++: primary implementation and benchmarks.
- TeX: manuscript and report source.
