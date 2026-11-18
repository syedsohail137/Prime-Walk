// Block 1: Headers & Fast 32/64-bit Modular Number Theory
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

using u64 = uint64_t;
using u32 = uint32_t;
using i64 = int64_t;

static inline u64 mod_mul(u64 a, u64 b, u64 mod) {
    return (u64)((__uint128_t)a * b % mod);
}

static u64 mod_pow(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        d >>= 1;
    }
    return r;
}

static bool is_prime_32(u32 n) {
    if (n < 2) return false;
    static const u32 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u32 p : small) {
        if (n % p == 0) return n == p;
    }
    u32 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    for (u32 a : {2u, 7u, 61u}) {
        if (a % n == 0) continue;
        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul(x, x, n);
            if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
    }
    return true;
}
// Block 2: Tonelli-Shanks & Cornacchia Solver (R = e^2 + f^2)
static u64 tonelli(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return n;
    if (mod_pow(n, (p - 1) / 2, p) != 1) return 0;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);

    u64 q = p - 1;
    int s = 0;
    while ((q & 1) == 0) { q >>= 1; ++s; }

    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;
    u64 c = mod_pow(z, q, p);
    u64 x = mod_pow(n, (q + 1) / 2, p);
    u64 t = mod_pow(n, q, p);
    int m = s;

    while (t != 1) {
        u64 t2 = mod_mul(t, t, p);
        int i = 1;
        while (i < m) {
            if (t2 == 1) break;
            t2 = mod_mul(t2, t2, p);
            ++i;
        }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        x = mod_mul(x, b, p);
        c = mod_mul(b, b, p);
        t = mod_mul(t, c, p);
        m = i;
    }
    return x;
}

static bool sum2_prime(u32 R, u32 &e, u32 &f) {
    if (R % 4 != 1) return false;
    u64 r = tonelli(R - 1, R);
    if (r == 0) return false;

    u64 aa = R, bb = r;
    while (bb * bb > R) {
        u64 tmp = aa % bb;
        aa = bb;
        bb = tmp;
    }
    e = (u32)bb;
    u64 rem = (u64)R - (u64)e * e;
    u64 ff = (u64)sqrt((long double)rem);
    while ((ff + 1) * (ff + 1) <= rem) ++ff;
    while (ff * ff > rem) --ff;
    if (ff * ff != rem) return false;
    f = (u32)ff;
    return e > 0 && f > 0;
}
// Block 3: Exact Haar-Measure S^3+ Prime Solver + S4 (24-Element) Symmetry Orbit
// Plastic constant R2 quasi-random generator
static inline pair<double, double> plastic_r2(size_t n) {
    const double phi2 = 1.32471795724474602596;
    const double a1 = 1.0 / phi2;
    const double a2 = 1.0 / (phi2 * phi2);
    double u1 = fmod(0.5 + a1 * (double)(n + 1), 1.0);
    double u2 = fmod(0.5 + a2 * (double)(n + 1), 1.0);
    return {u1, u2};
}

// Inverse Beta(1/2, 1/2) CDF: maps u in [0,1] -> sin^2(pi/2 * u)
static inline double inv_beta_half(double u) {
    double s = sin(0.5 * M_PI * u);
    return s * s;
}

// Finds p = c^2 + d^2 + e^2 + f^2 steered to exact Haar target (u1, u2) on S^3_+
static bool haar_4square(u32 p, size_t step_idx, array<double, 4> &s3_pt) {
    auto [u1, u2] = plastic_r2(step_idx);

    // On S^3_+, r^2 = x1^2 + x2^2 ~ Beta(1, 1) = Uniform(0, 1)
    // Inner split x1^2 / r^2 ~ Beta(1/2, 1/2)
    double r2 = 0.02 + 0.96 * u1; // keep interior away from degenerate 0-boundary
    double split = inv_beta_half(u2);

    double target_c = sqrt((double)p * r2 * split);
    double target_d = sqrt((double)p * r2 * (1.0 - split));

    i64 c0 = (i64)llround(target_c);
    i64 d0 = (i64)llround(target_d);

    // Enforce parity so c^2 + d^2 = 0 mod 4 (if p=1 mod 4) or 2 mod 4 (if p=3 mod 4)
    if (p % 4 == 1) {
        if (c0 % 2 != 0) c0 += (c0 < target_c ? 1 : -1);
        if (d0 % 2 != 0) d0 += (d0 < target_d ? 1 : -1);
    } else {
        if (c0 % 2 == 0) c0 += (c0 < target_c ? 1 : -1);
        if (d0 % 2 == 0) d0 += (d0 < target_d ? 1 : -1);
    }
    if (c0 <= 0) c0 += 2;
    if (d0 <= 0) d0 += 2;

    // Expanding diamond spiral around (c0, d0) in steps of 2
    for (int rad = 0; rad <= 60; ++rad) {
        for (int dc = -rad; dc <= rad; ++dc) {
            int dd_abs = rad - abs(dc);
            int signs[2] = {-dd_abs, dd_abs};
            int n_signs = (dd_abs == 0) ? 1 : 2;
            for (int si = 0; si < n_signs; ++si) {
                i64 c = c0 + 2LL * dc;
                i64 d = d0 + 2LL * signs[si];
                if (c <= 0 || d <= 0) continue;

                i64 Rll = (i64)p - c * c - d * d;
                if (Rll <= 4 || Rll % 4 != 1) continue;

                u32 R = (u32)Rll;
                if (!is_prime_32(R)) continue;

                u32 e, f;
                if (!sum2_prime(R, e, f)) continue;

                // Apply 24-element S4 coordinate permutation group (Hurwitz octant symmetry)
                // Removes Cornacchia's e < f ordering bias
                array<double, 4> raw = {
                    (double)c / sqrt((double)p),
                    (double)d / sqrt((double)p),
                    (double)e / sqrt((double)p),
                    (double)f / sqrt((double)p)
                };
                array<int, 4> perm = {0, 1, 2, 3};
                size_t orbit_id = step_idx % 24;
                for (size_t k = 0; k < orbit_id; ++k) {
                    next_permutation(perm.begin(), perm.end());
                }
                for (int j = 0; j < 4; ++j) s3_pt[j] = raw[perm[j]];
                return true;
            }
        }
    }
    return false;
}
// Block 4: Riesz s=2 Energy Benchmark on S^3_+ (Haar Prime Walk vs Halton vs Monte Carlo)
static double halton_seq(size_t index, int base) {
    double f = 1.0, r = 0.0;
    size_t i = index + 1;
    while (i > 0) {
        f /= (double)base;
        r += f * (double)(i % base);
        i /= base;
    }
    return r;
}

// Maps 3D [0,1]^3 -> exact Haar uniform point on positive 3-sphere S^3_+
static array<double, 4> map_cube_to_s3_plus(double u1, double u2, double u3) {
    // Dirichlet(1/2, 1/2, 1/2, 1/2) via stick-breaking Beta variables:
    // w1 ~ Beta(1/2, 3/2), w2 ~ Beta(1/2, 1), w3 ~ Beta(1/2, 1/2)
    // Equivalently, two independent Beta(1/2, 1/2) angles + Uniform(0,1) radius^2
    double r2 = u1;
    double a1 = inv_beta_half(u2);
    double a2 = inv_beta_half(u3);
    return {
        sqrt(r2 * a1),
        sqrt(r2 * (1.0 - a1)),
        sqrt((1.0 - r2) * a2),
        sqrt((1.0 - r2) * (1.0 - a2))
    };
}

// Computes normalized Riesz s=2 energy: (2 / N(N-1)) * sum_{i < j} 1 / ||x_i - x_j||^2
static double riesz_s2_energy(const vector<array<double, 4>> &pts) {
    size_t N = pts.size();
    double total_energy = 0.0;
    #pragma omp parallel for reduction(+:total_energy) schedule(dynamic, 16)
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double d2 = 0.0;
            for (int k = 0; k < 4; ++k) {
                double diff = pts[i][k] - pts[j][k];
                d2 += diff * diff;
            }
            if (d2 > 1e-18) total_energy += 1.0 / d2;
        }
    }
    return (2.0 * total_energy) / ((double)N * (double)(N - 1));
}

int main() {
    const size_t TARGET_N = 1024;
    vector<array<double, 4>> prime_pts;
    prime_pts.reserve(TARGET_N);

    // Stream primes > 10,000 and construct Haar-mapped S^3_+ points
    size_t accepted = 0;
    for (u32 p = 10007; accepted < TARGET_N; p += 2) {
        if (!is_prime_32(p)) continue;
        array<double, 4> pt;
        if (haar_4square(p, accepted, pt)) {
            prime_pts.push_back(pt);
            accepted++;
        }
    }

    size_t N = prime_pts.size();
    vector<array<double, 4>> halton_pts(N), mc_pts(N);
    mt19937_64 rng(42);
    uniform_real_distribution<double> unif(0.0, 1.0);

    for (size_t i = 0; i < N; ++i) {
        halton_pts[i] = map_cube_to_s3_plus(
            halton_seq(i, 2), halton_seq(i, 3), halton_seq(i, 5)
        );
        mc_pts[i] = map_cube_to_s3_plus(unif(rng), unif(rng), unif(rng));
    }

    double e_prime  = riesz_s2_energy(prime_pts);
    double e_halton = riesz_s2_energy(halton_pts);
    double e_mc     = riesz_s2_energy(mc_pts);

    cout << fixed << setprecision(6);
    cout << "=== S^3_+ HAAR-MEASURE QMC BENCHMARK (N = " << N << ") ===\n";
    cout << "Halton QMC (Exact Haar)     : Riesz s=2 Energy = " << e_halton << "\n";
    cout << "Haar Prime Walk (S4 Orbit)  : Riesz s=2 Energy = " << e_prime  << "\n";
    cout << "Monte Carlo (Exact Haar)    : Riesz s=2 Energy = " << e_mc     << "\n";
    cout << "Ratio (Prime Walk / Halton) : " << (e_prime / e_halton) << "x\n";
    return 0;
}
