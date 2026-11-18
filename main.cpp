// Block 1: Headers & Fast 32-bit Modular Number Theory
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
// Block 3: Unbiased Full-Range Haar Prime Walk & Corrected SO(3) Integrands
static inline pair<double, double> plastic_r2(size_t n) {
    const double phi2 = 1.32471795724474602596;
    const double a1 = 1.0 / phi2;
    const double a2 = 1.0 / (phi2 * phi2);
    return {
        fmod(0.5 + a1 * (double)(n + 1), 1.0),
        fmod(0.5 + a2 * (double)(n + 1), 1.0)
    };
}

static inline double inv_beta_half(double u) {
    double s = sin(0.5 * M_PI * u);
    return s * s;
}

// Unbiased Haar 4-square decomposition: full [0,1] range, zero cap clipping, wide fallback
static bool haar_4square_unbiased(u32 p, size_t step_idx, array<double, 4> &s3_pt) {
    auto [u1, u2] = plastic_r2(step_idx);
    double r2 = u1; // Full radial coverage [0, 1] (no cap truncation)
    double split = inv_beta_half(u2);

    double target_c = sqrt((double)p * r2 * split);
    double target_d = sqrt((double)p * r2 * (1.0 - split));

    i64 c0 = (i64)llround(target_c);
    i64 d0 = (i64)llround(target_d);

    if (p % 4 == 1) {
        if (c0 % 2 != 0) c0 += (c0 < target_c ? 1 : -1);
        if (d0 % 2 != 0) d0 += (d0 < target_d ? 1 : -1);
    } else {
        if (c0 % 2 == 0) c0 += (c0 < target_c ? 1 : -1);
        if (d0 % 2 == 0) d0 += (d0 < target_d ? 1 : -1);
    }
    if (c0 <= 0) c0 += 2;
    if (d0 <= 0) d0 += 2;

    // Wide diamond spiral (rad <= 250) ensures virtually 100% acceptance without boundary bias
    for (int rad = 0; rad <= 250; ++rad) {
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

                double inv_sqp = 1.0 / sqrt((double)p);
                array<double, 4> raw = {
                    (double)c * inv_sqp, (double)d * inv_sqp,
                    (double)e * inv_sqp, (double)f * inv_sqp
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

// Exact Haar SO(3) integrals: E[Tr(R)] = 0.0, E[Tr(R^2)] = 1.0, E[v^T R^T M R v] = 2.0
static array<double, 3> eval_so3_octant_avg(const array<double, 4> &q_pos) {
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    const double vx = 1.0 / sqrt(14.0), vy = 2.0 / sqrt(14.0), vz = 3.0 / sqrt(14.0);

    for (int mask = 0; mask < 8; ++mask) {
        double w = q_pos[0];
        double x = (mask & 1) ? -q_pos[1] : q_pos[1];
        double y = (mask & 2) ? -q_pos[2] : q_pos[2];
        double z = (mask & 4) ? -q_pos[3] : q_pos[3];

        double r00 = 1.0 - 2.0 * (y*y + z*z), r01 = 2.0 * (x*y - z*w),       r02 = 2.0 * (x*z + y*w);
        double r10 = 2.0 * (x*y + z*w),       r11 = 1.0 - 2.0 * (x*x + z*z), r12 = 2.0 * (y*z - x*w);
        double r20 = 2.0 * (x*z - y*w),       r21 = 2.0 * (y*z + x*w),       r22 = 1.0 - 2.0 * (x*x + y*y);

        double trR = r00 + r11 + r22;
        sum1 += trR;

        double trR2 = (r00*r00 + r01*r10 + r02*r20)
                    + (r10*r01 + r11*r11 + r12*r21)
                    + (r20*r02 + r21*r12 + r22*r22);
        sum2 += trR2;

        double u0 = r00*vx + r01*vy + r02*vz;
        double u1 = r10*vx + r11*vy + r12*vz;
        double u2 = r20*vx + r21*vy + r22*vz;
        sum3 += 1.0 * u0 * u0 + 2.0 * u1 * u1 + 3.0 * u2 * u2;
    }
    return {sum1 / 8.0, sum2 / 8.0, sum3 / 8.0};
}
// Block 4: 3D Sobol' + Halton + MC Benchmark with Corrected Exact Integrals {0.0, 1.0, 2.0}
struct Sobol3D {
    uint32_t dir[3][32];
    Sobol3D() {
        // Dim 1: van der Corput base 2
        for (int i = 0; i < 32; ++i) dir[0][i] = 1U << (31 - i);
        // Dim 2: primitive poly x^2 + x + 1 (deg 2, coeff 1), m = {1, 3}
        uint32_t m2[2] = {1, 3};
        for (int i = 0; i < 2; ++i) dir[1][i] = m2[i] << (31 - i);
        for (int i = 2; i < 32; ++i) {
            dir[1][i] = dir[1][i - 2] ^ (dir[1][i - 2] >> 2) ^ dir[1][i - 1];
        }
        // Dim 3: primitive poly x^3 + x + 1 (deg 3, coeff 1), m = {1, 3, 7}
        uint32_t m3[3] = {1, 3, 7};
        for (int i = 0; i < 3; ++i) dir[2][i] = m3[i] << (31 - i);
        for (int i = 3; i < 32; ++i) {
            dir[2][i] = dir[2][i - 3] ^ (dir[2][i - 3] >> 3) ^ dir[2][i - 2];
        }
    }
    array<double, 3> point(uint32_t idx) const {
        uint32_t g = (idx + 1) ^ ((idx + 1) >> 1);
        uint32_t x = 0, y = 0, z = 0;
        for (int bit = 0; bit < 32; ++bit) {
            if (g & (1U << bit)) {
                x ^= dir[0][bit];
                y ^= dir[1][bit];
                z ^= dir[2][bit];
            }
        }
        const double inv2_32 = 1.0 / 4294967296.0;
        return {x * inv2_32, y * inv2_32, z * inv2_32};
    }
};

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

static array<double, 4> map_cube_to_s3_plus(double u1, double u2, double u3) {
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

int main() {
    const size_t MAX_N = 2048;
    vector<array<double, 4>> prime_pts, sobol_pts, halton_pts;
    prime_pts.reserve(MAX_N);
    sobol_pts.reserve(MAX_N);
    halton_pts.reserve(MAX_N);

    Sobol3D sobol;
    size_t accepted = 0, attempted = 0;

    // Start at p > 10^7 (sqrt(p) > 3162) so integer parity quantization is < 0.03%
    for (u32 p = 10000019; accepted < MAX_N; p += 2) {
        if (!is_prime_32(p)) continue;
        attempted++;
        array<double, 4> pt;
        if (haar_4square_unbiased(p, accepted, pt)) {
            prime_pts.push_back(pt);
            auto sb = sobol.point((uint32_t)accepted);
            sobol_pts.push_back(map_cube_to_s3_plus(sb[0], sb[1], sb[2]));
            halton_pts.push_back(map_cube_to_s3_plus(
                halton_seq(accepted, 2), halton_seq(accepted, 3), halton_seq(accepted, 5)
            ));
            accepted++;
        }
    }

    cout << "Primes attempted: " << attempted << " | Accepted: " << accepted
         << " (" << (100.0 * accepted / attempted) << "% acceptance)\n\n";

    // 1. Detailed 3-Function Comparison at N = 1024 with Corrected Exact Integrals
    const size_t N_REF = 1024;
    const int MC_TRIALS = 64;
    array<double, 3> p_sum = {0, 0, 0}, s_sum = {0, 0, 0}, h_sum = {0, 0, 0};

    for (size_t i = 0; i < N_REF; ++i) {
        auto vp = eval_so3_octant_avg(prime_pts[i]);
        auto vs = eval_so3_octant_avg(sobol_pts[i]);
        auto vh = eval_so3_octant_avg(halton_pts[i]);
        for (int f = 0; f < 3; ++f) {
            p_sum[f] += vp[f];
            s_sum[f] += vs[f];
            h_sum[f] += vh[f];
        }
    }

    array<double, 3> mc_sq_err = {0, 0, 0};
    array<double, 3> exact = {0.0, 1.0, 2.0}; // Corrected: E[Tr(R)] = 0, E[Tr(R^2)] = 1, E[f3] = 2
    for (int t = 0; t < MC_TRIALS; ++t) {
        mt19937_64 rng(1000 + t);
        uniform_real_distribution<double> unif(0.0, 1.0);
        array<double, 3> m_sum = {0, 0, 0};
        for (size_t i = 0; i < N_REF; ++i) {
            auto vm = eval_so3_octant_avg(map_cube_to_s3_plus(unif(rng), unif(rng), unif(rng)));
            for (int f = 0; f < 3; ++f) m_sum[f] += vm[f];
        }
        for (int f = 0; f < 3; ++f) {
            double err = (m_sum[f] / N_REF) - exact[f];
            mc_sq_err[f] += err * err;
        }
    }

    cout << fixed << setprecision(6);
    cout << "=== UNBIASED SO(3) INTEGRATION ERROR (N = " << N_REF << " Base Samples) ===\n";
    const char* names[3] = {"f1(R) = Tr(R)       ", "f2(R) = Tr(R^2)     ", "f3(R) = v^T R^T M Rv"};
    for (int f = 0; f < 3; ++f) {
        double err_p  = fabs((p_sum[f] / N_REF) - exact[f]);
        double err_s  = fabs((s_sum[f] / N_REF) - exact[f]);
        double err_h  = fabs((h_sum[f] / N_REF) - exact[f]);
        double err_mc = sqrt(mc_sq_err[f] / MC_TRIALS);
        cout << names[f]
             << " | Sobol: " << err_s
             << " | Halton: " << err_h
             << " | Prime Walk: " << err_p
             << " | MC RMSE: " << err_mc << "\n";
    }

    // 2. Convergence Scaling on f3(R) from N = 128 to 2048
    cout << "\n=== UNBIASED CONVERGENCE SCALING ON f3(R) (N = 128 -> 2048) ===\n";
    vector<size_t> scales = {128, 256, 512, 1024, 2048};
    vector<double> logN, logEp, logEs, logEh, logEm;

    for (size_t N : scales) {
        double sp = 0.0, ss = 0.0, sh = 0.0;
        for (size_t i = 0; i < N; ++i) {
            sp += eval_so3_octant_avg(prime_pts[i])[2];
            ss += eval_so3_octant_avg(sobol_pts[i])[2];
            sh += eval_so3_octant_avg(halton_pts[i])[2];
        }
        double ep = fabs((sp / N) - 2.0);
        double es = fabs((ss / N) - 2.0);
        double eh = fabs((sh / N) - 2.0);

        double m_sq = 0.0;
        for (int t = 0; t < MC_TRIALS; ++t) {
            mt19937_64 rng(5000 + t);
            uniform_real_distribution<double> unif(0.0, 1.0);
            double sm = 0.0;
            for (size_t i = 0; i < N; ++i) {
                sm += eval_so3_octant_avg(map_cube_to_s3_plus(unif(rng), unif(rng), unif(rng)))[2];
            }
            double e = (sm / N) - 2.0;
            m_sq += e * e;
        }
        double em = sqrt(m_sq / MC_TRIALS);

        cout << "N = " << setw(4) << N
             << " | Sobol: " << es
             << " | Halton: " << eh
             << " | Prime Walk: " << ep
             << " | MC RMSE: " << em << "\n";

        logN.push_back(log((double)N));
        logEs.push_back(log(max(1e-12, es)));
        logEh.push_back(log(max(1e-12, eh)));
        logEp.push_back(log(max(1e-12, ep)));
        logEm.push_back(log(max(1e-12, em)));
    }

    auto calc_slope = [&](const vector<double> &y) {
        double mx = 0, my = 0;
        for (size_t i = 0; i < logN.size(); ++i) { mx += logN[i]; my += y[i]; }
        mx /= logN.size(); my /= logN.size();
        double num = 0, den = 0;
        for (size_t i = 0; i < logN.size(); ++i) {
            num += (logN[i] - mx) * (y[i] - my);
            den += (logN[i] - mx) * (logN[i] - mx);
        }
        return -num / den;
    };

    cout << "\nEmpirical Convergence Exponent alpha (Error ~ O(N^-alpha)):\n";
    cout << "Sobol' QMC Slope      : alpha = " << calc_slope(logEs) << "\n";
    cout << "Halton QMC Slope      : alpha = " << calc_slope(logEh) << "\n";
    cout << "Haar Prime Walk Slope : alpha = " << calc_slope(logEp) << "\n";
    cout << "Monte Carlo Slope     : alpha = " << calc_slope(logEm) << "\n";
    return 0;
}
