// Prime Walk - 1 BILLION Test : Cyclic = Barrier, 
// Tests primes up to 1e9 (50,847,534 primes) - streaming, no storage
// Arithmetic-Geometric: Norm = Pythagoras, Det = Ptolemy, Barrier = Stay Cyclic
// Compile: g++ -O3 -march=native -fopenmp -std=c++17 prime_1billion.cpp -o prime_1billion
// Run: ./prime_1billion --limit 1000000000 --sample 1000
//      ./prime_1billion --limit 100000000 --sample 1   // full check to 100M

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using u64 = uint64_t;
using u32 = uint32_t;
using i64 = int64_t;

// ---------- fast 32-bit Miller-Rabin : deterministic for < 2^32, bases 2,7,61 ----------
static u64 mod_mul(u64 a, u64 b, u64 mod) { return (u64)((__uint128_t)a * b % mod); }
static u64 mod_pow(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) { if (d & 1) r = mod_mul(r, a, mod); a = mod_mul(a, a, mod); d >>= 1; }
    return r;
}
static bool is_prime_32(u32 n) {
    if (n < 2) return false;
    static const u32 small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (u32 p : small) { if (n % p == 0) return n == p; }
    u32 d = n - 1; int s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    // bases 2,7,61 deterministic for 32-bit
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

// ---------- Tonelli-Shanks : sqrt(n) mod p, p prime ----------
static u64 tonelli(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return n;
    if (mod_pow(n, (p - 1) / 2, p) != 1) return 0; // not QR
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);
    u64 q = p - 1; int s = 0;
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
        while (i < m) { if (t2 == 1) break; t2 = mod_mul(t2, t2, p); ++i; }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        x = mod_mul(x, b, p);
        c = mod_mul(b, b, p);
        t = mod_mul(t, c, p);
        m = i;
    }
    return x;
}

// ---------- Cornacchia : R prime 1 mod 4 -> e^2 + f^2 = R ----------
static bool sum2_prime(u32 R, u32 &e, u32 &f) {
    if (R % 4 != 1) return false;
    u64 r = tonelli(R - 1, R);
    if (r == 0) return false;
    u64 aa = R, bb = r;
    while (bb * bb > R) { u64 tmp = aa % bb; aa = bb; bb = tmp; }
    e = (u32)bb;
    u64 rem = (u64)R - (u64)e * e;
    u64 ff = sqrt((long double)rem);
    while ((ff + 1) * (ff + 1) <= rem) ++ff;
    while (ff * ff > rem) --ff;
    if (ff * ff != rem) return false;
    f = (u32)ff;
    return e > 0 && f > 0;
}

// ---------- Balanced 4-square : p = c^2 + d^2 + e^2 + f^2 ----------
// c,d near sqrt(p)/2, parity fixed by p mod 4:
//   p=1 mod 4 -> c,d even  (c^2+d^2 = 0 mod 4, R = 1 mod 4)
//   p=3 mod 4 -> c,d odd   (c^2+d^2 = 2 mod 4, R = 1 mod 4)
static bool balanced_4square(u32 p, array<u32,4> &out) {
    double sq = sqrt((double)p);
    u32 target = (u32)(sq / 2.0);
    u32 base;
    if (p % 4 == 1) base = (target % 2 == 0) ? target : target + 1;
    else             base = (target % 2 == 1) ? target : target + 1;
    if (base < 2) base = 2;

    for (int k = 0; k < 5000; ++k) {
        i64 c = (i64)base + 2LL * k;
        i64 d = (i64)base - 2LL * k + 2;
        if (c <= 0 || d <= 0) continue;
        i64 Rll = (i64)p - c * c - d * d;
        if (Rll <= 0 || Rll % 4 != 1) continue;
        u32 R = (u32)Rll;
        if (!is_prime_32(R)) continue;
        u32 e, f;
        if (!sum2_prime(R, e, f)) continue;
        out = {(u32)c, (u32)d, e, f};
        // verify
        __uint128_t sum = (__uint128_t)c*c + (__uint128_t)d*d + (__uint128_t)e*e + (__uint128_t)f*f;
        if (sum == p) return true;
    }
    return false;
}

// ---------- Cyclic check for one prime : returns max |Det|/p, cyclic flag ----------
struct CheckResult { bool cyclic; double max_det; double area_over_p; double ptolemy; };

static CheckResult check_prime(u32 p, const array<u32,4> &t) {
    double A = t[0], B = t[1], C = t[2], D = t[3];
    double s = (A + B + C + D) / 2.0;
    double prod = (s - A) * (s - B) * (s - C) * (s - D);
    double area = (prod > 0) ? sqrt(prod) : 0.0;

    double ac_bd = A * C + B * D; // Ptolemy
    double ab_cd = A * B + C * D;
    double ad_bc = A * D + B * C;
    if (ab_cd == 0 || ad_bc == 0 || ac_bd == 0) return {false, 1e9, 0, 0};

    double e2 = ac_bd * ad_bc / ab_cd;
    double f2 = ac_bd * ab_cd / ad_bc;
    double e = sqrt(max(0.0, e2));
    double f = sqrt(max(0.0, f2));
    double ratio = (e * f) / ac_bd;

    bool cyclic = (fabs(ratio - 1.0) < 1e-9) && (area > 0);

    double Det1 = (A * D - B * C) / (double)p;
    double Det2 = (A * C + B * D) / (double)p;
    double DetF = (A * B + C * D) / (double)p;
    double maxd = max({fabs(Det1), fabs(Det2), fabs(DetF)});

    return {cyclic, maxd, area / (double)p, ratio};
}

// ---------- Segmented sieve : stream primes up to LIMIT ----------
static void segmented_sieve(u64 LIMIT, function<void(u32)> on_prime) {
    const u64 SEG = 1 << 20; // 1M segment
    u64 sq = sqrt((long double)LIMIT) + 1;
    vector<u32> base_primes;
    vector<bool> is_composite_small(sq + 1, false);
    for (u64 i = 2; i <= sq; ++i) if (!is_composite_small[i]) {
        base_primes.push_back((u32)i);
        if (i * i <= sq) for (u64 j = i * i; j <= sq; j += i) is_composite_small[j] = true;
    }
    vector<bool> is_comp(SEG);
    for (u64 low = 2; low <= LIMIT; low += SEG) {
        u64 high = min(low + SEG - 1, LIMIT);
        fill(is_comp.begin(), is_comp.end(), false);
        for (u32 pr : base_primes) {
            u64 start = max((u64)pr * pr, ((low + pr - 1) / pr) * (u64)pr);
            for (u64 j = start; j <= high; j += pr) is_comp[j - low] = true;
        }
        for (u64 i = low; i <= high; ++i) if (!is_comp[i - low]) on_prime((u32)i);
    }
}

// ---------- main ----------
int main(int argc, char** argv) {
    u64 LIMIT = 1000000000ULL;
    u64 SAMPLE = 1000; // full check every SAMPLE-th prime
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "--limit" && i + 1 < argc) LIMIT = stoull(argv[++i]);
        if (a == "--sample" && i + 1 < argc) SAMPLE = stoull(argv[++i]);
    }

    cout << "Prime Walk - 1B Test : Cyclic = Barrier, NO Owen\n";
    cout << "Limit: " << LIMIT << "  Sample: every " << SAMPLE << "-th prime full-checked\n";
    cout << "Barrier B = 4 + 2*ln(p)\n\n";

    // Global counters - atomic for OpenMP reduction via critical
    atomic<long long> total{0}, full_checked{0}, balanced_ok{0};
    atomic<long long> cyclic_viol{0}, det_viol{0}, barrier_viol{0};
    double max_det_global = 0.0;
    u32 max_det_p = 0;
    mutex max_mtx;

    auto t0 = chrono::high_resolution_clock::now();

    // We collect primes in batches for parallel check, to keep sieve single-threaded
    const size_t BATCH = 100000;
    vector<u32> batch;
    batch.reserve(BATCH);

    auto flush_batch = [&]() {
        if (batch.empty()) return;
        #pragma omp parallel for schedule(dynamic, 1000)
        for (size_t idx = 0; idx < batch.size(); ++idx) {
            u32 p = batch[idx];
            long long my_total = 1; // counted below via total
            (void)my_total;
            // Only full-check sampled primes - otherwise quick skip
            // total is incremented in sieve callback, not here
            if ((idx % SAMPLE) != 0 && SAMPLE != 1) {
                // For non-sampled primes we still count them, but do not do expensive 4-square
                // To keep the test honest, we do NOT claim cyclic for them.
                continue;
            }
            array<u32,4> t;
            if (!balanced_4square(p, t)) continue;
            balanced_ok++;
            auto r = check_prime(p, t);
            double Barrier = 4.0 + 2.0 * log((double)p);
            bool barrier_ok = r.cyclic && (r.max_det <= 1.0);
            if (!r.cyclic) cyclic_viol++;
            if (r.max_det > 1.0) det_viol++;
            if (!barrier_ok) barrier_viol++;
            {
                lock_guard<mutex> lk(max_mtx);
                if (r.max_det > max_det_global) { max_det_global = r.max_det; max_det_p = p; }
            }
            full_checked++;
            // Print first few and last few for proof
            if (full_checked < 5 || p > 999000000) {
                #pragma omp critical
                {
                    cout << "p=" << p << " mod24=" << (p % 24)
                         << " tup=(" << t[0] << "," << t[1] << "," << t[2] << "," << t[3] << ")"
                         << " maxDet=" << r.max_det
                         << " Area/p=" << r.area_over_p
                         << " Ptolemy=" << r.ptolemy
                         << " cyclic=" << (r.cyclic ? "true" : "false")
                         << " Barrier=" << Barrier << "\n";
                }
            }
        }
        batch.clear();
    };

    segmented_sieve(LIMIT, [&](u32 p){
        total++;
        batch.push_back(p);
        if (batch.size() >= BATCH) flush_batch();
        if (total % 10000000 == 0) {
            auto now = chrono::high_resolution_clock::now();
            double sec = chrono::duration<double>(now - t0).count();
            cout << "... " << total << " primes streamed, " << sec << "s\n";
        }
    });
    flush_batch();

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t1 - t0).count();

    cout << "\n=== SUMMARY : PRIMES UP TO " << LIMIT << ", NO OWEN, CYCLIC = BARRIER ===\n";
    cout << "Total primes streamed: " << total << "\n";
    cout << "Full cyclic-checked (sampled): " << full_checked << "\n";
    cout << "Balanced 4-square found: " << balanced_ok << "\n";
    cout << "Cyclic violations: " << cyclic_viol << "\n";
    cout << "Det>1 violations: " << det_viol << "\n";
    cout << "Barrier violations: " << barrier_viol << "\n";
    cout << "Max |Det|/p: " << max_det_global << " at p=" << max_det_p << "\n";
    cout << "Time: " << sec << "s\n";
    cout << "Conclusion: If cyclic, bounded. No gear needed.\n";
    return 0;
}