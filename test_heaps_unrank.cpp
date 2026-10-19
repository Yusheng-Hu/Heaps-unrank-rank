// =============================================================================
//  test_heaps_unrank.cpp -- verification suite for the heaps_unrank library
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Both implementations are exercised side by side:
//
//    main     (Theta(n^2))  heaps_unrank.hpp       heaps::unrank      / heaps::rank
//    baseline (Theta(n^3))  heaps_unrank_o_n3.hpp  heaps::unrank_o_n3 / heaps::rank_o_n3
//
//  The suite checks them against independent ground truth and against each
//  other:
//
//   1. step_map(k) equals the closed form F_k of the paper (Theorem 5.2),
//      for k = 1..24.
//   2. The step table has the structure the main implementation relies on
//      (Corollary following the cycle-dichotomy lemma): for odd k it is the
//      transposition (0, k-1); for even k it is a single k-cycle that follows
//      the cycle order [0, k-3, ..., 1, k-2, k-1].
//   3. Exhaustive comparison with the iterative Heap algorithm (ground truth)
//      for every code of order n = 1..9 -- 362,880 codes at n = 9 -- for BOTH
//      implementations.  With --full, n = 10 (3,628,800 codes) is included.
//   4. Exhaustive round trip rank(unrank(C)) == C for every code of order
//      n = 1..8, for BOTH implementations, plus the agreement
//      unrank(C) == unrank_o_n3(C) on every code.
//   5. Randomised round trips at larger orders (n up to 100), for both
//      implementations, together with a validity check that unrank really
//      returns a permutation of 0..n-1.
//
//  Exit status is 0 when every check passes, 1 when a check fails, 2 on a
//  usage error.
//
//  Usage:  test_heaps_unrank [--quick | --full | --help]
//            (default)  checks 1-5 with n <= 9        (~20 s)
//            --quick    smaller random sample          (~3 s)
//            --full     additionally n = 10           (~40 s)
// =============================================================================

#include "heaps_unrank.hpp"
#include "heaps_unrank_o_n3.hpp"

#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

int g_checks = 0;
int g_failures = 0;

void report(bool ok, const std::string& what) {
    g_checks++;
    if (!ok) g_failures++;
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
}

void section(const char* title) {
    std::printf("\n%s\n", title);
}

// ---- ground truth: the iterative Heap algorithm, exposed as a stream --------
class HeapStream {
public:
    explicit HeapStream(int n) : A_((size_t)n), c_((size_t)n, 0), n_(n), i_(0) {
        for (int i = 0; i < n; i++) A_[(size_t)i] = i;
    }
    const std::vector<int>& current() const { return A_; }
    // Advance to the next permutation in Heap's order (no-op at the end).
    void advance() {
        if (i_ >= n_) return;
        for (;;) {
            if (c_[(size_t)i_] < i_) {
                if (i_ % 2 == 0) std::swap(A_[0], A_[(size_t)i_]);
                else             std::swap(A_[(size_t)c_[(size_t)i_]], A_[(size_t)i_]);
                c_[(size_t)i_]++;
                i_ = 0;
                return;
            }
            c_[(size_t)i_] = 0;
            i_++;
            if (i_ >= n_) return;
        }
    }

private:
    std::vector<int> A_;
    std::vector<int> c_;
    int n_;
    int i_;
};

// ---- the paper's closed form for the final arrangement F_m (Theorem 5.2) ----
std::vector<int> final_arrangement(int m) {
    if (m == 1) return std::vector<int>{0};
    if (m == 2) return std::vector<int>{1, 0};
    if (m == 3) return std::vector<int>{2, 1, 0};
    std::vector<int> f((size_t)m);
    if (m % 2 != 0) {                 // odd : F_m = [m-1, 1, 2, ..., m-2, 0]
        f[0] = m - 1;
        for (int i = 1; i < m - 1; i++) f[(size_t)i] = i;
        f[(size_t)(m - 1)] = 0;
    } else {                          // even: F_m = [m-3, m-2, 1, ..., m-4, m-1, 0]
        f[0] = m - 3;
        f[1] = m - 2;
        for (int i = 2; i < m - 2; i++) f[(size_t)i] = i - 1;
        f[(size_t)(m - 2)] = m - 1;
        f[(size_t)(m - 1)] = 0;
    }
    return f;
}

// ---- the cycle order used by the main implementation (even order only) ------
// cyc = [0, m-3, m-4, ..., 1, m-2, m-1];  cycle_pos: index -> position
int cycle_pos(int m, int i) {
    if (i == 0) return 0;
    if (i <= m - 3) return m - 2 - i;
    return i;
}
// cycle_index: position -> index in the cycle order (inverse of cycle_pos)
int cycle_index(int m, int p) {
    if (p == 0) return 0;
    if (p <= m - 3) return m - 2 - p;
    return p;
}

// ---- high-first Lehmer code counter: the lowest digit moves first ----------
void next_code(std::vector<int>& C, int n) {
    if (n < 2) return;                       // n = 1 has a single code
    int j = n - 2;
    C[(size_t)j]++;
    while (j > 0 && C[(size_t)j] > n - j - 1) { C[(size_t)j] = 0; j--; C[(size_t)j]++; }
}

long long factorial(int n) {
    long long f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

bool is_permutation(const std::vector<int>& p, int n) {
    std::vector<char> seen((size_t)n, 0);
    for (int i = 0; i < n; i++) {
        if (p[(size_t)i] < 0 || p[(size_t)i] >= n) return false;
        if (seen[(size_t)p[(size_t)i]]) return false;
        seen[(size_t)p[(size_t)i]] = 1;
    }
    return true;
}

// =============================================================================
// 1. step_map against the closed form
// =============================================================================
void test_step_map_closed_form() {
    section("1. step_map(k) against the closed form F_k (paper, Theorem 5.2)");
    bool ok = true;
    for (int k = 1; k <= 24; k++) {
        if (heaps::step_map(k) != final_arrangement(k)) {
            std::printf("    mismatch at k = %d\n", k);
            ok = false;
        }
    }
    report(ok, "step_map(k) == F_k for k = 1..24");
}

// =============================================================================
// 2. the structural dichotomy the main implementation relies on
// =============================================================================
void test_step_dichotomy() {
    section("2. step table: transposition (odd k) or one cycle in cycle order (even k)");
    bool ok = true;
    for (int m = 3; m <= 24; m++) {
        const std::vector<int> t = heaps::step_map(m);
        if (m % 2 != 0) {
            // exactly the transposition (0, m-1)
            for (int i = 0; i < m; i++) {
                const int want = (i == 0) ? m - 1 : (i == m - 1 ? 0 : i);
                if (t[(size_t)i] != want) { ok = false; }
            }
        } else {
            // one m-cycle following the cycle order: t[cyc[i]] == cyc[(i+1) mod m]
            for (int i = 0; i < m; i++) {
                const int from = cycle_pos(m, i);
                const int to = cycle_pos(m, (i + 1) % m);
                if (t[(size_t)from] != to) { ok = false; }
            }
            // and the two closed forms are mutually inverse
            for (int i = 0; i < m; i++) {
                if (cycle_index(m, cycle_pos(m, i)) != i) { ok = false; }
            }
        }
    }
    report(ok, "dichotomy of the step table holds for m = 3..24");
}

// =============================================================================
// 3. exhaustive: every code of order n against the Heap sequence, both versions
// =============================================================================
void test_exhaustive_against_sequence(int n) {
    HeapStream heap(n);
    std::vector<int> C((size_t)n, 0), dense((size_t)n), impl((size_t)n);
    const long long total = factorial(n);
    long long r = 0;
    bool ok = true;
    for (; r < total; r++, next_code(C, n)) {
        heaps::unrank_o_n3(n, C.data(), dense.data());
        heaps::unrank(n, C.data(), impl.data());
        if (dense != heap.current()) {
            std::printf("    baseline mismatch at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
        if (impl != heap.current()) {
            std::printf("    main mismatch at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
        if (impl != dense) {
            std::printf("    implementations disagree at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
        heap.advance();
    }
    report(ok, "n = " + std::to_string(n) + ": " + std::to_string(r) +
                   " codes match the iterative Heap sequence (both implementations)");
}

// =============================================================================
// 4. exhaustive round trip rank(unrank(C)) == C, both versions
// =============================================================================
void test_exhaustive_roundtrip(int n) {
    std::vector<int> C((size_t)n, 0), dense((size_t)n), impl((size_t)n);
    std::vector<int> back1((size_t)n), back2((size_t)n);
    const long long total = factorial(n);
    long long r = 0;
    bool ok = true;
    for (; r < total; r++, next_code(C, n)) {
        heaps::unrank_o_n3(n, C.data(), dense.data());
        heaps::unrank(n, C.data(), impl.data());
        if (dense != impl) {
            std::printf("    implementations disagree at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
        heaps::rank_o_n3(n, dense.data(), back1.data());
        heaps::rank(n, impl.data(), back2.data());
        if (C != back1 || C != back2) {
            std::printf("    round trip mismatch at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
    }
    report(ok, "n = " + std::to_string(n) + ": round trip for " + std::to_string(r) +
                   " codes (both implementations)");
}

// =============================================================================
// 5. randomised round trips at larger orders, both versions
// =============================================================================
void test_random_roundtrips(int trials) {
    section("5. randomised round trips at larger orders (fixed seed 20261019)");
    std::mt19937 rng(20261019u);
    const int sizes[] = {12, 16, 20, 32, 48, 64, 100};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
        const int n = sizes[s];
        std::vector<int> C((size_t)n), dense((size_t)n), impl((size_t)n);
        std::vector<int> back1((size_t)n), back2((size_t)n);
        bool ok = true;
        for (int t = 0; t < trials && ok; t++) {
            for (int j = 0; j <= n - 2; j++) C[(size_t)j] = (int)(rng() % (unsigned)(n - j));
            C[(size_t)(n - 1)] = 0;

            heaps::unrank_o_n3(n, C.data(), dense.data());
            heaps::unrank(n, C.data(), impl.data());
            if (!is_permutation(dense, n) || !is_permutation(impl, n)) {
                std::printf("    unrank returned a non-permutation at n = %d\n", n);
                ok = false;
                break;
            }
            if (dense != impl) {
                std::printf("    implementations disagree at n = %d\n", n);
                ok = false;
                break;
            }
            heaps::rank_o_n3(n, dense.data(), back1.data());
            heaps::rank(n, impl.data(), back2.data());
            if (C != back1 || C != back2) {
                std::printf("    round trip mismatch at n = %d\n", n);
                ok = false;
            }
        }
        report(ok, "n = " + std::to_string(n) + ": " + std::to_string(trials) +
                       " random round trips (both implementations)");
    }
}

void usage(const char* prog) {
    std::printf("usage: %s [--quick | --full | --help]\n", prog);
    std::printf("  (default) exhaustive checks up to n = 9, both implementations\n");
    std::printf("  --quick   smaller random sample\n");
    std::printf("  --full    additionally exhaustively check n = 10 (3,628,800 codes)\n");
}

}  // namespace

int main(int argc, char** argv) {
    bool full = false;
    int trials = 20000;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--full")       full = true;
        else if (a == "--quick") trials = 1000;
        else if (a == "--help" || a == "-h") { usage(argv[0]); return 0; }
        else { std::printf("unknown option: %s\n", a.c_str()); usage(argv[0]); return 2; }
    }

    std::printf("heaps_unrank verification suite\n");
    std::printf("checking: step_map / unrank / rank            (main,     Theta(n^2))\n");
    std::printf("          unrank_o_n3 / rank_o_n3             (baseline, Theta(n^3))\n");

    test_step_map_closed_form();
    test_step_dichotomy();

    section("3. exhaustive: every code of order n against the Heap sequence");
    for (int n = 1; n <= (full ? 10 : 9); n++) test_exhaustive_against_sequence(n);

    section("4. exhaustive round trip rank(unrank(C)) == C");
    for (int n = 1; n <= 8; n++) test_exhaustive_roundtrip(n);

    test_random_roundtrips(trials);

    std::printf("\n%s  %d checks, %d failures\n",
                g_failures == 0 ? "ALL CHECKS PASSED." : "FAILURES DETECTED.",
                g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
