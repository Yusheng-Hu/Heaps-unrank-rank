// =============================================================================
//  test_heaps_unrank.cpp -- verification suite for the heaps_unrank library
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.0 (2026-10-15)
//  License: MIT (see LICENSE)
//
//  The suite checks the implementation against independent ground truth:
//
//   1. step_map(k) equals the closed form F_k of the paper (Theorem 5.2),
//      for k = 1..24.
//   2. Exhaustive comparison with the iterative Heap algorithm (the ground
//      truth generator) for every code of order n = 1..9 -- 362,880 codes at
//      n = 9.  With --full, n = 10 (3,628,800 codes) is included as well.
//   3. Exhaustive round trip rank(unrank(C)) == C for every code of order
//      n = 1..8 (all n! codes).
//   4. Randomised round trips at larger orders (n up to 100), together with a
//      validity check that unrank really returns a permutation of 0..n-1.
//
//  Exit status is 0 when every check passes, 1 when a check fails, 2 on a
//  usage error.
//
//  Usage:  test_heaps_unrank [--quick | --full | --help]
//            (default)  checks 1-4 with n <= 9        (~15 s)
//            --quick    smaller random sample          (~2 s)
//            --full     additionally n = 10           (~25 s)
// =============================================================================

#include "heaps_unrank.hpp"

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
    explicit HeapStream(int n) : A_(n), c_(n, 0), n_(n), i_(0) {
        for (int i = 0; i < n; i++) A_[i] = i;
    }
    const std::vector<int>& current() const { return A_; }
    // Advance to the next permutation in Heap's order (no-op at the end).
    void advance() {
        if (i_ >= n_) return;
        for (;;) {
            if (c_[i_] < i_) {
                if (i_ % 2 == 0) std::swap(A_[0], A_[i_]);
                else             std::swap(A_[c_[i_]], A_[i_]);
                c_[i_]++;
                i_ = 0;
                return;
            }
            c_[i_] = 0;
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
    std::vector<int> f(m);
    if (m % 2 != 0) {                 // odd : F_m = [m-1, 1, 2, ..., m-2, 0]
        f[0] = m - 1;
        for (int i = 1; i < m - 1; i++) f[i] = i;
        f[m - 1] = 0;
    } else {                          // even: F_m = [m-3, m-2, 1, ..., m-4, m-1, 0]
        f[0] = m - 3;
        f[1] = m - 2;
        for (int i = 2; i < m - 2; i++) f[i] = i - 1;
        f[m - 2] = m - 1;
        f[m - 1] = 0;
    }
    return f;
}

// ---- high-first Lehmer code counter: the lowest digit moves first ----------
void next_code(std::vector<int>& C, int n) {
    if (n < 2) return;                       // n = 1 has a single code
    int j = n - 2;
    C[j]++;
    while (j > 0 && C[j] > n - j - 1) { C[j] = 0; j--; C[j]++; }
}

long long factorial(int n) {
    long long f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

bool is_permutation(const std::vector<int>& p, int n) {
    std::vector<char> seen((size_t)n, 0);
    for (int i = 0; i < n; i++) {
        if (p[i] < 0 || p[i] >= n || seen[(size_t)p[i]]) return false;
        seen[(size_t)p[i]] = 1;
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
// 2. exhaustive: every code of order n against the Heap sequence
// =============================================================================
void test_exhaustive_against_sequence(int n) {
    HeapStream heap(n);
    std::vector<int> C((size_t)n, 0), got((size_t)n);
    const long long total = factorial(n);
    long long r = 0;
    bool ok = true;
    for (; r < total; r++, next_code(C, n)) {
        heaps::unrank(n, C.data(), got.data());
        if (got != heap.current()) {
            std::printf("    mismatch at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
        heap.advance();
    }
    report(ok, "n = " + std::to_string(n) + ": " + std::to_string(r) +
                   " codes match the iterative Heap sequence");
}

// =============================================================================
// 3. exhaustive round trip rank(unrank(C)) == C
// =============================================================================
void test_exhaustive_roundtrip(int n) {
    std::vector<int> C((size_t)n, 0), back((size_t)n), perm((size_t)n);
    const long long total = factorial(n);
    long long r = 0;
    bool ok = true;
    for (; r < total; r++, next_code(C, n)) {
        heaps::unrank(n, C.data(), perm.data());
        heaps::rank(n, perm.data(), back.data());
        if (C != back) {
            std::printf("    round trip mismatch at n = %d, code index %lld\n", n, r);
            ok = false;
            break;
        }
    }
    report(ok, "n = " + std::to_string(n) + ": round trip for " + std::to_string(r) + " codes");
}

// =============================================================================
// 4. randomised round trips at larger orders
// =============================================================================
void test_random_roundtrips(int trials) {
    section("4. randomised round trips at larger orders (fixed seed 20261015)");
    std::mt19937 rng(20261015u);
    const int sizes[] = {12, 16, 20, 32, 48, 64, 100};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
        const int n = sizes[s];
        std::vector<int> C((size_t)n), back((size_t)n), perm((size_t)n);
        bool ok = true;
        for (int t = 0; t < trials && ok; t++) {
            for (int j = 0; j <= n - 2; j++) C[(size_t)j] = (int)(rng() % (unsigned)(n - j));
            C[(size_t)(n - 1)] = 0;
            heaps::unrank(n, C.data(), perm.data());
            if (!is_permutation(perm, n)) {
                std::printf("    unrank returned a non-permutation at n = %d\n", n);
                ok = false;
                break;
            }
            heaps::rank(n, perm.data(), back.data());
            if (C != back) {
                std::printf("    round trip mismatch at n = %d\n", n);
                ok = false;
            }
        }
        report(ok, "n = " + std::to_string(n) + ": " + std::to_string(trials) +
                       " random round trips, output is a valid permutation");
    }
}

void usage(const char* prog) {
    std::printf("usage: %s [--quick | --full | --help]\n", prog);
    std::printf("  (default) exhaustive checks up to n = 9\n");
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
    std::printf("checking: step_map / unrank / rank\n");

    test_step_map_closed_form();

    section("2. exhaustive: every code of order n against the Heap sequence");
    for (int n = 1; n <= (full ? 10 : 9); n++) test_exhaustive_against_sequence(n);

    section("3. exhaustive round trip rank(unrank(C)) == C");
    for (int n = 1; n <= 8; n++) test_exhaustive_roundtrip(n);

    test_random_roundtrips(trials);

    std::printf("\n%s  %d checks, %d failures\n",
                g_failures == 0 ? "ALL CHECKS PASSED." : "FAILURES DETECTED.",
                g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
