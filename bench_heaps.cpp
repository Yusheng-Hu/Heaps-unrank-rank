// =============================================================================
//  bench_heaps.cpp -- wall-clock comparison of the two implementations
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Times one unrank call of the baseline (dense, Theta(n^3)) implementation
//  against the main (Theta(n^2)) one for a range of orders, using the same
//  random code for both so the comparison is like with like.
//
//  This program measures wall-clock time only.  The operation counts reported
//  in the paper come from an instrumented build; see the paper, Section 6.
//
//  The ratio grows with n, in the direction implied by the O(k) factor that the
//  main implementation saves at every level; for small n it is dominated by
//  constants and cache effects, which is why the operation counts are the more
//  reliable indicator of the asymptotic order.
// =============================================================================

#include "heaps_unrank.hpp"
#include "heaps_unrank_o_n3.hpp"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <typename F>
double best_of(int rounds, F f) {
    double best = 1e100;
    for (int r = 0; r < rounds; r++) {
        const Clock::time_point t0 = Clock::now();
        f();
        const Clock::time_point t1 = Clock::now();
        const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        if (ns < best) best = ns;
    }
    return best;
}

}  // namespace

int main() {
    std::printf("heaps_unrank benchmark -- single unrank call, best of 5\n");
    std::printf("%6s %16s %16s %9s\n", "n", "baseline (ns)", "main (ns)", "ratio");

    std::mt19937 rng(20261019u);
    const int sizes[] = {10, 20, 40, 60, 80, 100, 150, 200};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
        const int n = sizes[s];
        std::vector<int> C((size_t)n), perm((size_t)n);
        for (int j = 0; j <= n - 2; j++) C[(size_t)j] = (int)(rng() % (unsigned)(n - j));
        C[(size_t)(n - 1)] = 0;
        heaps::unrank_o_n3(n, C.data(), perm.data());   // warm up both
        heaps::unrank(n, C.data(), perm.data());
        const double tb = best_of(5, [&] { heaps::unrank_o_n3(n, C.data(), perm.data()); });
        const double tm = best_of(5, [&] { heaps::unrank(n, C.data(), perm.data()); });
        std::printf("%6d %16.1f %16.1f %8.2fx\n", n, tb, tm, tb / tm);
    }
    return 0;
}
