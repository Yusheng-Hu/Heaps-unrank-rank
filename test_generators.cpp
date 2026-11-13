// =============================================================================
//  test_generators.cpp -- consistency of the two resumable-generation routes
//
//  Checks (all against heaps::unrank() and against an independent classical
//  iteration of Heap's algorithm):
//
//   1. code_from_index(n, r) is the high-first Lehmer code of index r, and
//      heaps::unrank(n, code_from_index(n, r)) equals the r-th permutation of
//      the classical iterative Heap sequence, for every r at n = 1..7.
//   2. Route A: restoring the counters with restore_context(n, a) and running
//      one step of the classical loop yields exactly the permutation of index
//      a + 1, for every a at n = 2..6 (866 states in total).
//   3. Route B: FactorialHeapGenerator(n, 0) enumerated by step() reproduces the
//      same sequence, and its code() equals code_from_index(n, index()) at
//      every step, for n = 1..7.
//   4. The two routes agree digit by digit: restore_context(n, a)[k] ==
//      code_from_index(n, a)[n-1-k] for k = 1..n-1.
//   5. Random restarts at n = 8..12 land on the right permutation, and step()
//      reports the end of the sequence correctly.
//
//  Exit status: 0 when every check passes, 1 otherwise.
// =============================================================================

#include "heaps_generator.hpp"
#include "heaps_unrank.hpp"

#include <cstdio>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool ok, const char* what, int n, unsigned long long index) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        if (g_failures <= 10) {
            std::printf("  FAIL: %s  (n=%d, index=%llu)\n", what, n,
                        static_cast<unsigned long long>(index));
        }
    }
}

unsigned long long factorial_u(int n) {
    unsigned long long f = 1;
    for (int i = 2; i <= n; ++i) {
        f *= static_cast<unsigned long long>(i);
    }
    return f;
}

// Independent ground truth: the classical iterative Heap algorithm, identity
// emitted first.
std::vector<std::vector<int>> classical_sequence(int n) {
    std::vector<int> a(static_cast<std::size_t>(n));
    std::iota(a.begin(), a.end(), 0);
    std::vector<int> c(static_cast<std::size_t>(n), 0);
    std::vector<std::vector<int>> seq;
    seq.push_back(a);
    const unsigned long long total = factorial_u(n);
    int i = 0;
    while (static_cast<unsigned long long>(seq.size()) < total) {
        if (c[static_cast<std::size_t>(i)] < i) {
            if (i % 2 == 0) {
                std::swap(a[0], a[static_cast<std::size_t>(i)]);
            } else {
                std::swap(a[static_cast<std::size_t>(c[static_cast<std::size_t>(i)])],
                          a[static_cast<std::size_t>(i)]);
            }
            seq.push_back(a);
            ++c[static_cast<std::size_t>(i)];
            i = 0;
        } else {
            c[static_cast<std::size_t>(i)] = 0;
            ++i;
            if (i >= n) {
                break;
            }
        }
    }
    return seq;
}

void test_code_and_sequence() {
    for (int n = 1; n <= 8; ++n) {
        const auto seq = classical_sequence(n);
        for (unsigned long long r = 0; r < seq.size(); ++r) {
            const auto code = heaps::code_from_index(n, r);
            const auto perm = heaps::unrank(n, code);
            check(perm == seq[static_cast<std::size_t>(r)],
                  "unrank(code_from_index(r)) matches the classical sequence", n, r);
        }
    }
}

void test_route_a() {
    for (int n = 2; n <= 8; ++n) {
        const auto seq = classical_sequence(n);
        for (unsigned long long a = 0; a + 1 < seq.size(); ++a) {
            const auto code = heaps::code_from_index(n, a);
            std::vector<int> perm = seq[static_cast<std::size_t>(a)];
            std::vector<int> c = heaps::restore_context(n, a);
            int i = 0;
            bool emitted = false;
            for (int guard = 0; guard < 1000 && !emitted; ++guard) {
                if (c[static_cast<std::size_t>(i)] < i) {
                    if (i % 2 == 0) {
                        std::swap(perm[0], perm[static_cast<std::size_t>(i)]);
                    } else {
                        std::swap(perm[static_cast<std::size_t>(
                                      c[static_cast<std::size_t>(i)])],
                                  perm[static_cast<std::size_t>(i)]);
                    }
                    emitted = true;
                } else {
                    c[static_cast<std::size_t>(i)] = 0;
                    ++i;
                    if (i >= n) {
                        break;
                    }
                }
            }
            check(emitted && perm == seq[static_cast<std::size_t>(a + 1)],
                  "route A: one step after restore yields index a+1", n, a);
            // the two routes describe the same digits
            for (int k = 1; k <= n - 1; ++k) {
                check(heaps::restore_context(n, a)[static_cast<std::size_t>(k)] ==
                          code[static_cast<std::size_t>(n - 1 - k)],
                      "routes A and B agree digit by digit", n, a);
            }
        }
    }
}

void test_route_b() {
    for (int n = 1; n <= 8; ++n) {
        const auto seq = classical_sequence(n);
        heaps::FactorialHeapGenerator gen(n, 0);
        for (unsigned long long r = 0; r < seq.size(); ++r) {
            check(gen.index() == r, "route B: index advances by one", n, r);
            check(gen.current() == seq[static_cast<std::size_t>(r)],
                  "route B: current permutation matches", n, r);
            check(gen.code() == heaps::code_from_index(n, r),
                  "route B: code equals code_from_index", n, r);
            if (r + 1 < seq.size()) {
                check(gen.step(), "route B: step reports success", n, r);
            }
        }
        check(!gen.step(), "route B: step reports the end of the sequence", n,
              static_cast<unsigned long long>(seq.size()));
    }
}

void test_random_restarts() {
    std::mt19937_64 rng(20261020ULL);
    for (int n = 8; n <= 12; ++n) {
        const unsigned long long total = factorial_u(n);
        for (int t = 0; t < 50; ++t) {
            const unsigned long long r = rng() % total;
            heaps::FactorialHeapGenerator gen(n, r);
            const auto code = heaps::code_from_index(n, r);
            check(gen.current() == heaps::unrank(n, code),
                  "random restart lands on the right permutation", n, r);
            check(gen.code() == code, "random restart has the right code", n, r);
            if (r + 1 < total) {
                check(gen.step(), "random restart: step works", n, r);
                const auto next = heaps::code_from_index(n, r + 1);
                check(gen.current() == heaps::unrank(n, next),
                      "random restart: next permutation matches", n, r);
            }
        }
    }
}

void test_errors() {
    bool threw = false;
    try {
        heaps::FactorialHeapGenerator gen(21, 0);
        (void)gen;
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "order above 20 is rejected", 21, 0);

    threw = false;
    try {
        (void)heaps::restore_context(4, 24);
    } catch (const std::out_of_range&) {
        threw = true;
    }
    check(threw, "restore_context rejects an out-of-range index", 4, 24);

    threw = false;
    try {
        heaps::FactorialHeapGenerator gen(4, 24);
        (void)gen;
    } catch (const std::out_of_range&) {
        threw = true;
    }
    check(threw, "index out of range is rejected", 4, 24);
}

}  // namespace

int main() {
    std::printf("test_generators -- Heap order as a resumable machine\n\n");
    try {
        test_code_and_sequence();
        test_route_a();
        test_route_b();
        test_random_restarts();
        test_errors();
    } catch (const std::exception& e) {
        std::printf("  FAIL: unexpected exception: %s\n", e.what());
        return 1;
    }
    std::printf("\n  checks: %d   failures: %d\n", g_checks, g_failures);
    std::printf("%s\n", g_failures == 0 ? "  ALL CHECKS PASSED" : "  CHECKS FAILED");
    return g_failures == 0 ? 0 : 1;
}
