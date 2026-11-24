// =============================================================================
//  test_variant.cpp -- verification of rank / unrank for the swap-after-call
//                      variant of Heap's recursion (heaps_unrank_variant.hpp)
//
//  Checks:
//
//   1. base_map(k, q) reproduces the arrangement on which block q of level k
//      actually starts, for k = 2..8 (ground truth: an instrumented run of the
//      variant generator).
//   2. Exhaustive, in both directions, for n = 1..9 (409,110 permutations by
//      default): the factoradic digits of the emission index must unrank to the
//      emitted permutation, and ranking that permutation must give the digits
//      back.  The digit counter is independent of both implementations.
//   3. For n <= 3 the variant order coincides with Heap's order, so the two
//      implementations must agree there.
//   4. Round trips at n = 64, 128, 256 on random codes.
//
//  Options: --quick  exhaustive up to n = 8  ( 46,232 permutations, ~0.05 s)
//           (none)   exhaustive up to n = 9  (409,110 permutations, ~0.4 s)
//           --full   exhaustive up to n = 11 (43,954,710 permutations, ~40 s)
//
//  Exit status: 0 when every check passes, 1 otherwise.
// =============================================================================

#include "heaps_unrank.hpp"
#include "heaps_unrank_variant.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) ++g_failures;
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
}

void section(const char* title) { std::printf("\n%s\n", title); }

// Ground truth: the variant generator, emitting every permutation of k elements
// in order.  `emit` is called once per output, with A holding that output.
template <class Emit>
void enumerate(int k, std::vector<int>& A, Emit& emit) {
    if (k == 1) {
        emit();
        return;
    }
    if (k % 2 != 0) {
        for (int i = 0; i < k; ++i) {
            enumerate(k - 1, A, emit);
            std::swap(A[0], A[static_cast<std::size_t>(k - 1)]);
        }
    } else {
        for (int i = 0; i < k; ++i) {
            enumerate(k - 1, A, emit);
            std::swap(A[static_cast<std::size_t>(i)], A[static_cast<std::size_t>(k - 1)]);
        }
    }
}

bool is_permutation(const std::vector<int>& p) {
    std::vector<char> seen(p.size(), 0);
    for (int v : p) {
        if (v < 0 || static_cast<std::size_t>(v) >= p.size()) return false;
        if (seen[static_cast<std::size_t>(v)]) return false;
        seen[static_cast<std::size_t>(v)] = 1;
    }
    return true;
}

long long factorial(int n) {
    long long f = 1;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

// --- 1 ----------------------------------------------------------------------
void check_base_map(int k) {
    const long long block = factorial(k - 1);
    std::vector<int> A(static_cast<std::size_t>(k));
    std::iota(A.begin(), A.end(), 0);

    long long index = 0;
    std::vector<int> map;
    bool ok = true;
    auto emit = [&]() {
        if (ok && index % block == 0) {
            const int q = static_cast<int>(index / block);
            heaps::variant::base_map(k, q, map);
            for (int i = 0; i < k; ++i) {
                if (A[static_cast<std::size_t>(i)] != map[static_cast<std::size_t>(i)]) {
                    ok = false;
                    std::printf("    k=%d q=%d slot %d: generator has %d, base_map %d\n", k,
                                q, i, A[static_cast<std::size_t>(i)],
                                map[static_cast<std::size_t>(i)]);
                    break;
                }
            }
        }
        ++index;
    };
    enumerate(k, A, emit);
    check(ok, "base_map reproduces the block starts at level k = " + std::to_string(k));
}

// --- 2 ----------------------------------------------------------------------
// Exhaustive in both directions, streaming (nothing is stored, so n = 11 fits).
void check_exhaustive(int n) {
    std::vector<int> A(static_cast<std::size_t>(n));
    std::iota(A.begin(), A.end(), 0);

    std::vector<int> C(static_cast<std::size_t>(n), 0);      // independent digit counter
    std::vector<int> got(static_cast<std::size_t>(n));
    std::vector<int> back(static_cast<std::size_t>(n));
    heaps::variant::Workspace workspace;

    long long count = 0;
    long long bad = 0;
    const auto start = std::chrono::steady_clock::now();

    auto emit = [&]() {
        heaps::variant::unrank(n, C.data(), got.data(), workspace);
        if (got != A) {
            ++bad;
            if (bad <= 3) {
                std::printf("    n=%lld unrank mismatch at index %lld\n",
                            static_cast<long long>(n), count);
            }
        }
        heaps::variant::rank(n, A.data(), back.data(), workspace);
        if (back != C) {
            ++bad;
            if (bad <= 3) {
                std::printf("    n=%lld rank mismatch at index %lld\n",
                            static_cast<long long>(n), count);
            }
        }
        ++count;
        // mixed-radix increment: d_k has radix k, most significant digit first
        for (int k = 2; k <= n; ++k) {
            if (++C[static_cast<std::size_t>(n - k)] < k) break;
            C[static_cast<std::size_t>(n - k)] = 0;
        }
    };
    enumerate(n, A, emit);

    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const long long expected = factorial(n);
    const bool ok = (bad == 0) && (count == expected);
    if (count != expected) {
        std::printf("    n=%d emitted %lld outputs, expected %lld\n", n, count, expected);
    }
    check(ok, "n = " + std::to_string(n) + ": " + std::to_string(count) +
                  " permutations, " + std::to_string(bad) + " mismatches (" +
                  std::to_string(seconds).substr(0, 5) + " s)");
}

// --- 3 ----------------------------------------------------------------------
void check_agrees_with_heap_order(int n) {
    std::vector<int> A(static_cast<std::size_t>(n));
    std::iota(A.begin(), A.end(), 0);

    bool ok = true;
    long long index = 0;
    std::vector<int> C(static_cast<std::size_t>(n), 0);
    std::vector<int> mine, heap;
    auto emit = [&]() {
        if (ok) {
            mine = heaps::variant::unrank(n, C);
            heap = heaps::unrank(n, C);
            if (mine != heap) {
                ok = false;
                std::printf("    n=%d index %lld: variant and Heap orders differ\n", n, index);
            }
        }
        ++index;
        for (int k = 2; k <= n; ++k) {
            if (++C[static_cast<std::size_t>(n - k)] < k) break;
            C[static_cast<std::size_t>(n - k)] = 0;
        }
    };
    enumerate(n, A, emit);
    check(ok, "n = " + std::to_string(n) + ": variant order equals Heap's order");
}

// --- 4 ----------------------------------------------------------------------
void check_round_trips() {
    std::mt19937 rng(20261124u);
    for (int n : {64, 128, 256}) {
        bool ok = true;
        for (int trial = 0; trial < 20 && ok; ++trial) {
            std::vector<int> C(static_cast<std::size_t>(n), 0);
            for (int k = 2; k <= n; ++k) {
                C[static_cast<std::size_t>(n - k)] =
                    static_cast<int>(rng() % static_cast<unsigned>(k));
            }
            const std::vector<int> perm = heaps::variant::unrank(n, C);
            if (!is_permutation(perm)) {
                ok = false;
                std::printf("    n=%d: unrank returned a non-permutation\n", n);
                break;
            }
            if (heaps::variant::rank(n, perm) != C) {
                ok = false;
                std::printf("    n=%d: round trip mismatch\n", n);
            }
        }
        check(ok, "round trips at n = " + std::to_string(n));
    }
}

void usage(const char* prog) {
    std::printf("usage: %s [--quick | --full | --help]\n", prog);
    std::printf("  (default) exhaustive checks up to n = 9   (409,110 permutations)\n");
    std::printf("  --quick   exhaustive checks up to n = 8   ( 46,232 permutations)\n");
    std::printf("  --full    exhaustive checks up to n = 11  (43,954,710 permutations)\n");
}

}  // namespace

int main(int argc, char** argv) {
    int n_max = 9;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--quick") {
            n_max = 8;
        } else if (a == "--full") {
            n_max = 11;
        } else if (a == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            std::printf("unknown option: %s\n", a.c_str());
            usage(argv[0]);
            return 2;
        }
    }

    std::printf("heaps_unrank_variant verification suite\n");
    std::printf("checking: unrank_variant / rank_variant     (swap-after-call order, Theta(n^2))\n");
    std::printf("          exhaustive up to n = %d\n", n_max);

    section("Block starts (closed form W(k, q) vs. the generator)");
    for (int k = 2; k <= 8; ++k) check_base_map(k);

    section("Exhaustive check of unrank and rank");
    for (int n = 1; n <= n_max; ++n) check_exhaustive(n);

    section("Agreement with Heap's order (holds for n <= 3 only)");
    for (int n = 1; n <= 3; ++n) check_agrees_with_heap_order(n);

    section("Round trips at larger n");
    check_round_trips();

    std::printf("\n%s  %d checks, %d failures\n",
                g_failures == 0 ? "ALL CHECKS PASSED." : "FAILURES DETECTED.", g_checks,
                g_failures);
    return g_failures == 0 ? 0 : 1;
}
