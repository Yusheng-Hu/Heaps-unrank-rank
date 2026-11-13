// =============================================================================
//  heaps_generator.cpp -- resumable generation for the Heap order (Section 5)
//
//  See heaps_generator.hpp for the two routes and for the conventions.
// =============================================================================

#include "heaps_generator.hpp"

#include "heaps_unrank.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace heaps {

namespace {

// n! as unsigned long long; throws when it does not fit.
unsigned long long factorial_checked(int n) {
    if (n < 0) {
        throw std::invalid_argument("heaps: negative order");
    }
    if (n > 20) {
        throw std::invalid_argument("heaps: order above 20 (n! does not fit)");
    }
    unsigned long long f = 1;
    for (int i = 2; i <= n; ++i) {
        f *= static_cast<unsigned long long>(i);
    }
    return f;
}

}  // namespace

std::vector<int> code_from_index(int n, unsigned long long index) {
    if (n < 0) {
        throw std::invalid_argument("heaps: negative order");
    }
    const unsigned long long total = factorial_checked(n);
    if (index >= total) {
        throw std::out_of_range("heaps: index out of range");
    }
    std::vector<int> code(static_cast<std::size_t>(n), 0);
    if (n <= 1) {
        return code;
    }
    // weight (n-1-i)! for position i; build the weights from the top
    std::vector<unsigned long long> w(static_cast<std::size_t>(n), 1);
    for (int i = n - 2; i >= 0; --i) {
        w[static_cast<std::size_t>(i)] = w[static_cast<std::size_t>(i + 1)] *
                                         static_cast<unsigned long long>(n - 1 - i);
    }
    unsigned long long rest = index;
    for (int i = 0; i < n - 1; ++i) {
        const unsigned long long wi = w[static_cast<std::size_t>(i)];
        code[static_cast<std::size_t>(i)] = static_cast<int>(rest / wi);
        rest %= wi;
    }
    code[static_cast<std::size_t>(n - 1)] = 0;
    return code;
}

std::vector<int> restore_context(int n, unsigned long long index) {
    const std::vector<int> code = code_from_index(n, index);
    std::vector<int> c(static_cast<std::size_t>(n), 0);
    // c[k] is the digit of weight k!, i.e. code[n-1-k]
    for (int k = 1; k <= n - 1; ++k) {
        c[static_cast<std::size_t>(k)] = code[static_cast<std::size_t>(n - 1 - k)];
    }
    return c;
}

FactorialHeapGenerator::FactorialHeapGenerator(int n, unsigned long long index)
    : n_(n), total_(factorial_checked(n)), index_(index),
      perm_(static_cast<std::size_t>(n)), code_(static_cast<std::size_t>(n)) {
    if (index >= total_) {
        throw std::out_of_range("heaps: index out of range");
    }
    code_ = code_from_index(n, index);
    perm_ = unrank(n, code_);
}

bool FactorialHeapGenerator::step() {
    if (index_ + 1 >= total_) {
        return false;
    }
    // (1) mixed-radix carry: the digit of level k (weight k!) lives at position
    // n-1-k of the high-first code, and it carries when it equals k.
    int k = 1;
    int pos = n_ - 1 - k;
    while (k <= n_ - 1 && code_[static_cast<std::size_t>(pos)] == k) {
        code_[static_cast<std::size_t>(pos)] = 0;
        ++k;
        pos = n_ - 1 - k;
    }
    if (k > n_ - 1) {
        return false;  // sequence exhausted (defensive; index_ check covers it)
    }
    const int old = code_[static_cast<std::size_t>(pos)];
    code_[static_cast<std::size_t>(pos)] = old + 1;
    // (2) transposition selected by the carry position
    if (k % 2 == 0) {
        std::swap(perm_[0], perm_[static_cast<std::size_t>(k)]);
    } else {
        std::swap(perm_[static_cast<std::size_t>(old)], perm_[static_cast<std::size_t>(k)]);
    }
    ++index_;
    return true;
}

}  // namespace heaps
