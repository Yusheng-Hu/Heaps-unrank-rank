// =============================================================================
//  heaps_generator.hpp -- Heap's generation order as a resumable machine
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.2 (2026-10-20)
//  License: MIT (see LICENSE)
//
//  Companion implementation for the paper
//    "Ranking and Unranking Permutations in Heap's Order", Section 5.
//
//  The paper distinguishes two routes that both make the Heap sequence
//  resumable from an arbitrary index.  Their output sequences coincide; they
//  differ in where the change lives.
//
//  Route A -- keep the classical iterative loop and restore its counter array.
//      restore_context(n, index) returns the counter array c the classical
//      loop needs in order to continue with the permutation of index
//      `index + 1`, given that the arrays are read through heaps::unrank().
//      The classical counters are advanced *after* an emission, so at the
//      moment of an emission they describe the previous index; that off-by-one
//      convention is why the restore needs the code of `index`, not of
//      `index + 1`.  The scan position must be set to 0.
//
//  Route B -- drive the loop by the factorial counter instead.
//      FactorialHeapGenerator keeps the state (array, factorial code) and
//      advances the code by the standard mixed-radix carry rule; the digit
//      where the carry stops decides which transposition to apply.  Its state
//      is then the very object that heaps::unrank() / heaps::rank() consume,
//      and a restart is a copy rather than a reconstruction.
//
//  Both routes emit exactly the sequence of the classical Heap algorithm; the
//  test suite compares them against heaps::unrank() and against each other.
// =============================================================================

#ifndef HEAPS_GENERATOR_HPP
#define HEAPS_GENERATOR_HPP

#include <vector>

namespace heaps {

/// High-first Lehmer code of the given index at order n (length n, C[n-1] = 0),
/// i.e. the code whose value is index,
///
///     C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1!  ==  index.
///
/// The factorial n! must fit in unsigned long long (n <= 20).
std::vector<int> code_from_index(int n, unsigned long long index);

/// Route A: counter array for the classical iterative loop.
///
/// Returns c of length n with c[0] = 0 and, for k = 1 .. n-1,
/// c[k] = C[n-1-k] where C = code_from_index(n, index).  Restarting the
/// classical loop with the array heaps::unrank(n, C), the counters returned
/// here and the scan position 0 emits the permutation of index index + 1.
std::vector<int> restore_context(int n, unsigned long long index);

/// Route B: generator driven by the factorial counter.
class FactorialHeapGenerator {
public:
    /// Starts at the permutation of the given index.
    FactorialHeapGenerator(int n, unsigned long long index);

    /// Current permutation (Heap's sequence, index()).
    const std::vector<int>& current() const noexcept { return perm_; }

    /// Factorial code of the current permutation (same layout as unrank()).
    const std::vector<int>& code() const noexcept { return code_; }

    /// Index of the current permutation.
    unsigned long long index() const noexcept { return index_; }

    /// Number of permutations of this order.
    unsigned long long total() const noexcept { return total_; }

    /// Advances to the next permutation of the sequence.  Returns false (and
    /// leaves the state unchanged) when the end of the sequence is reached.
    bool step();

private:
    int n_;
    unsigned long long total_;
    unsigned long long index_;
    std::vector<int> perm_;
    std::vector<int> code_;
};

}  // namespace heaps

#endif  // HEAPS_GENERATOR_HPP
