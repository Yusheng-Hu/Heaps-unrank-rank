// =============================================================================
//  heaps_unrank.hpp -- Rank / Unrank for the Heap permutation order
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Companion implementation for the paper
//    "Ranking and Unranking Permutations in Heap's Order"
//
//  It realises the per-step O(1) block transition of the paper (Section 4.5)
//  and the corollary following the cycle-dichotomy lemma: the source-index
//  table of one in-block step at order m -- which equals the paper's final
//  arrangement F_m, the same table that heaps::step_map() returns -- is
//     * the transposition (0, m-1)              when m is odd,
//     * a single m-cycle                        when m is even.
//  One step can therefore be applied with O(1) element moves instead of O(m),
//  which lowers the cost of each level from O(m^2) to O(m) and the whole
//  algorithms from Theta(n^3) to Theta(n^2).  Space stays O(n) and no table of
//  factorial size is used.
//
//  The reference (dense) implementation is heaps_unrank.hpp; it is kept as the
//  baseline, and the test suite verifies that both produce identical results.
// =============================================================================

#ifndef HEAPS_UNRANK_HPP
#define HEAPS_UNRANK_HPP

#include <vector>

namespace heaps {

/// Closed-form in-block step table (paper, Algorithm 1).
///
/// Returns the source-index table t of one complete in-block traversal at
/// order k, i.e. the operator Lambda_{F_k}: applying it to an array X of
/// length k means  X[i] <- X[t[i]]  for i = 0 .. k-1.  The table coincides
/// with the paper's final arrangement F_k.  The main implementation exploits
/// its structure directly and does not need the table; it is used by the
/// baseline implementation (heaps_unrank_o_n3.hpp) and by the test suite.
std::vector<int> step_map(int k);

/// Implicit unrank -- Theta(n^2) time, O(n) space.
///
/// Same convention as heaps::unrank(): `C` is the high-first Lehmer code of
/// length n,
///
///     C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!,
///     0 <= C[i] <= n-i-1,   C[n-1] = 0,
///
/// and `perm` (room for n entries) receives the corresponding permutation of
/// Heap's sequence.
void unrank(int n, const int* C, int* perm);

/// Implicit rank -- Theta(n^2) time, O(n) space.
///
/// Inverse of unrank(): `perm` must be a permutation of 0..n-1 and `C`
/// must have room for n entries.
void rank(int n, const int* perm, int* C);

/// Convenience overloads operating on std::vector.
std::vector<int> unrank(int n, const std::vector<int>& C);
std::vector<int> rank(int n, const std::vector<int>& perm);

}  // namespace heaps

#endif  // HEAPS_UNRANK_HPP
