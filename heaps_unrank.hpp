// =============================================================================
//  heaps_unrank.hpp -- Rank / Unrank for the Heap permutation generation order
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.0 (2026-10-15)
//  License: MIT (see LICENSE)
//
//  Reference implementation for the paper
//    "Ranking and Unranking Permutations in Heap's Order"
//  covering Algorithm 1 (step_map), Algorithm 2 (unrank), Algorithm 3 (rank).
// =============================================================================

#ifndef HEAPS_UNRANK_HPP
#define HEAPS_UNRANK_HPP

#include <vector>

namespace heaps {

/// Algorithm 1 -- in-block net mapping.
///
/// Returns the source-index table t of one complete in-block traversal at
/// order k, i.e. the operator Lambda_{F_k} of the paper: applying it to an
/// array X of length k means
///
///     X[i] <- X[t[i]]   for i = 0 .. k-1.
///
/// The table coincides with the paper's final arrangement F_k, so it can be
/// written down in O(k) time.
std::vector<int> step_map(int k);

/// Algorithm 2 -- unrank.
///
/// High-first Lehmer code C of length n  ->  Heap permutation perm of length n.
/// The code is interpreted as
///
///     C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!,
///     0 <= C[i] <= n-i-1,   C[n-1] = 0.
///
/// `perm` must have room for n entries.
void unrank(int n, const int* C, int* perm);

/// Algorithm 3 -- rank.
///
/// Heap permutation perm of length n  ->  high-first Lehmer code C of length n.
/// `perm` must be a permutation of 0 .. n-1; `C` must have room for n entries.
void rank(int n, const int* perm, int* C);

/// Convenience overloads operating on std::vector.
std::vector<int> unrank(int n, const std::vector<int>& C);
std::vector<int> rank(int n, const std::vector<int>& perm);

}  // namespace heaps

#endif  // HEAPS_UNRANK_HPP
