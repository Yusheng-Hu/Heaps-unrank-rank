// =============================================================================
//  heaps_unrank_o_n3_o_n3.hpp -- Baseline (dense, Theta(n^3)) Rank / Unrank_o_n3
//  for the Heap permutation generation order
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Baseline implementation for the paper (kept for reference and
//  cross-checking; the main implementation is heaps_unrank_o_n3.hpp)
//    "Ranking and Unrank_o_n3ing Permutations in Heap's Order"
//  covering Algorithm 2 (unrank_o_n3_o_n3) and Algorithm 3 (rank_o_n3_o_n3).
//  The closed-form helper step_map() is declared in heaps_unrank_o_n3.hpp.
// =============================================================================

#ifndef HEAPS_UNRANK_O_N3_HPP
#define HEAPS_UNRANK_O_N3_HPP

#include <vector>

namespace heaps {

/// Algorithm 2 -- unrank.
///
/// High-first Lehmer code C of length n  ->  Heap permutation perm of length n.
/// The code is interpreted as
///
///     C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!,
///     0 <= C[i] <= n-i-1,   C[n-1] = 0.
///
/// `perm` must have room for n entries.
void unrank_o_n3(int n, const int* C, int* perm);

/// Algorithm 3 -- rank.
///
/// Heap permutation perm of length n  ->  high-first Lehmer code C of length n.
/// `perm` must be a permutation of 0 .. n-1; `C` must have room for n entries.
void rank_o_n3(int n, const int* perm, int* C);

/// Convenience overloads operating on std::vector.
std::vector<int> unrank_o_n3(int n, const std::vector<int>& C);
std::vector<int> rank_o_n3(int n, const std::vector<int>& perm);

}  // namespace heaps

#endif  // HEAPS_UNRANK_O_N3_HPP
