// =============================================================================
//  heaps_unrank_variant.hpp -- Rank / Unrank for the swap-after-call variant
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.2 (2026-11-24)
//  License: MIT (see LICENSE)
//
//  Heap's recursion has a second common shape, in which the swap is performed
//  *after* each recursive call while the loop index runs over 0 .. k-1:
//
//      gen(k):
//          if k == 1: emit; return
//          for i = 0 .. k-1:
//              gen(k-1)
//              if k is odd: swap(A[0], A[k-1]) else swap(A[i], A[k-1])
//
//  That loop still emits every permutation of n elements exactly once, but the
//  resulting order is NOT Heap's order: the two sequences agree for n <= 3 and
//  differ from n = 4 on (the fourth block hands back a different arrangement).
//  This module supplies rank / unrank for that order.
//
//  The structure is the same as in heaps_unrank.hpp -- factoradic digits
//  d_k in [0, k) stored at C[n-k], one level per digit, Theta(n^2) time and
//  O(n) space -- but the per-level map is even simpler.  If block q of level k
//  is entered with arrangement A, it starts with W(k, q) . A, where
//
//      k odd  : W(k, q) = right rotation of the first k slots by q positions,
//      k even : W(k, q) = (0, k-1)(1, k-1) ... (q-1, k-1).
//
//  Both are given in closed form by base_map(), and both their inverses are
//  closed form as well, so neither direction simulates the q steps of a block.
//
//  Note: this order does not have Heap's minimal-change property -- at some
//  block boundaries consecutive outputs differ by more than one transposition.
// =============================================================================

#ifndef HEAPS_UNRANK_VARIANT_HPP
#define HEAPS_UNRANK_VARIANT_HPP

#include <vector>

namespace heaps {
namespace variant {

/// Scratch space reused across calls, so that the exhaustive tests can run
/// millions of unrank/rank calls without allocating.
struct Workspace {
    std::vector<int> p;       // working arrangement
    std::vector<int> map;     // W(k, q) of the current level
    std::vector<int> next;    // composition buffer
    std::vector<int> labels;  // working arrangement of rank()
    std::vector<int> pos;     // value -> slot
    std::vector<int> npos;    // value -> slot, next level
};

/// Source-index table W(k, q) of block q at level k: the block starts with
/// A_start[i] = A_enter[W[i]] for i = 0 .. k-1.  `map` receives k entries.
void base_map(int k, int q, std::vector<int>& map);

/// Unrank -- Theta(n^2) time, O(n) space.
///
/// `C` is the high-first Lehmer code of length n,
///
///     C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!,
///     0 <= C[i] <= n-i-1,   C[n-1] = 0,
///
/// and `perm` (room for n entries) receives the permutation of that index in
/// the variant's sequence.
void unrank(int n, const int* C, int* perm);

/// Rank -- Theta(n^2) time, O(n) space; inverse of unrank().  `perm` must be a
/// permutation of 0 .. n-1 and `C` must have room for n entries.
void rank(int n, const int* perm, int* C);

/// Same as above, reusing caller-provided scratch space.
void unrank(int n, const int* C, int* perm, Workspace& workspace);
void rank(int n, const int* perm, int* C, Workspace& workspace);

/// Convenience overloads operating on std::vector.
std::vector<int> unrank(int n, const std::vector<int>& C);
std::vector<int> rank(int n, const std::vector<int>& perm);

}  // namespace variant
}  // namespace heaps

#endif  // HEAPS_UNRANK_VARIANT_HPP
