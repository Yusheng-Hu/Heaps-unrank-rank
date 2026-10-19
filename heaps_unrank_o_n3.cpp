// =============================================================================
//  heaps_unrank_o_n3_o_n3.cpp -- Baseline (dense, Theta(n^3)) Rank / Unrank_o_n3
//  for the Heap permutation generation order
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Baseline implementation for the paper (kept for reference and
//  cross-checking; the main implementation is heaps_unrank_o_n3.cpp)
//    "Ranking and Unrank_o_n3ing Permutations in Heap's Order"
//  Algorithms 2 (unrank_o_n3_o_n3) and 3 (rank_o_n3_o_n3); step_map() is defined
//  in heaps_unrank_o_n3.cpp.
//
//  Code convention (high-first Lehmer code C):
//      C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1!,   C[n-1] = 0,
//      0 <= C[i] <= n-i-1.
//  The paper uses the low-order-first code d, related by d[i] = C[n-2-i];
//  in particular the digit read at level k is C[n-k] = d[k-2].
//
//  Self-similarity: the block structure of Heap(k) is the accumulation of the
//  same in-block rearrangement rule (step_map) level after level.  The
//  trailing-element sequence is not independent information: it is exactly the
//  sequence of old values overwritten while the mapping S is advanced, so it
//  can be derived on the fly and never has to be precomputed (O(n) space, no
//  factorial-size lookup table).
// =============================================================================

#include "heaps_unrank_o_n3.hpp"
#include "heaps_unrank.hpp"   // step_map()

#include <algorithm>
#include <cassert>

namespace heaps {

// -----------------------------------------------------------------------------
// Algorithm 2: unrank.
//
// Bottom-up construction.  Starting from the seed permutation [0] at k = 1,
// the permutation is extended level by level up to k = n.  At level k the
// digit b = C[n-k] is the index of the block that is entered; the position
// mapping S is advanced from the base of block 0 to the base of block b by
// repeating the pair (in-block traversal, boundary interchange).
// -----------------------------------------------------------------------------
void unrank_o_n3(int n, const int* C, int* perm) {
    assert(n >= 1);
    assert(C != nullptr && perm != nullptr);
    perm[0] = 0;                        // level k = 1: the seed permutation

    std::vector<int> S, tmp, result;
    for (int k = 2; k <= n; k++) {
        const int size = k - 1;         // length of the permutation built so far
        const int b = C[n - k];         // block index at this level

        // level-invariant: the step table and the scratch buffers are reused
        const std::vector<int> t = step_map(size);

        S.resize(size);                 // S: position -> value, identity initially
        for (int i = 0; i < size; i++) S[i] = i;
        tmp.resize(size);
        int last = k - 1;               // trailing value at the base of block 0

        for (int j = 0; j < b; j++) {
            for (int i = 0; i < size; i++) tmp[i] = S[t[i]];
            S.swap(tmp);                // tmp now holds the previous mapping
            const int pos = (k % 2 == 0) ? j : 0;
            std::swap(S[pos], last);    // the value overwritten here becomes the
                                        // new trailing value (derived on the fly)
        }

        // translate the symbolic permutation of the previous level into values
        result.resize(size + 1);
        for (int i = 0; i < size; i++) result[i] = S[perm[i]];
        result[size] = last;
        for (int i = 0; i <= size; i++) perm[i] = result[i];
    }
}

// -----------------------------------------------------------------------------
// Algorithm 3: rank.
//
// Top-down reduction.  At level k, the trailing sequence is derived step by
// step until the trailing entry of the current permutation matches; the number
// of steps is the digit C[n-k].  The trailing element is then stripped off with
// S^{-1} and the algorithm descends one level.
// -----------------------------------------------------------------------------
void rank_o_n3(int n, const int* perm, int* C) {
    assert(n >= 1);
    assert(perm != nullptr && C != nullptr);
    for (int i = 0; i < n; i++) C[i] = 0;   // C[n-1] is always 0

    std::vector<int> p(perm, perm + n);     // working copy, mutated below

    assert(([&] {                          // the input must be a permutation of 0..n-1
        std::vector<char> seen(n, 0);
        for (int i = 0; i < n; i++) {
            if (p[i] < 0 || p[i] >= n || seen[p[i]]) return false;
            seen[p[i]] = 1;
        }
        return true;
    })());

    std::vector<int> S, tmp;
    for (int k = n; k >= 2; k--) {
        const int size = k - 1;
        const std::vector<int> t = step_map(size);

        S.resize(size);
        for (int i = 0; i < size; i++) S[i] = i;
        tmp.resize(size);
        int last = k - 1;                    // trailing value at the base of block 0

        // derive the trailing sequence until it meets the current trailing entry
        int b = 0;
        while (p[k - 1] != last && b < size) {  // "b < size" is only a safety net:
                                                // for valid input the loop is
                                                // bounded by size (paper, Lemma 5.1)
            for (int i = 0; i < size; i++) tmp[i] = S[t[i]];
            S.swap(tmp);
            const int pos = (k % 2 == 0) ? b : 0;
            std::swap(S[pos], last);
            b++;
        }
        C[n - k] = b;

        // inverse mapping: look each value up in S (S is a bijection, so S^{-1}
        // recovers the lower-order sub-permutation)
        for (int i = 0; i < size; i++) {
            int idx = 0;
            while (S[idx] != p[i]) idx++;
            p[i] = idx;
        }
    }
}

std::vector<int> unrank_o_n3(int n, const std::vector<int>& C) {
    std::vector<int> perm(n);
    unrank_o_n3(n, C.data(), perm.data());
    return perm;
}

std::vector<int> rank_o_n3(int n, const std::vector<int>& perm) {
    std::vector<int> C(n);
    rank_o_n3(n, perm.data(), C.data());
    return C;
}

}  // namespace heaps
