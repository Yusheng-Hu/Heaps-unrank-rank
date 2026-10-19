// =============================================================================
//  heaps_unrank.cpp -- Rank / Unrank for the Heap permutation order
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.1 (2026-10-19)
//  License: MIT (see LICENSE)
//
//  Companion implementation for the paper
//    "Ranking and Unranking Permutations in Heap's Order"
//
//  Difference from the reference implementation (heaps_unrank.cpp): the only
//  change is *how one step is applied*.  The reference version materialises
//  the step table and rearranges the whole mapping entry by entry, which costs
//  O(m) per step; here the step is applied directly:
//
//    * m odd  : the step is the transposition (0, m-1)  -> two swaps;
//    * m even : the step is a single m-cycle            -> advance an offset
//               along the cycle order and swap one entry.
//
//  The cycle order for even m is
//      cyc = [0, m-3, m-4, ..., 2, 1, m-2, m-1],
//  and cycle_index()/cycle_pos() map between a position and its index in cyc in
//  closed form, so no auxiliary table is needed.  Everything else -- the code
//  convention, the block/trailing-value bookkeeping, the level order -- is
//  identical to the reference implementation.
// =============================================================================

#include "heaps_unrank.hpp"

#include <algorithm>
#include <cassert>

namespace heaps {
namespace {

// Position p  ->  index of p in the cycle order cyc (even order m only).
//
// cyc = [0, m-3, m-4, ..., 2, 1, m-2, m-1], so
//     p = 0            -> index 0
//     1 <= p <= m-3    -> index m-2-p        (the descending run)
//     p = m-2, m-1     -> index p            (the tail, already in place)
int cycle_index(int m, int p) {
    if (p == 0) return 0;
    if (p <= m - 3) return m - 2 - p;
    return p;
}

// Index i in the cycle order  ->  position (inverse of cycle_index).
int cycle_pos(int m, int i) {
    if (i == 0) return 0;
    if (i <= m - 3) return m - 2 - i;
    return i;
}

}  // namespace

// -----------------------------------------------------------------------------
// unrank
//
// Bottom-up, level k = 2 .. n.  At level k the digit b = C[n-k] is the block
// index; the working mapping S is advanced b times, one O(1) step each.  Odd
// order stores S as a position -> value array; even order stores it in cycle
// order (W) together with an offset, so that the rotation of the m-cycle is a
// single integer increment.
// -----------------------------------------------------------------------------
void unrank(int n, const int* C, int* perm) {
    assert(n >= 1);
    assert(C != nullptr && perm != nullptr);
    perm[0] = 0;                          // level k = 1: the seed permutation

    std::vector<int> p(1, 0), S, W, nxt;
    for (int k = 2; k <= n; k++) {
        const int m = k - 1;              // order of the step taken at level k
        const int b = C[n - k];           // block index at this level
        int last = k - 1;                 // trailing value at the base of block 0
        nxt.resize((size_t)k);

        if (m % 2 != 0) {
            // odd order: one step == swap positions 0 and m-1, in O(1)
            S.resize((size_t)m);
            for (int i = 0; i < m; i++) S[(size_t)i] = i;
            for (int j = 0; j < b; j++) {
                std::swap(S[0], S[(size_t)(m - 1)]);
                std::swap(S[(size_t)j], last);   // boundary interchange
            }
            for (int i = 0; i < m; i++) nxt[(size_t)i] = S[(size_t)p[(size_t)i]];
            nxt[(size_t)m] = last;
        } else {
            // even order: one step == advance the offset along the m-cycle
            W.resize((size_t)m);
            for (int i = 0; i < m; i++) W[(size_t)i] = cycle_pos(m, i);
            int off = 0;                        // current rotation of the cycle
            for (int j = 0; j < b; j++) {
                off = (off + 1 == m) ? 0 : off + 1;
                std::swap(W[(size_t)off], last);   // boundary interchange
            }
            for (int i = 0; i < m; i++) {
                int ci = cycle_index(m, p[(size_t)i]) + off;
                if (ci >= m) ci -= m;
                nxt[(size_t)i] = W[(size_t)ci];
            }
            nxt[(size_t)m] = last;
        }
        p.swap(nxt);
    }
    std::copy(p.begin(), p.end(), perm);
}

// -----------------------------------------------------------------------------
// rank
//
// Top-down, level k = n .. 2: advance the mapping until the trailing entry of
// the current permutation matches (the number of advances is the digit b), then
// strip the trailing element with the inverse mapping.  The inverse is built in
// one O(m) pass instead of being looked up by a linear scan per entry.
// -----------------------------------------------------------------------------
void rank(int n, const int* perm, int* C) {
    assert(n >= 1);
    assert(perm != nullptr && C != nullptr);
    for (int i = 0; i < n; i++) C[i] = 0;     // C[n-1] is always 0

    std::vector<int> p(perm, perm + n);       // working copy
    std::vector<int> S, W, inv, nxt;

    assert(([&] {                             // input must be a permutation
        std::vector<char> seen((size_t)n, 0);
        for (int i = 0; i < n; i++) {
            if (p[(size_t)i] < 0 || p[(size_t)i] >= n || seen[(size_t)p[(size_t)i]]) return false;
            seen[(size_t)p[(size_t)i]] = 1;
        }
        return true;
    })());

    for (int k = n; k >= 2; k--) {
        const int m = k - 1;
        const int target = p[(size_t)(k - 1)];   // trailing entry to reproduce
        int last = k - 1;
        int b = 0;

        if (m % 2 != 0) {
            S.resize((size_t)m);
            for (int i = 0; i < m; i++) S[(size_t)i] = i;
            while (target != last && b < m) {    // "b < m" is a safety net only
                std::swap(S[0], S[(size_t)(m - 1)]);
                std::swap(S[(size_t)b], last);
                b++;
            }
            inv.assign((size_t)m + 1, -1);       // S may hold the value m
            for (int i = 0; i < m; i++) inv[(size_t)S[(size_t)i]] = i;
            nxt.resize((size_t)m);
            for (int i = 0; i < m; i++) nxt[(size_t)i] = inv[(size_t)p[(size_t)i]];
        } else {
            W.resize((size_t)m);
            for (int i = 0; i < m; i++) W[(size_t)i] = cycle_pos(m, i);
            int off = 0;
            while (target != last && b < m) {
                off = (off + 1 == m) ? 0 : off + 1;
                std::swap(W[(size_t)off], last);
                b++;
            }
            inv.assign((size_t)m + 1, -1);
            for (int i = 0; i < m; i++) {
                int ci = i + off;
                if (ci >= m) ci -= m;
                inv[(size_t)W[(size_t)ci]] = cycle_pos(m, i);
            }
            nxt.resize((size_t)m);
            for (int i = 0; i < m; i++) nxt[(size_t)i] = inv[(size_t)p[(size_t)i]];
        }
        C[n - k] = b;
        p.swap(nxt);
    }
}

std::vector<int> unrank(int n, const std::vector<int>& C) {
    std::vector<int> perm((size_t)n);
    unrank(n, C.data(), perm.data());
    return perm;
}

std::vector<int> rank(int n, const std::vector<int>& perm) {
    std::vector<int> C((size_t)n);
    rank(n, perm.data(), C.data());
    return C;
}

// -----------------------------------------------------------------------------
// Algorithm 1: step_map(k)
//
// Source-index table of the net position mapping of one complete in-block
// traversal at order k (the paper's Lambda_{F_k}).  Closed form:
//     k=1 -> [0]                      (identity)
//     k=2 -> [1, 0]
//     k=3 -> [2, 1, 0]
//     k=4 -> [1, 2, 3, 0]
//     k=5 -> [4, 1, 2, 3, 0]
//     k=6 -> [3, 4, 1, 2, 5, 0]
// -----------------------------------------------------------------------------
std::vector<int> step_map(int k) {
    assert(k >= 1);
    if (k == 1) return std::vector<int>{0};
    if (k == 2) return std::vector<int>{1, 0};  // needed separately: the even
                                                // formula below would give t[0] = -1
    std::vector<int> t(k);
    if (k % 2 != 0) {               // odd : t = [k-1, 1, 2, ..., k-2, 0]
        t[0] = k - 1;
        for (int i = 1; i < k - 1; i++) t[i] = i;
        t[k - 1] = 0;
    } else {                        // even: t = [k-3, k-2, 1, ..., k-4, k-1, 0]
        t[0] = k - 3;
        t[1] = k - 2;
        for (int i = 2; i < k - 2; i++) t[i] = i - 1;
        t[k - 2] = k - 1;
        t[k - 1] = 0;
    }
    return t;
}

}  // namespace heaps
