// =============================================================================
//  heaps_unrank_variant.cpp -- Rank / Unrank for the swap-after-call variant
//
//  Author : YUSHENG HU <dr.huyusheng@gmail.com>
//  Version: 2.2 (2026-11-24)
//  License: MIT (see LICENSE)
//
//  The order.  The variant generator (see the header) emits the same n!
//  permutations as Heap's algorithm but in a different order.  Its block
//  structure is the same: at level k the block index b = C[n-k] selects one of
//  k blocks of size (k-1)!, and block q starts with W(k, q) applied to the
//  arrangement that entered the level.  The direction matters:
//
//      A_emit = W(n, d_n) . W(n-1, d_{n-1}) . ... . W(2, d_2) . A_enter,
//
//  so unrank has to start at k = n (outermost level first) and compose
//  p <- p . W(k, b) going down, while rank peels k = n, n-1, ..., 2 off the top.
//
//  Cost.  W(k, q) is closed form (base_map), and so is its inverse, so a level
//  costs O(k) in either direction and the whole algorithms cost Theta(n^2).
// =============================================================================

#include "heaps_unrank_variant.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace heaps {
namespace variant {

// -----------------------------------------------------------------------------
// W(k, q): source-index table of the arrangement on which block q of level k
// starts, relative to the arrangement that entered the level.
//
//   k odd  : right rotation by q, i.e. W[i] = (i - q) mod k;
//   k even : the transpositions (i, k-1) for i = 0 .. q-1, in that order.
// -----------------------------------------------------------------------------
void base_map(int k, int q, std::vector<int>& map) {
    assert(k >= 1);
    assert(q >= 0 && q < k);
    map.resize(static_cast<std::size_t>(k));
    if (k % 2 != 0) {
        for (int i = 0; i < k; ++i) {
            map[static_cast<std::size_t>(i)] = ((i - q) % k + k) % k;
        }
    } else {
        for (int i = 0; i < k; ++i) {
            map[static_cast<std::size_t>(i)] = i;
        }
        for (int i = 0; i < q; ++i) {
            std::swap(map[static_cast<std::size_t>(i)],
                      map[static_cast<std::size_t>(k - 1)]);
        }
    }
}

// -----------------------------------------------------------------------------
// unrank: level k = n down to 2, composing p <- p . W(k, b) each time.
// Composition: (p . W)[i] = p[W[i]], which is why the table maps destination
// slot to source slot.
// -----------------------------------------------------------------------------
void unrank(int n, const int* C, int* perm, Workspace& workspace) {
    assert(n >= 1);
    assert(C != nullptr && perm != nullptr);

    std::vector<int>& p = workspace.p;
    p.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        p[static_cast<std::size_t>(i)] = i;
    }

    for (int k = n; k >= 2; --k) {
        const int b = C[n - k];
        base_map(k, b, workspace.map);
        std::vector<int>& next = workspace.next;
        next.assign(p.begin(), p.end());
        for (int i = 0; i < k; ++i) {
            next[static_cast<std::size_t>(i)] =
                p[static_cast<std::size_t>(workspace.map[static_cast<std::size_t>(i)])];
        }
        p.swap(next);
    }

    for (int i = 0; i < n; ++i) {
        perm[i] = p[static_cast<std::size_t>(i)];
    }
}

void unrank(int n, const int* C, int* perm) {
    Workspace workspace;
    unrank(n, C, perm, workspace);
}

// -----------------------------------------------------------------------------
// rank: level k = n down to 2.
//
// The digit is read off the trailing entry t = A[k-1]: the arrangement of block
// b is the only one whose last slot holds that value, and
//
//   k odd  : W(k, b) sends entry k-1-b to slot k-1, so b = (k-1-t) mod k;
//   k even : W(k, 0) keeps k-1, and W(k, b) with b >= 1 sends entry b-1 there,
//            so b = 0 if t == k-1 and b = t+1 otherwise.
//
// Peeling one level replaces every value v < k by the pre-image under W(k, b),
// i.e. it permutes the *labels*, not the slots.  Doing that through pos[]
// (value -> slot) costs O(k) instead of an O(n) search per entry.
// -----------------------------------------------------------------------------
void rank(int n, const int* perm, int* C, Workspace& workspace) {
    assert(n >= 1);
    assert(perm != nullptr && C != nullptr);

    std::vector<int>& labels = workspace.labels;
    labels.assign(perm, perm + n);
    for (int i = 0; i < n; ++i) {
        C[i] = 0;                              // C[n-1] is always 0
    }

    std::vector<int>& pos = workspace.pos;
    std::vector<int>& npos = workspace.npos;
    pos.resize(static_cast<std::size_t>(n));
    npos.resize(static_cast<std::size_t>(n));
    for (int slot = 0; slot < n; ++slot) {
        pos[static_cast<std::size_t>(labels[static_cast<std::size_t>(slot)])] = slot;
    }

    for (int k = n; k >= 2; --k) {
        const int t = labels[static_cast<std::size_t>(k - 1)];
        const int b = (k % 2 != 0) ? ((k - 1 - t) % k + k) % k
                                   : (t == k - 1 ? 0 : t + 1);
        C[n - k] = b;

        for (int v = 0; v < k; ++v) {           // closed-form inverse of W(k, b)
            int u;
            if (k % 2 != 0) {
                u = (v + b) % k;                // inverse of a right rotation
            } else if (b == 0) {
                u = v;
            } else if (b == 1) {
                u = (v == 0) ? k - 1 : (v == k - 1 ? 0 : v);
            } else if (v == 0) {
                u = 1;
            } else if (v <= b - 2) {
                u = v + 1;
            } else if (v == b - 1) {
                u = k - 1;
            } else if (v == k - 1) {
                u = 0;
            } else {
                u = v;
            }
            npos[static_cast<std::size_t>(u)] = pos[static_cast<std::size_t>(v)];
        }

        for (int v = 0; v < k; ++v) {
            pos[static_cast<std::size_t>(v)] = npos[static_cast<std::size_t>(v)];
            labels[static_cast<std::size_t>(pos[static_cast<std::size_t>(v)])] = v;
        }
    }
}

void rank(int n, const int* perm, int* C) {
    Workspace workspace;
    rank(n, perm, C, workspace);
}

std::vector<int> unrank(int n, const std::vector<int>& C) {
    std::vector<int> perm(static_cast<std::size_t>(n));
    unrank(n, C.data(), perm.data());
    return perm;
}

std::vector<int> rank(int n, const std::vector<int>& perm) {
    std::vector<int> C(static_cast<std::size_t>(n));
    rank(n, perm.data(), C.data());
    return C;
}

}  // namespace variant
}  // namespace heaps
