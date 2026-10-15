# heaps_unrank — Rank / Unrank for the Heap permutation generation order

[![CI](https://github.com/USER/REPO/actions/workflows/ci.yml/badge.svg)](https://github.com/USER/REPO/actions/workflows/ci.yml)

Rank and unrank algorithms for **Heap's permutation generation order** (B. R. Heap,
*Permutations by interchanges*, 1963). Given an index, `unrank` returns the
corresponding permutation of Heap's sequence; given a permutation, `rank` returns
its index. A detailed index for the Heap order is obtained directly, with

* **O(n³) time**, **O(n) space**,
* **no precomputed table** of any kind (nothing of factorial size is stored),
* a single small closed-form helper (`step_map`) shared by both directions.

**Author:** YUSHENG HU <dr.huyusheng@gmail.com>
**Version:** 2.0
**License:** MIT — see [`LICENSE`](LICENSE)

Reference implementation for the paper *"Bidirectional Mapping Algorithms,
Operator-Algebraic Models, and Canonical Inductive Proofs for Heap's
Permutations"*; it corresponds to Algorithm 1 (`step_map`), Algorithm 2
(`unrank`) and Algorithm 3 (`rank`).

---

## Index convention

The index is a **high-first Lehmer code** `C` of length `n`:

```
index = C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!
0 <= C[i] <= n-i-1,        C[n-1] = 0
```

The paper writes the same index in low-order-first form `d`, related by
`d[i] = C[n-2-i]`; the digit read at level `k` is `C[n-k] = d[k-2]`.
The code space has `n!` elements, matching the `n!` permutations.

## What the code does

The `n!` permutation space is split into `n` blocks of `(n-1)!` permutations.
Inside a block the first `n-1` positions run through the lower-order Heap
permutations while the trailing element stays fixed, and at a block boundary a
single interchange updates the trailing element. Two facts make the algorithms
work:

* the net effect of one complete in-block traversal is a fixed permutation,
  available in closed form — this is `step_map(k)`, the paper's `Λ` operator;
* the trailing value of the next block is exactly the value overwritten by that
  interchange, so it is **derived on the fly** and never has to be tabulated.

## Files

| File | Purpose |
|---|---|
| `heaps_unrank.hpp` | public interface (Algorithm 1/2/3) |
| `heaps_unrank.cpp` | implementation |
| `test_heaps_unrank.cpp` | verification suite (ground-truth comparisons, exhaustive and randomised) |
| `Makefile` | simple build / test targets |
| `CMakeLists.txt` | CMake build with a CTest test |
| `.github/workflows/ci.yml` | GitHub Actions pipeline |

## Build and run

With `make`:

```bash
make            # build and run the standard verification suite  (~15 s)
make full       # full suite: also exhaustively checks n = 10    (~25 s)
make quick      # fast smoke test                                (~2 s)
```

With CMake:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Directly:

```bash
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -I. \
    heaps_unrank.cpp test_heaps_unrank.cpp -o test_heaps_unrank
./test_heaps_unrank --full
```

## Using the library

```cpp
#include "heaps_unrank.hpp"
#include <cstdio>

int main() {
    const int n = 4;
    int C[n] = {0, 0, 0, 0};        // index  (high-first Lehmer code)

    int perm[n];
    heaps::unrank(n, C, perm);      // index  -> permutation

    int back[n];
    heaps::rank(n, perm, back);     // permutation -> index

    for (int i = 0; i < n; i++) std::printf("%d ", perm[i]);
    std::printf("\n");
}
```

`std::vector` overloads are provided as well:

```cpp
std::vector<int> perm = heaps::unrank(n, C);        // C is a std::vector<int>
std::vector<int> C    = heaps::rank(n, perm);
```

## What the verification suite checks

| # | Check | Scale |
|---|---|---|
| 1 | `step_map(k)` equals the closed form `F_k` of the paper (Theorem 5.2) | `k = 1..24` |
| 2 | **Exhaustive**: every code is compared with the iterative Heap algorithm (ground truth) | `n = 1..9` — 362,880 codes at `n = 9`; `--full` adds `n = 10` — 3,628,800 codes |
| 3 | **Exhaustive** round trip `rank(unrank(C)) == C` | `n = 1..8` (all `n!` codes) |
| 4 | Randomised round trips plus a validity check that the output is a permutation of `0..n-1` | `n ∈ {12,16,20,32,48,64,100}`, 20,000 trials each |

The program exits with status 0 only when every check passes, and prints a
`PASS`/`FAIL` line per check. Check 2 is the strongest evidence: all `n!` codes
of order `n` are matched one by one against the permutation produced by Heap's
own iterative algorithm.

Compilation is warning-free under `-Wall -Wextra -Wpedantic`, and the suite runs
clean under AddressSanitizer and UndefinedBehaviorSanitizer.

## Complexity

| Direction | Time | Space | Precomputation |
|---|---|---|---|
| `unrank` | `Θ(n³)` | `O(n)` | none |
| `rank` | `Θ(n³)` | `O(n)` | none |

The bound comes from the `n-1` levels, each advancing at most `k-1` times at a
cost of `Θ(k)` per advance. Both directions keep only three linear arrays; no
factorial-size lookup table is used. See the paper, Section 7, for the full
analysis.

## GitHub Actions

The pipeline in `.github/workflows/ci.yml` runs three jobs:

| Job | What it does |
|---|---|
| `build & test (os)` | configures and builds with CMake, then runs the **full** suite through CTest on `ubuntu-latest`, `macos-latest` and `windows-latest` |
| `strict warnings (cxx)` | compiles library and tests with `-Wall -Wextra -Wpedantic -Werror` under `g++` and `clang++` |
| `ASan + UBSan` | rebuilds with `-fsanitize=address,undefined` and runs the suite |

**How it runs**

* **Automatically** on every push to `main`/`master` and on every pull request.
* **Manually**: open the repository on GitHub, go to the **Actions** tab, select
  **CI** in the left-hand list, click **Run workflow**, choose a branch and
  confirm. (This is the `workflow_dispatch` trigger.)
* **Results**: the same Actions tab shows each run; a red cross marks a failing
  check, and the failing step's log is available under the job. The test program
  prints one `PASS`/`FAIL` line per check, so the log identifies the failure
  directly.

**Badge**: replace `USER/REPO` in the badge at the top of this README with your
GitHub account and repository name.

**Reproducing a CI job locally**

```bash
# job "build & test"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure

# job "strict warnings"
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. -c heaps_unrank.cpp -o /tmp/hu.o
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. -c test_heaps_unrank.cpp -o /tmp/t.o

# job "ASan + UBSan"
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_unrank.cpp test_heaps_unrank.cpp -o test_san && ./test_san
```

## License

Released under the MIT License; see [`LICENSE`](LICENSE).

## References

1. B. R. Heap. *Permutations by interchanges.* The Computer Journal 6(3):293–294, 1963.
2. B. R. Heap's method appears in R. Sedgewick, *Permutation generation methods*,
   Computing Surveys 9:137–164, 1977, and in D. E. Knuth, *TAOCP Vol. 4, Fascicle 2*, 2005.
3. YUSHENG HU. *Bidirectional Mapping Algorithms, Operator-Algebraic Models, and
   Canonical Inductive Proofs for Heap's Permutations* (the companion paper).
