# heaps_unrank -- Rank / Unrank for the Heap permutation generation order

[![CI](https://github.com/Yusheng-Hu/Heaps-unrank-rank/actions/workflows/ci.yml/badge.svg)](https://github.com/Yusheng-Hu/Heaps-unrank-rank/actions/workflows/ci.yml)
[![GitHub stars](https://img.shields.io/github/stars/Yusheng-Hu/Heaps-unrank-rank?style=social)](https://github.com/Yusheng-Hu/Heaps-unrank-rank)
![Field: Combinatorial Algorithms](https://img.shields.io/badge/Field-Combinatorial_Algorithms-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: C++](https://img.shields.io/badge/Language-C%2B%2B-red.svg)](#)
[![Complexity: O(n^2) / O(n^3)](https://img.shields.io/badge/Complexity-O(n%5E2)%20%2F%20O(n%5E3)-brightgreen.svg)](#)
[![Space: O(n)](https://img.shields.io/badge/Space-O(n)-green.svg)](#)

Rank and unrank algorithms for **Heap's permutation generation order** (B. R. Heap,
*Permutations by interchanges*, 1963). Given an index, `unrank` returns the
corresponding permutation of Heap's sequence; given a permutation, `rank` returns
its index.

Two implementations are provided side by side, sharing one code convention and
one test suite:

| Implementation | Header / API | Time | Space | Role |
|---|---|---|---|---|
| **Main** | `heaps_unrank.hpp` -- `unrank`, `rank` | **Theta(n^2)** | O(n) | the recommended one |
| **Baseline** (dense) | `heaps_unrank_o_n3.hpp` -- `unrank_o_n3`, `rank_o_n3` | Theta(n^3) | O(n) | the paper's dense version, kept as an independent cross-check |

Neither uses a precomputed table of any kind, and both build on the same single
closed-form helper (`step_map`, the paper's in-block step operator).

- **Repository:** <https://github.com/Yusheng-Hu/Heaps-unrank-rank>
- **Companion paper:** *Ranking and Unranking Permutations in Heap's Order* -- PDF
  in [`paper/heap_paper_en.pdf`](paper/heap_paper_en.pdf)
- **Author:** YUSHENG HU <dr.huyusheng@gmail.com>
- **Version:** 2.1
- **License:** MIT -- see [`LICENSE`](LICENSE)

The baseline implements Algorithms 2-3 of the paper as they are stated there,
applying the whole step table at every step. The main implementation keeps the
same construction but replaces the application of one step by the per-step O(1)
form of Section 4.5 and the corollary that follows the cycle-dichotomy lemma;
that is the only difference, and it is what lowers the total cost from
Theta(n^3) to Theta(n^2). The baseline files are deliberately left unchanged, so
the two implementations remain independent evidence for each other.

---

## Index convention

The index is a **high-first Lehmer code** `C` of length `n`:

```
index = C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!
0 <= C[i] <= n-i-1,        C[n-1] = 0
```

The paper writes the same index in low-order-first form `d`, related by
`d[i] = C[n-2-i]`; the digit read at level `k` is `C[n-k] = d[k-2]`.
The code space has `n!` elements, matching the `n!` permutations. Both
implementations take and return exactly the same convention, so they are
interchangeable in a program.

## What the code does

The `n!` permutation space is split into `n` blocks of `(n-1)!` permutations.
Inside a block the first `n-1` positions run through the lower-order Heap
permutations while the trailing element stays fixed, and at a block boundary a
single interchange updates the trailing element. Three facts make the algorithms
work:

* the net effect of one complete in-block traversal is a fixed permutation,
  available in closed form -- this is `step_map(k)`, the paper's `Lambda`
  operator, and it equals the final arrangement `F_k`;
* the trailing value of the next block is exactly the value overwritten by that
  interchange, so it is **derived on the fly** and never has to be tabulated;
* **`step_map(k)` has only two shapes** -- the transposition `(0, k-1)` for odd
  `k`, and a single `k`-cycle for even `k` -- so one step can be applied with
  O(1) element moves instead of O(k). This is what separates the main
  implementation from the baseline one, and it is what lowers the total cost
  from Theta(n^3) to Theta(n^2).

## Files

| File | Purpose |
|---|---|
| `heaps_unrank.hpp` | public interface of the **main** implementation (Theta(n^2)) + `step_map` |
| `heaps_unrank.cpp` | its implementation (per-step O(1) block transition) |
| `heaps_unrank_o_n3.hpp` | public interface of the **baseline** implementation (dense, Theta(n^3)) |
| `heaps_unrank_o_n3.cpp` | its implementation (Algorithms 2-3 applied densely) |
| `heaps_generator.hpp` | resumable generation (paper, Section 5): `restore_context`, `FactorialHeapGenerator` |
| `heaps_generator.cpp` | the two routes that start the Heap sequence at an arbitrary index |
| `test_generators.cpp` | consistency test of both routes against `heaps::unrank()` (n <= 8) |
| `test_heaps_unrank.cpp` | verification suite, exercises **both** implementations |
| `bench_heaps.cpp` | wall-clock comparison of the two implementations |
| `paper/heap_paper_en.pdf` | the companion paper (current version) |
| `Makefile` | simple build / test / bench targets |
| `CMakeLists.txt` | CMake build with two CTest tests |
| `.github/workflows/ci.yml` | GitHub Actions pipeline |

## Build and run

With `make`:

```bash
make            # build and run the standard suite + generator suite  (~20 s)
make full       # full suite: also exhaustively checks n = 10          (~40 s)
make quick      # quick smoke test                                     (~3 s)
make gen        # resumable-generation consistency test               (~0.3 s)
make bench      # wall-clock benchmark of both implementations
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
    heaps_unrank.cpp heaps_unrank_o_n3.cpp test_heaps_unrank.cpp -o test_heaps_unrank
./test_heaps_unrank --full
```

## Using the library

The two implementations are interchangeable and return identical results; use
the main one, and reach for the baseline only if you want the dense version for
reference or cross-checking.

```cpp
#include "heaps_unrank.hpp"        // main:     Theta(n^2)
#include "heaps_unrank_o_n3.hpp"   // baseline: Theta(n^3)
#include <cstdio>

int main() {
    const int n = 4;
    int C[n] = {0, 0, 0, 0};        // index (high-first Lehmer code)

    int main_perm[n], base_perm[n];
    heaps::unrank(n, C, main_perm);         // index -> permutation, Theta(n^2)
    heaps::unrank_o_n3(n, C, base_perm);    // index -> permutation, Theta(n^3)

    int back[n];
    heaps::rank(n, main_perm, back);        // permutation -> index, Theta(n^2)

    for (int i = 0; i < n; i++) std::printf("%d ", main_perm[i]);
    std::printf("\n");
    // main_perm and base_perm agree for every code, and rank(main_perm) == C
}
```

`std::vector` overloads are provided for both:

```cpp
std::vector<int> perm  = heaps::unrank(n, C);            // C is a std::vector<int>
std::vector<int> C2    = heaps::rank(n, perm);
std::vector<int> perm2 = heaps::unrank_o_n3(n, C);       // baseline version
```

The closed-form step table is available as `heaps::step_map(k)` (declared in
`heaps_unrank.hpp`); the baseline implementation and the test suite use it, the
main implementation does not need it.

## Resumable generation

Section 5 of the paper asks what the state of the iterative algorithm is, whether
`unrank` supplies it, and whether the machine can be started at an arbitrary
index without replaying the sequence from the identity. Two routes are provided;
their output sequences coincide.

The counters of the classical loop are a *state record* rather than the factorial
code: they are advanced after an emission (so at emission time they describe the
previous index), and the carry position is found by rescanning from `i = 0`.

```cpp
// Route A -- keep the classical loop, restore its counter array.
//   A = heaps::unrank(n, heaps::code_from_index(n, a));   // start array
//   c = heaps::restore_context(n, a);                     // its counters
//   i = 0;  // then run the classical loop: the first emission is index a + 1
std::vector<int> c = heaps::restore_context(n, a);

// Route B -- drive the loop by the factorial counter.
heaps::FactorialHeapGenerator gen(n, a);   // state = (array, code of a)
while (gen.step()) { /* gen.current(), gen.index(), gen.code() */ }
```

Route A does not touch the generator, at the price of keeping the state and the
index in step by conversion. Route B keeps them identical: `gen.code()` is the
code of the permutation currently held, and a restart is a copy. Progress per
permutation is the same in both routes (amortised O(1)); route B is not faster,
it is self-describing. Route B is the classical mixed-radix scheme known from
Arndt's `perm_heap2` (see References); the paper states that it is not new, and
says why the two routes are not interchangeable.

`test_generators` compares both routes with `heaps::unrank()` for every index at
n <= 8 -- 595,011 checks, ~0.3 s.

## What the verification suite checks

Everything is checked for **both** implementations, and the two are compared
against each other as well as against independent ground truth.

| # | Check | Scale |
|---|---|---|
| 1 | `step_map(k)` equals the closed form `F_k` of the paper | `k = 1..24` |
| 2 | the step table has the structure the main implementation relies on: transposition `(0,k-1)` for odd `k`, one `k`-cycle along the cycle order for even `k` | `k = 3..24` |
| 3 | **Exhaustive**: every code is compared with the iterative Heap algorithm (ground truth) *and* the two implementations must agree | `n = 1..9` -- 362,880 codes at `n = 9`; `--full` adds `n = 10` -- 3,628,800 codes |
| 4 | **Exhaustive** round trip `rank(unrank(C)) == C` for both implementations, plus `unrank(C) == unrank_o_n3(C)` | `n = 1..8` (all `n!` codes) |
| 5 | Randomised round trips for both implementations, plus a validity check that the output is a permutation of `0..n-1` | `n in {12,16,20,32,48,64,100}`, 20,000 trials each |

The generator suite is a separate program (`test_generators`, 595,011 checks);
it is described under *Resumable generation*.

The program exits with status 0 only when every check passes, and prints a
`PASS`/`FAIL` line per check. Check 3 is the strongest evidence: all `n!` codes
of order `n` are matched one by one against the permutation produced by Heap's
own iterative algorithm, for both implementations.

Measured on a development machine (clang 22, `-O2`): **26 checks, 0 failures in
~19 s** for the default run, **27 checks in ~34 s** with `--full`, and ~3 s with
`--quick`. Compilation is warning-free under `-Wall -Wextra -Wpedantic -Werror`,
and both suites run clean under AddressSanitizer and UndefinedBehaviorSanitizer.

## Complexity

| Direction | Baseline (dense, `_o_n3`) | Main | Space | Precomputation |
|---|---|---|---|---|
| `unrank` | Theta(n^3) | **Theta(n^2)** | `O(n)` | none |
| `rank` | Theta(n^3) | **Theta(n^2)** | `O(n)` | none |

The baseline applies the whole `step_map` table at every step, which costs
`Theta(k)` per step and `Theta(k^2)` per level. The main implementation applies
the same step in O(1) -- two swaps for odd `k-1`, one offset increment for even
`k-1` -- so a level costs `Theta(k)` and the total is `Theta(n^2)`. Both keep
only linear working arrays; no factorial-size lookup table is used. The paper
proves this in Section 8 (time) and Section 8.2 (space), and verifies the
resulting sequences in Section 7.

Sample wall-clock timings from `make bench` (single unrank call, best of five,
one random code per size, clang 22 `-O2`; machine-specific, for orientation only):

| `n` | baseline (ns) | main (ns) | ratio |
|---:|---:|---:|---:|
| 10 | 2 200 | 1 500 | 1.47x |
| 40 | 14 300 | 4 700 | 3.04x |
| 80 | 74 000 | 14 300 | 5.17x |
| 100 | 154 200 | 22 800 | 6.76x |
| 200 | 940 000 | 82 500 | 11.39x |

The ratio grows with `n`, in the direction implied by the saved `O(k)` factor --
at `n = 200` the quadratic version is already more than eleven times faster --
while at small `n` it is dominated by constants and cache effects. The paper
deliberately keeps its analysis asymptotic and does not report machine
timings.

## GitHub Actions

The pipeline in `.github/workflows/ci.yml` runs three jobs:

| Job | What it does |
|---|---|
| `build & test (os)` | configures and builds with CMake, then runs **both** CTest tests (the full verification suite and the generator suite) on `ubuntu-latest`, `macos-latest` and `windows-latest` |
| `strict warnings (cxx)` | compiles every source file, including the two test programs, with `-Wall -Wextra -Wpedantic -Werror` under `g++` and `clang++` |
| `ASan + UBSan` | rebuilds both suites with `-fsanitize=address,undefined` and runs them |

**How it runs**

* **Automatically** on every push to `main`/`master` and on every pull request.
* **Manually**: open the repository on GitHub, go to the **Actions** tab, select
  **CI** in the left-hand list, click **Run workflow**, choose a branch and
  confirm. (This is the `workflow_dispatch` trigger.)
* **Results**: the same Actions tab shows each run; a red cross marks a failing
  check, and the failing step's log is available under the job. The test programs
  print one `PASS`/`FAIL` line per check, so the log identifies the failure
  directly.

**Badge**: the badge at the top of this README points at
`.github/workflows/ci.yml`; it shows `passing` while the latest run of the
workflow on the default branch succeeded.

**Reproducing a CI job locally**

```bash
# job "build & test"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure

# job "strict warnings"
for f in heaps_unrank.cpp heaps_unrank_o_n3.cpp heaps_generator.cpp \
         bench_heaps.cpp test_heaps_unrank.cpp test_generators.cpp; do
  g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. -c "$f" -o "/tmp/${f%.cpp}.o"
done

# job "ASan + UBSan"
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_unrank.cpp heaps_unrank_o_n3.cpp test_heaps_unrank.cpp -o test_san && ./test_san
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_generator.cpp heaps_unrank.cpp heaps_unrank_o_n3.cpp test_generators.cpp \
    -o test_gen_san && ./test_gen_san
```

## License

Released under the MIT License; see [`LICENSE`](LICENSE).

## References

1. B. R. Heap. *Permutations by interchanges.* The Computer Journal 6(3):293-298, 1963.
2. B. R. Heap's method appears in R. Sedgewick, *Permutation generation methods*,
   Computing Surveys 9(2):137-164, 1977, and in D. E. Knuth, *TAOCP Vol. 4, Fascicle 2*, 2005.
3. J. Arndt. *Matters Computational: Ideas, Algorithms, Source Code.* Springer, Berlin,
   Heidelberg, 2011 -- `perm_heap2`, the mixed-radix generator of route B.
4. YUSHENG HU. *Ranking and Unranking Permutations in Heap's Order* (the companion paper).
