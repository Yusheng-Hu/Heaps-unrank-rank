# heaps_unrank -- Rank / Unrank for the Heap permutation generation order

[![CI](https://github.com/Yusheng-Hu/Heaps-unrank-rank/actions/workflows/ci.yml/badge.svg)](https://github.com/Yusheng-Hu/Heaps-unrank-rank/actions/workflows/ci.yml)
[![GitHub stars](https://img.shields.io/github/stars/Yusheng-Hu/Heaps-unrank-rank?style=social)](https://github.com/Yusheng-Hu/Heaps-unrank-rank)
![Field: Combinatorial Algorithms](https://img.shields.io/badge/Field-Combinatorial_Algorithms-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: C++](https://img.shields.io/badge/Language-C%2B%2B-red.svg)](#)
[![Time: Theta(n^2)](https://img.shields.io/badge/Time-Theta(n%5E2)-brightgreen.svg)](#)
[![Space: O(n)](https://img.shields.io/badge/Space-O(n)-green.svg)](#)

Rank and unrank for **Heap's permutation generation order** (B. R. Heap,
*Permutations by interchanges*, 1963). `unrank` maps an index in `0 .. n!-1` to
the permutation standing at that position of Heap's sequence; `rank` is its
inverse. Both take **Theta(n^2)** time and **O(n)** working space, and neither
uses a precomputed table of any kind.

A dense Theta(n^3) variant of the same construction is kept in
`heaps_unrank_o_n3.*` and is exercised by the test suite as an independent
cross-check.

The library also indexes a **second** recursion shape of the same family, in
which the swap comes after each recursive call (`heaps_unrank_variant.*`, see
*The swap-after-call order* below). It emits the same `n!` permutations in a
different order, also in **Theta(n^2)** time.

- **Repository:** <https://github.com/Yusheng-Hu/Heaps-unrank-rank>
- **Companion paper:** *Ranking and Unranking Permutations in Heap's Order* -- PDF
  in [`paper/heap_paper_en.pdf`](paper/heap_paper_en.pdf)
- **Author:** YUSHENG HU <dr.huyusheng@gmail.com>
- **Version:** 2.2
- **License:** MIT -- see [`LICENSE`](LICENSE)

---

## Using the library

```cpp
#include "heaps_unrank.hpp"
#include <cstdio>

int main() {
    const int n = 4;
    int C[n] = {0, 0, 0, 0};        // index: high-first Lehmer code, see below

    int perm[n];
    heaps::unrank(n, C, perm);      // index -> permutation

    int back[n];
    heaps::rank(n, perm, back);     // permutation -> index, back == C

    for (int i = 0; i < n; i++) std::printf("%d ", perm[i]);
    std::printf("\n");
}
```

`std::vector` overloads are provided as well:

```cpp
std::vector<int> perm = heaps::unrank(n, C);   // C is a std::vector<int>
std::vector<int> C2   = heaps::rank(n, perm);
```

The closed form of one in-block step is exposed as `heaps::step_map(k)`
(declared in `heaps_unrank.hpp`); it is the paper's `Step(k)` table.

## Index convention

The index is a **high-first Lehmer code** `C` of length `n`:

```
index = C[0]*(n-1)! + C[1]*(n-2)! + ... + C[n-2]*1! + C[n-1]*0!
0 <= C[i] <= n-i-1,        C[n-1] = 0
```

The paper writes the same index in low-order-first form `d`, related by
`d[i] = C[n-2-i]`; the digit read at level `k` is `C[n-k] = d[k-2]`.
The code space has `n!` elements, matching the `n!` permutations.

## How the algorithm works

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
  O(1) element moves instead of O(k). A level therefore costs `Theta(k)` and a
  whole call costs `Theta(n^2)`.

## Files

| File | Purpose |
|---|---|
| `heaps_unrank.hpp` | public interface: `unrank`, `rank`, `step_map`, vector overloads |
| `heaps_unrank.cpp` | implementation (per-step O(1) block transition, paper Section 4.5) |
| `heaps_unrank_variant.hpp` | public interface of the swap-after-call order: `base_map`, `unrank`, `rank` |
| `heaps_unrank_variant.cpp` | implementation of that order (closed-form `W(k, q)`, Theta(n^2)) |
| `heaps_generator.hpp` | resumable generation (paper, Section 5): `restore_context`, `FactorialHeapGenerator` |
| `heaps_generator.cpp` | the two routes that start the Heap sequence at an arbitrary index |
| `test_heaps_unrank.cpp` | verification suite |
| `test_variant.cpp` | verification suite of the swap-after-call order |
| `test_generators.cpp` | consistency test of both restart routes against `heaps::unrank()` (n <= 8) |
| `bench_heaps.cpp` | wall-clock benchmark of a single `unrank` call |
| `paper/heap_paper_en.pdf` | the companion paper (current version) |
| `Makefile` | simple build / test / bench targets |
| `CMakeLists.txt` | CMake build with three CTest tests |
| `.github/workflows/ci.yml` | GitHub Actions pipeline |

## Build and run

With `make`:

```bash
make            # build and run the standard suite + generator suite  (~20 s)
make full       # full suite: also exhaustively checks n = 10          (~34 s)
make quick      # quick smoke test                                     (~3 s)
make gen        # resumable-generation consistency test               (~0.3 s)
make variant    # swap-after-call order, exhaustive to n = 9            (~0.4 s)
make variant-full  # same, exhaustive to n = 11                        (~47 s)
make bench      # wall-clock benchmark
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

## What the verification suite checks

| # | Check | Scale |
|---|---|---|
| 1 | `step_map(k)` equals the closed form `F_k` of the paper | `k = 1..24` |
| 2 | the step table has the structure the implementation relies on: transposition `(0,k-1)` for odd `k`, one `k`-cycle along the cycle order for even `k` | `k = 3..24` |
| 3 | **Exhaustive**: every code is compared with the permutation produced by Heap's own iterative algorithm (independent ground truth) | `n = 1..9` -- 362,880 codes at `n = 9`; `--full` adds `n = 10` -- 3,628,800 codes |
| 4 | **Exhaustive** round trip `rank(unrank(C)) == C` for every code | `n = 1..8` (all `n!` codes) |
| 5 | Randomised round trips, plus a validity check that the output is a permutation of `0..n-1` | `n in {12,16,20,32,48,64,100}`, 20,000 trials each |

The generator suite is a separate program (`test_generators`) and the
swap-after-call suite another (`test_variant`); they are described under
*Resumable generation* and *The swap-after-call order*.

The program exits with status 0 only when every check passes, and prints a
`PASS`/`FAIL` line per check. Check 3 is the strongest evidence: all `n!` codes
of order `n` are matched one by one against the sequence that Heap's algorithm
itself emits.

Measured on a development machine (clang 22, `-O2`): **26 checks, 0 failures in
~19 s** for the default run, **27 checks in ~34 s** with `--full`, and ~3 s with
`--quick`. Compilation is warning-free under `-Wall -Wextra -Wpedantic -Werror`,
and both suites run clean under AddressSanitizer and UndefinedBehaviorSanitizer.

## Complexity

| Direction | Time | Space | Precomputation |
|---|---|---|---|
| `unrank` | **Theta(n^2)** | `O(n)` | none |
| `rank` | **Theta(n^2)** | `O(n)` | none |

The block transition is applied in O(1) rather than by walking the whole
`step_map` table, so a level costs `Theta(k)` and a whole call costs
`Theta(n^2)`. Only linear working arrays are kept; no factorial-size lookup
table is used. The paper proves the bound in Section 8 (time) and Section 8.2
(space), and verifies the resulting sequences in Section 7.

Sample wall-clock timings from `make bench` (one `unrank` call, best of five,
one random code per size, clang 22 `-O2`; machine-specific, for orientation
only):

| `n` | 10 | 20 | 40 | 60 | 80 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|
| `unrank`, us | 2.9 | 3.6 | 6.2 | 10.7 | 15.0 | 22.7 | 49.1 | 110.6 |

Doubling `n` from 100 to 200 multiplies the time by about five, which is the
quadratic growth the analysis predicts. The paper deliberately keeps its
analysis asymptotic and reports no machine timings.

## The swap-after-call order

Heap's recursion has a second common shape, with the swap placed after each
recursive call while the loop index runs over `0 .. k-1`:

```
gen(k):
    if k == 1: emit; return
    for i = 0 .. k-1:
        gen(k-1)
        if k is odd: swap(A[0], A[k-1]) else swap(A[i], A[k-1])
```

It emits the same `n!` permutations, but in a different order: the two
sequences agree for `n <= 3` and separate from `n = 4` on. `heaps_unrank_variant`
indexes that order with the same index convention as the main module:

```cpp
#include "heaps_unrank_variant.hpp"

std::vector<int> perm = heaps::variant::unrank(n, C);
std::vector<int> C2   = heaps::variant::rank(n, perm);   // C2 == C
```

The construction is unchanged -- one factoradic digit per level, `n` blocks of
`(n-1)!` at level `n` -- and the per-level map is even simpler. Block `q` of
level `k` starts with `W(k, q)` applied to the arrangement that entered the
level, and `W(k, q)` is closed form (`heaps::variant::base_map`); `rank` uses
its inverse, which is closed form too and is exposed as
`heaps::variant::inverse_map`. The two cases are:

| `k` | `W(k, q)` |
|---|---|
| odd | right rotation of the first `k` slots by `q` positions |
| even | the transpositions `(0,k-1)(1,k-1)...(q-1,k-1)`, in that order |

`W(k, q)` and its inverse are both closed form, so neither direction replays the
`q` steps of a block: a level costs `Theta(k)`, a call `Theta(n^2)`, and no
precomputation is used -- the same bounds as for Heap's order.

Two differences are worth knowing:

* the digit is read off the trailing element `t` of the current permutation as
  `b = (k-1-t) mod k` for odd `k`, and `b = 0` or `b = t+1` for even `k`;
* this order is **not minimal-change**: at some block boundaries consecutive
  outputs differ by more than one transposition (three places at `n = 4`,
  nineteen at `n = 5`). Heap's order keeps that property, which is why the
  swap-after-call shape is a separate order rather than the same algorithm.

`test_variant` checks the closed form of `W(k, q)` against the generator for
`k = 2..8` and that `W(k, q)` and `inverse_map(k, q)` are mutual inverses for
`k = 2..12` (649 combinations), then verifies both directions exhaustively in
one streaming pass: every emission is unranked from an independent mixed-radix
counter and ranked back, up to `n = 9` by default (409,110 permutations) and
`n = 11` with `--full` (43,954,710 permutations, ~47 s), plus round trips at
`n = 64, 128, 256`. The fast run reports **25 checks, 0 failures** in ~0.4 s and
the `--full` run **27 checks, 0 failures** in ~47 s. Two of the checks are
controls: the two orders must agree at `n <= 3` and must *disagree* from
`n = 4` on, so a test that accidentally compared one implementation with itself
would fail.

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
`n <= 8` -- 595,011 checks, ~0.3 s.

## GitHub Actions

The pipeline in `.github/workflows/ci.yml` runs three jobs:

| Job | What it does |
|---|---|
| `build & test (os)` | configures and builds with CMake, then runs all three CTest tests (verification, generator, swap-after-call) on `ubuntu-latest`, `macos-latest` and `windows-latest` |
| `strict warnings (cxx)` | compiles every source file, including the test programs, with `-Wall -Wextra -Wpedantic -Werror` under `g++` and `clang++` |
| `ASan + UBSan` | rebuilds all three suites with `-fsanitize=address,undefined` and runs them |

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
         bench_heaps.cpp test_heaps_unrank.cpp test_generators.cpp test_variant.cpp; do
  g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. -c "$f" -o "/tmp/${f%.cpp}.o"
done

# job "ASan + UBSan"
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_unrank.cpp heaps_unrank_o_n3.cpp test_heaps_unrank.cpp -o test_san && ./test_san
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_generator.cpp heaps_unrank.cpp heaps_unrank_o_n3.cpp test_generators.cpp \
    -o test_gen_san && ./test_gen_san
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -I. \
    heaps_unrank_variant.cpp heaps_unrank.cpp test_variant.cpp -o test_var_san && ./test_var_san
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
