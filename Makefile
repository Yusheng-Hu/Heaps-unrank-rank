# Makefile for heaps_unrank (rank / unrank for the Heap permutation order)
#
# Two implementations are built side by side:
#   heaps_unrank.cpp       main                   Theta(n^2)
#   heaps_unrank_o_n3.cpp  baseline (dense)       Theta(n^3)
#
# Targets:
#   make          build and run the standard verification suite
#   make test     same as above
#   make full     run the full suite (adds the exhaustive n = 10 check, ~40 s)
#   make quick    run a quick smoke test
#   make bench    build and run the wall-clock benchmark
#   make clean    remove the build directory

CXX      ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
BUILD    := build
ifeq ($(OS),Windows_NT)
  EXE    := .exe
endif
BIN      := $(BUILD)/test_heaps_unrank$(EXE)
BENCH    := $(BUILD)/bench_heaps$(EXE)
SOURCES  := heaps_unrank.cpp heaps_unrank_o_n3.cpp test_heaps_unrank.cpp
HEADERS  := heaps_unrank.hpp heaps_unrank_o_n3.hpp

.PHONY: all test full quick bench clean

all: test

$(BIN): $(SOURCES) $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -I. -o $@ $(SOURCES)

$(BENCH): bench_heaps.cpp heaps_unrank.cpp heaps_unrank_o_n3.cpp $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -I. -o $@ bench_heaps.cpp heaps_unrank.cpp heaps_unrank_o_n3.cpp

test: $(BIN)
	./$(BIN)

full: $(BIN)
	./$(BIN) --full

quick: $(BIN)
	./$(BIN) --quick

bench: $(BENCH)
	./$(BENCH)

clean:
	rm -rf $(BUILD)
