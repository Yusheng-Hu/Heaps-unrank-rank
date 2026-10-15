# Makefile for heaps_unrank (rank / unrank for the Heap permutation order)
#
# Targets:
#   make          build and run the standard verification suite
#   make test     same as above
#   make full     run the full suite (adds the exhaustive n = 10 check, ~30 s)
#   make quick    run a fast smoke test
#   make clean    remove the build directory

CXX      ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
BUILD    := build
ifeq ($(OS),Windows_NT)
  EXE    := .exe
endif
BIN      := $(BUILD)/test_heaps_unrank$(EXE)
SOURCES  := heaps_unrank.cpp test_heaps_unrank.cpp
HEADERS  := heaps_unrank.hpp

.PHONY: all test full quick clean

all: test

$(BIN): $(SOURCES) $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -I. -o $@ $(SOURCES)

test: $(BIN)
	./$(BIN)

full: $(BIN)
	./$(BIN) --full

quick: $(BIN)
	./$(BIN) --quick

clean:
	rm -rf $(BUILD)
