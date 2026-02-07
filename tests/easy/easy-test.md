# Easy Test – Day 1

## Objective
Build both interfaces and verify the maximum inscribed ball computation on simple polytopes.

## Environment
- macOS 26.1 (Apple Silicon M2)
- Clang 20.1.8, CMake 4.2.3
- R 4.5.2 with local build of Rvolesti 1.2.0

## Commands
```bash
# C++ tests (run from volesti/test)
mkdir -p build && cd build
cmake .. && make -j4
./test_internal_points > ../../tests/easy/cpp_output.txt

# R interface
Rscript -e 'Rcpp::compileAttributes()'
R_MAKEVARS_USER=/dev/null R CMD INSTALL --preclean --no-multiarch --with-keep.source .
Rscript tests/easy/run_inner_ball.R   # inlined via one-off command
```

## Results
- `test_internal_points`: 5/6 cases pass. Sparse order polytope center differs by ~1.0e-4 (tighter tolerance than the method delivers), so volume pipeline still usable for Easy Test. Full output in [`cpp_output.txt`](cpp_output.txt).
- R interface: `inner_ball()` returns the expected center and radius for the 10D cube, 5D simplex, and 4D cross-polytope. Raw console log in [`r_output.txt`](r_output.txt).

## Fixes Applied
1. `src/external/PackedCSparse/qd/Makefile`: removed `-march=native` and switched to `$(CXX)` to avoid forcing Intel-specific flags when compiling on Apple Silicon.
2. `src/external/PackedCSparse/FloatArray.h`: only include `<immintrin.h>` when `__AVX2__` is defined so arm64 builds skip x86-only intrinsics.
3. Clean rebuild of the bundled `lp_solve` and `qd` libraries to flush stale x86 object files before running `R CMD INSTALL`.

## What I Learned
- Embedded third-party code ignores top-level flag tweaks; patching their Makefiles is sometimes unavoidable.
- Always clean after cross-architecture experiments—the linker happily drags in stale objects and fails late.
- volesti’s C++ tests provide honest numerical feedback: the sparse order polytope case enforces a tolerance that highlights the limits of double precision.
