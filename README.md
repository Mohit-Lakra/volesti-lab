# volesti-lab

Prep repo for the GeomScale GSoC 2026 idea “Exclude lp_solve from volesti”. Everything here is work-in-progress notes, logs, and tests that show I can build the stack and understand where the pain points are.

## Current Focus
- ✅ Day 1 (Easy Test groundwork): built C++ + R interfaces, captured max-ball outputs.
- 🔄 Day 2+: start solver research for the Medium Test.

## Day 1 Snapshot
- **C++ side:** `cmake .. && make -j4` inside `volesti/test`; kept 24 targets by commenting out the four CRHMC executables that choke on Clang 20. `./test_internal_points` now runs end-to-end (5/6 doctests pass, sparse order polytope fails only because tolerance is 1e-6 vs the achieved 1.0e-4).
- **R side:** `Rscript -e 'Rcpp::compileAttributes()'` followed by `R_MAKEVARS_USER=/dev/null R CMD INSTALL --preclean --no-multiarch --with-keep.source .`. This needed two local patches so the bundled `qd` code stops injecting `-march=native`/`<immintrin.h>` when the compiler is targeting arm64.
- **Verification:** Stored the doctest log, an `inner_ball()` script, and console output for cube/simplex/cross polytopes under `tests/easy/` so mentors can replay the Easy Test quickly.

## Key Fixes
1. `src/external/PackedCSparse/qd/Makefile` – drop `-march=native`, compile `.cc` sources with `$(CXX)` to inherit the right standard/arch from R.
2. `src/external/PackedCSparse/FloatArray.h` – guard `<immintrin.h>` behind `#if defined(__AVX2__)` so ARM builds no longer drag in x86-only headers.
3. Force-clean both embedded libraries plus `R CMD INSTALL --preclean` after switching architectures; otherwise stale x86 objects poison the final shared object.

## Logs & Proof
- `logs/day1-build.md` – raw build diary (C++ focus).
- `logs/day1-summary.md` – condensed checklist + commands.
- `tests/easy/easy-test.md` – Easy Test recipe with links to [`cpp_output.txt`](tests/easy/cpp_output.txt) and [`r_output.txt`](tests/easy/r_output.txt).

## Next Steps
- Medium Test: catalog LP/ball solvers (CRAN + modern C++), sketch replacement plan.
- Hard Test: prototype a replacement (likely nloptr + Eigen glue) once the research doc is ready.

*Last updated: 7 Feb 2026*
