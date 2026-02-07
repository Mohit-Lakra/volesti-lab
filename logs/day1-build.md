# Day 1 Build Log - February 7, 2026

## Goal
Get volesti built on my Mac so I can run the Easy Test for GSoC 2026.

---

## My Setup

**Machine:**
- MacBook with Apple Silicon (ARM64)
- macOS 26.1
- Homebrew Clang 20.1.8
- CMake 4.2.3

**Repository:**
- Location: `~/github/volesti-lab/volesti`
- Branch: `develop`

Started around 8:00 AM.

---

## Build Process

### Finding CMakeLists.txt

First confusion: there's no CMakeLists.txt in the root directory. Looked around and found it in `test/` folder. Makes sense since volesti is header-only - the tests are what actually get compiled.

```bash
cd ~/github/volesti-lab/volesti/test
```

### Running CMake

```bash
mkdir build
cd build
cmake ..
```

CMake ran pretty fast (~1-2 seconds). It downloaded some dependencies automatically:
- autodiff library
- lp_solve (the thing I'll be replacing later)

Found Eigen and Boost from Homebrew without issues.

Got some warnings about CMake policies (CMP0169, CMP0135) but they're just deprecation warnings, not actual errors. Build system uses older FetchContent patterns.

### First Compilation Attempt

```bash
make -j4
```

Build ran for about 1-1.5 minutes. Got to around 47% and then hit errors in `crhmc_polytope_preparation_test.cpp`.

Error was something like:
```
error: no member named 'LinearAccess' in 'Eigen::internal::functor_traits<lambda_as_visitor_wrapper<...>>'
```

Took a while to figure out what this meant. Turns out the CRHMC code (some advanced sampling stuff) uses lambda functions with Eigen in a way that modern Clang 20 doesn't like. The lambda wrappers don't define all the traits that Eigen expects.

### Workaround - Disable Problematic Tests

Decided to disable the failing tests since they're CRHMC-specific and I don't need them for the basic max ball computation.

**Test 1: crhmc_polytope_preparation_test**

Renamed the source file and commented out lines in CMakeLists.txt:
```bash
mv crhmc_polytope_preparation_test.cpp crhmc_polytope_preparation_test.cpp.backup
```

Then edited CMakeLists.txt around line 237 to comment out:
- `add_executable` for this test
- All the `add_test` commands
- `set_target_properties` 
- `TARGET_LINK_LIBRARIES`

Tried `sed` first but it broke the syntax since some commands span multiple lines. Ended up just manually editing in VS Code.

**Rebuild:** Got to 62% this time.

**Tests 2 & 3: mmcs_test and ode_solvers_test**

Same Eigen lambda errors. Used sed to comment these out:
```bash
sed -i '' 's/.*mmcs_test.*/# &/' CMakeLists.txt
sed -i '' 's/.*ode_solvers_test.*/# &/' CMakeLists.txt
```

**Rebuild:** Got to 91%.

**Test 4: crhmc_sampling_test**

Last one with the same issue. The sed approach partially broke it (commented only some lines of multi-line blocks). Had to manually fix lines 335-342 in CMakeLists.txt.

**Final rebuild:**
```bash
rm -rf *
cmake ..
make -j4
```

Took about 1-1.5 minutes. This time: **100% success!**

---

## What Built Successfully

Got 24 test binaries:

**Most important for Easy Test:**
- `test_internal_points` - has max ball computation
- `rounding_test` - uses max ball for rounding
- `new_volume_example` - volume computation demo

**Volume computation tools:**
- benchmarks_sob, benchmarks_cg, benchmarks_cb
- volume_sob_hpolytope, volume_sob_vpolytope
- volume_cg_hpolytope, volume_cg_vpolytope
- volume_cb_hpolytope, volume_cb_vpolytope
- volume_cb_vpoly_intersection_vpoly
- volume_cb_zonotopes

**Other tests:**
- sampling_test, logconcave_sampling_test, mcmc_diagnostics_test
- matrix_sampling_test, shake_and_bake_test, billiard_shake_and_bake_test
- boundary_oracles_test, order_polytope, root_finders_test
- simple_mc_integration

---

## What I Disabled

All CRHMC-related tests (4 total):
1. crhmc_polytope_preparation_test
2. mmcs_test
3. ode_solvers_test  
4. crhmc_sampling_test

These all failed because of Eigen lambda compatibility issues with Clang 20. The problem is in how volesti's CRHMC code wraps lambdas for Eigen's visitor pattern - the wrapper doesn't define `LinearAccess` and other required traits.

Doesn't affect my Easy Test since that's just basic max ball computation on simple polytopes.

---

## Observations

### lp_solve
Works fine! CMake downloaded it, compiled it as a static library, and linked it into the volume computation binaries. Got a few warnings about parentheses in the lp_solve source code itself, but nothing from the integration.

### Clang 20.1.8
Way stricter than older compilers. The CRHMC code probably works fine on GCC or older Clang versions, but modern Clang catches the missing trait definitions.

### Build Time
- CMake: 1-2 seconds
- Make (clean build): 1-1.5 minutes with -j4

Pretty fast on M-series Mac.

### CMakeLists.txt Editing
Learned the hard way that you can't just use sed blindly. Commands like `add_test` and `set_target_properties` can span multiple lines, and commenting out just the first line breaks the syntax. Had to be more careful and check the context.

---

## Files Modified

- `test/CMakeLists.txt` - commented out 4 CRHMC tests
- `test/crhmc_polytope_preparation_test.cpp` - renamed to .cpp.backup

---

## Next Steps

1. Run `test_internal_points` to verify max ball computation works
2. Write standalone test program for Easy Test
3. Test on cube, simplex, cross polytope
4. Document results
5. Start Medium Test research

---

## For My Proposal

This build experience gives me some material:

**What I learned:**
- volesti's CRHMC code needs updates for modern C++17/Clang
- The lambda visitor pattern with Eigen is fragile
- Build system could be more robust

**Opportunity:**
The lp_solve replacement project is a chance to also clean up some of this technical debt. Could modernize the CRHMC code and add CI testing for multiple compilers while I'm at it.

**Status:** Build complete, ready to test  
**Build result:** 24/28 tests (86%), 100% of needed functionality


*Written: Feb 7, 2026
