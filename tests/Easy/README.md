# Easy Test — Building and Running volesti (C++ + R)

**Status:** Completed  
**Date:** Feb 9, 2026  
**Author:** Mohit Lakra  

---

## Goal

Build both the **C++** and **R** versions of **volesti**, then sanity-check the **inner ball (Chebyshev ball)** computation on a few simple polytopes. The intention here is not to stress-test the solver yet, but to confirm the pipeline works end-to-end on small, well-understood shapes.

---

## Setup

### 1) C++ Build

Build steps:

```bash
cd volesti
mkdir build && cd build
cmake ..
make
```

After building, I ran the volume example and got reasonable output, so the C++ side looks fine.

---

### 2) R Install

In R:

```r
install.packages("volesti")
library(volesti)
```

No issues during install or load.

---

## Tests I Ran (Inner Ball)

### Test 1 — 3D Unit Cube `[0,1]^3` (H-representation)

Defined the cube using 6 inequalities (H-representation), then called `inner_ball()`.

**Output:**
- **Center:** (0.5, 0.5, 0.5)  
- **Radius:** 0.5  

This matches the expected geometry: the largest inscribed sphere in a unit cube is centered at the cube center with radius 0.5. ✅

---

### Test 2 — 2D Square `[0,1]^2` (H-representation)

Same idea in 2D.

**Output:**
- **Center:** (0.5, 0.5)  
- **Radius:** 0.5  

Also exactly what you’d expect. ✅

---

### Test 3 — 5D Hypercube `[-1,1]^5` (generated)

Used `gen_cube(5, 'H')`.

**Output:**
- **Center:** (0, 0, 0, 0, 0)  
- **Radius:** 1  

Since the hypercube is centered at the origin and extends to ±1 in each coordinate direction, radius 1 is correct. ✅

---

### Test 4 — 3D Simplex (generated)

Used `gen_simplex(3, 'H')`.

**Output:**
- **Center:** (0.2113249, 0.2113249, 0.2113249)  
- **Radius:** 0.2113249  

I didn’t verify this analytically, but the result looks plausible for a simplex. ✅

---

## Results Summary

All test cases ran successfully and returned sensible values.

For these small examples, the existing LP-based implementation (via **LpSolve**) behaves correctly. The motivation of the GSoC project is to improve robustness for **high-dimensional** / **ill-conditioned** inputs where LpSolve can fail or return unstable results.

---

## Code (R Script)

> File: `tests/easy/easy_test.R`

```r
# Easy Test - Inner Ball Computation
# Mohit Lakra
# Feb 9, 2026

library(volesti)

cat("=== Testing volesti Inner Ball ===\n\n")

# Test 1: 3D cube [0,1]^3
cat("Test 1: 3D Unit Cube\n")
A1 <- matrix(c(-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1),
             ncol=3, byrow=TRUE)
b1 <- c(0, 1, 0, 1, 0, 1)
P1 <- Hpolytope(A = A1, b = b1)
result1 <- inner_ball(P1)
cat("Result:", result1, "\n")
cat("Center: (", result1[1], result1[2], result1[3], "), Radius:", result1[4], "\n\n")

# Test 2: 2D square
cat("Test 2: 2D Square\n")
A2 <- matrix(c(-1,0, 1,0, 0,-1, 0,1), ncol=2, byrow=TRUE)
b2 <- c(0, 1, 0, 1)
P2 <- Hpolytope(A = A2, b = b2)
result2 <- inner_ball(P2)
cat("Result:", result2, "\n\n")

# Test 3: 5D hypercube
cat("Test 3: 5D Cube\n")
P3 <- gen_cube(5, 'H')
result3 <- inner_ball(P3)
cat("Result:", result3, "\n\n")

# Test 4: 3D simplex
cat("Test 4: 3D Simplex\n")
P4 <- gen_simplex(3, 'H')
result4 <- inner_ball(P4)
cat("Result:", result4, "\n\n")

cat("All tests completed successfully!\n")
```

---

## Output

> File: `tests/easy/results.txt`

```txt
=== Testing volesti Inner Ball ===

Test 1: 3D Unit Cube
Result: 0.5 0.5 0.5 0.5
Center: ( 0.5 0.5 0.5 ), Radius: 0.5

Test 2: 2D Square
Result: 0.5 0.5 0.5

Test 3: 5D Cube
Result: 0 0 0 0 0 1

Test 4: 3D Simplex
Result: 0.2113249 0.2113249 0.2113249 0.2113249

All tests completed successfully!
```

---

## Next Step

Move to the **Medium test**:
- research alternative LP solvers that can replace/exclude LpSolve,
- evaluate which options are most realistic for integration into **volesti / Rvolesti** (C++ core + R interface),
- and identify what changes are needed in the build + API layer.

---
