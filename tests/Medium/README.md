# Medium Test — LP Solver Research

**Status:** Completed  
**Date:** Feb 9, 2026  

---

## Task

Find at least **1 R package** and **2 C++ libraries** for solving **linear programs (LPs)**.

**Exclusions:** nloptr, quadprog, solve.QP, ipop.

---

## What I Found

### R package: `highs`

The CRAN package `highs` wraps the **HiGHS** solver. From what I checked, this is the cleanest option to start with because:

- **MIT license** (this matters for integration)
- Generally reported as **very fast** among open-source LP solvers
- **SciPy uses HiGHS** as its default LP backend now
- Same solver exists as a **native C++ library** too (so it fits volesti + Rvolesti nicely)

Install:

```r
install.packages("highs")
```

---
