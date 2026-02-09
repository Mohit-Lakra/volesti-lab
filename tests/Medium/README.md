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

## C++ libraries

### 1) HiGHS (primary pick)

Repo: https://github.com/ERGO-Code/HiGHS

Why it looks like the best first integration target:

- License: **MIT**
- Active development (frequent commits / recent activity)
- Good performance reputation (especially on large sparse LPs)
- Supports presolve/scaling options and has a modern API

**Note:** I did read a few comments that “defaults” may not be ideal on numerically tricky instances. My takeaway is: enable presolve/scaling and be ready to tweak solver options (and add fallback behavior in volesti).

---

### 2) COIN-OR Clp (backup / reliability option)

Repo: https://github.com/coin-or/Clp

Why I kept this as a strong backup:

- License: **EPL** (permissive enough for most projects)
- Very mature / widely used in COIN-OR ecosystem
- Often described as “boring but dependable”
- Commonly used underneath other tools (CBC etc.)

Tradeoff: seems **slower than HiGHS** on a lot of modern benchmarks, but looks more “plug and play”.

---

