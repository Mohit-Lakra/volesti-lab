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

## Why not GLPK?

GLPK kept appearing in searches, but the **GPL** licensing is the main blocker for this project.

volesti is LGPL-3.0, and pulling in a GPL solver can create license headaches (especially depending on how it’s linked/distributed). So even if GLPK works technically, it’s not the safest choice for keeping the licensing clean.

*(If mentors want a precise licensing interpretation, I’ll treat this as “needs confirmation” and won’t assume anything silently.)*

---


## Quick comparison (my current view)

| Solver | Speed | Reliability | License | Verdict |
|--------|-------|-------------|---------|---------|
| HiGHS | Very fast | Good (may need tuning) | MIT | Best default pick |
| Clp | Fast | Very good | EPL | Solid fallback |
| GLPK | Medium | OK | GPL | License risk |
| LpSolve | — | Unreliable in high-dim cases | LGPL | What we’re replacing |

---

## Research process (what I actually did)

This wasn’t just “read one blog and decide” — I tried to cross-check most claims:

- Skimmed the official GitHub repos (activity, docs, issues)
- Checked SciPy docs / release notes around the switch to HiGHS
- Looked at a couple solver benchmark writeups (not perfect, but enough signal)
- Verified CRAN availability and basic install story for R

---

## Recommendation for the project

Start with **HiGHS** as the primary solver, and keep **Clp** as an optional fallback.

Why this combo works well for volesti/Rvolesti:

1. Performance: HiGHS is usually the fastest open-source option
2. Reliability: Clp gives a “known stable” backup path
3. Licensing: both are permissive enough to avoid GPL headaches
4. Ecosystem: workable story for both **C++ core** and **R interface**

---

## Next steps

Hard test is next: implement the inner ball LP formulation (current plan mentions `nloptr`, but since that’s excluded above, I’ll treat that as “hard test requirement” rather than a solver recommendation).

---

# Quick Solver Comparison

## Top picks

1. **HiGHS** — fast, MIT, active
2. **Clp** — stable, mature, good fallback

## Feature checklist

| Feature | HiGHS | Clp | GLPK | LpSolve |
|--------|-------|-----|------|--------|
| Speed | ⚡⚡⚡ | ⚡⚡ | ⚡ | 💀 |
| Reliability | Good* | Excellent | OK | Bad |
| License | MIT ✓ | EPL ✓ | GPL ✗ | LGPL |
| High-dim support | ✓✓ | ✓✓ | ✓ | ✗ |
| Active dev | ✓✓✓ | ✓ | Slow | Dead |
| R support | `highs` | via ROI | yes | yes |

\* HiGHS usually benefits from presolve/scaling (especially on harder instances).

---

## Benchmarks (notes-to-self)

I saw at least one “robustness style” benchmark on a set of hard LPs where Clp solved more out-of-the-box, and HiGHS needed option tweaks to behave well.  
I’m not treating those numbers as absolute truth — they’re just a reminder that defaults matter.

---

## Installation notes

```r
# R packages
install.packages("highs")

# For Clp access in R (via ROI ecosystem)
install.packages("ROI")
# (and potentially ROI.plugin.clp depending on setup)
```
