# GSoC 2026 Test Solutions - Exclude LpSolve Project

**Applicant:** Mohit Lakra  
**Organization:** GeomScale  
**Project:** Replace LpSolve in volesti  

## What's This?

This repo contains my solutions for the GSoC 2026 qualification tests. I'm applying to work on replacing the LpSolve library in volesti with something more reliable for high-dimensional polytopes.

## Current Status

- ✅ **Easy Test** - Done (Feb 9)
- ✅ **Medium Test** - Done (Feb 9)  
- 🔄 **Hard Test** - Starting Feb 10

## Quick Links

- [Easy Test Results](tests/easy/)
- [Medium Test Research](tests/medium/)
- [Hard Test (TBD)](tests/hard/)
- [My Timeline](timeline.md)

## Main Findings So Far

After researching alternative LP solvers, I think **HiGHS** is the way to go for volesti:
- It's MIT licensed (no GPL issues like GLPK)
- Way faster than the current LpSolve
- Already used in SciPy as their default
- Has both C++ and R interfaces which is perfect

**Clp** would be a good backup since it's super reliable even if a bit slower.

## Test Environment

- **OS:**  MacOS
- **R Version:** 4.5.2
- **volesti:** Built from source (develop branch)
- **Rvolesti:** Latest from CRAN

Last updated: Feb 9, 2026
