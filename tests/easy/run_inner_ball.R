#!/usr/bin/env Rscript
library(volesti)
cat("cube10\n")
print(inner_ball(gen_cube(10, "H")))
cat("simplex5\n")
print(inner_ball(gen_simplex(5, "H")))
cat("cross4\n")
print(inner_ball(gen_cross(4, "H")))
