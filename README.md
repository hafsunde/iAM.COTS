
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iAM.COTS

<!-- badges: start -->

<!-- badges: end -->

iAMCOTS implements an extended Children-of-Twins-and-Siblings (COTS) model with
indirect assortative mating (iAM). It includes a pair of twins/siblings, their spouses, and two offspring per nuclear family. The model is estimated in OpenMx and supports MZ, DZ, full-sibling, half-sibling, and unknown-zygosity family structures.

The model is described here:
>Sunde, H.F., Eilertsen, E.M. & Torvik, F.A. (2025) Understanding indirect assortative mating and its intergenerational consequences for educational attainment. *Nature Communications* https://doi.org/10.1038/s41467-025-60483-0

This repository was built after publication, mainly for personal use. See **Current status** below.

## Functions

The package currently exposes the following user-facing functions:

- `iam.cots.fun()` fits the main iAM-COTS OpenMx model from
  family-level input data.
- `get.varcomps()` extracts variance component estimates and
  confidence intervals from a fitted model.
- `get.POcor()` extracts the parent-offspring correlation
  decomposition from a fitted model.
- `iamcots.plot()` creates a summary plot of model results.
- `cots.spine.fun()` builds extended family units from cohort,
  ancestry, and twin-register data. The current selection procedure is
  deterministic, avoids re-using parent IDs across rows, and
  prioritizes family types in the order `MZ > DZ > UZ > HS > FS` (with
  `UZ` optional).
- `join.phenotypes()` attaches child and parent phenotypes to the
  spine created by `cots.spine.fun()`.
- `sex.invariant.fun()` remaps sex-specific parent columns into
  sex-invariant sibling/spouse columns for model input.
- `z.pheno()` applies pooled z-standardization across selected
  phenotype columns.
- `mx.WaldCI()` extracts Wald standard errors and confidence intervals
  from OpenMx objects.
- `combine_SE()` reshapes an estimate matrix and an SE matrix into long
  format.
- `combine_CI()` reshapes an estimate matrix and merges it with an
  OpenMx confidence interval table.

## Current status

The core modeling, extraction, plotting, and OpenMx utility functions
(`iam.cots.fun()`, `get.varcomps()`, `get.POcor()`, `iamcots.plot()`,
`mx.WaldCI()`, `combine_SE()`, and `combine_CI()`) have already been
used in the project in essentially their existing form, but the 
result-extraction helper functions may be needlessly hardcoded and brittle. 
These helper functions are mainly included because I have used them 
personally to help work with `iam.cots.fun()`.

By contrast, the data-preparation helper functions
(`cots.spine.fun()`, `join.phenotypes()`, `sex.invariant.fun()`, and
`z.pheno()`) have recently been
refurbished/refactored. They have unit tests and basic synthetic-data
checks, but they have **not yet been validated on real project data**. 
I originally used older tidyverse-based versions of these functions, 
which I have not uploaded here. 
Those helper functions should therefore still be treated as
provisional data preparation code until they have been cross-checked
against real data.
