
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

By contrast, recent data-preparation helper functions are still provisional and
have now been moved out of the package namespace to avoid masking similarly
named objects when loading `iAM.COTS`.

### New data-preparation script location

The following helpers are now available as scripts under `inst/scripts/` rather
than exported package functions:

- `cots.spine.fun()`
- `join.phenotypes()`
- `sex.invariant.fun()`
- `z.pheno()`
- internal support helpers used by those functions (for example
  `find.sibs()`, `fam.join()`, and `pair.reducer()`).

These scripts are included in the installed package for reference/reuse, while
the attached package API remains focused on the validated modeling and output
extraction functions.
See `NEWS.md` for a changelog entry listing these newly relocated helpers.
