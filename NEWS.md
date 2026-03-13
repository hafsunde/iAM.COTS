# iAM.COTS news

## iAM.COTS 0.0.0.9003

- Moved provisional data-preparation helpers out of `R/` and into
  `inst/scripts/` so they remain available in the package source/install but are
  not attached as package functions on load.
- Scripted helpers now located in `inst/scripts/` include:
  `cots.spine.fun()`, `join.phenotypes()`, `sex.invariant.fun()`, `z.pheno()`,
  and related support helpers (`find.sibs()`, `fam.join()`, `pair.reducer()`,
  etc.).
- Exported package API now focuses on model fitting, extraction, plotting, and
  OpenMx summary helpers.
