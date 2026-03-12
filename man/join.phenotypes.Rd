% Manual patch (Rscript/roxygen generation unavailable in this environment)
\name{join.phenotypes}
\alias{join.phenotypes}
\title{Attach phenotypes to a COTS spine table}
\usage{
join.phenotypes(
  spine_dt,
  pheno_dt,
  ph.ch = "ph.ch",
  ph.ma = "ph.ma",
  ph.pa = "ph.pa",
  proband = "w19_1011_lnr_k2_",
  extend = TRUE,
  pheno_ext = NULL,
  ph.ext = NULL,
  twins = NULL,
  NTR = FALSE,
  drop_na = NULL,
  missing_policy = c("available", "complete")
)
}
\arguments{
\item{spine_dt}{Spine table from \code{cots.spine.fun()}.}
\item{pheno_dt}{Phenotype source table.}
\item{ph.ch, ph.ma, ph.pa}{Column names for child/mother/father phenotypes.}
\item{proband}{Proband id column name.}
\item{extend}{Logical; add single-twin extensions if \code{TRUE}.}
\item{pheno_ext}{Phenotype table for twins without children.}
\item{ph.ext}{Phenotype column in \code{pheno_ext}.}
\item{twins}{Twin register table used when \code{extend = TRUE}.}
\item{NTR}{Logical; use NTR-specific twin columns in extension helper.}
\item{drop_na}{Deprecated compatibility argument.}
\item{missing_policy}{\code{"available"} (default) or \code{"complete"}.}
}
\value{
Data.table with columns \code{FAM1_F}, \code{FAM2_F}, \code{FAM1_M},
\code{FAM2_M}, \code{FAM1_CH1}, \code{FAM1_CH2}, \code{FAM2_CH1},
\code{FAM2_CH2}, \code{Gruppe}.
}
\description{
Joins child and parent phenotypes to a COTS spine, with sample-retaining
available-case behavior by default.
}
