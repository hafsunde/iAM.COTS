% Manual patch (Rscript/roxygen generation unavailable in this environment)
\name{find.sibs}
\alias{find.sibs}
\title{Identify full- and half-sibling pairs in a pedigree}
\usage{
find.sibs(
  pedigree,
  type,
  proband = "proband",
  mother = "mother",
  father = "father",
  sex = "sex",
  exclusive = FALSE,
  opposite.sex = FALSE
)
}
\arguments{
\item{pedigree}{Data.frame/data.table with pedigree columns.}
\item{type}{\code{"FS"}, \code{"HS"}, or \code{"both"}.}
\item{proband}{Column name for proband id.}
\item{mother}{Column name for mother id.}
\item{father}{Column name for father id.}
\item{sex}{Column name for sex.}
\item{exclusive}{Logical; if \code{TRUE}, each proband appears at most once.}
\item{opposite.sex}{Logical; if \code{FALSE}, keep only same-sex pairs.}
}
\value{
Data.table with columns \code{pair_id}, \code{relation}, \code{proband}, \code{sex}.
}
\description{
Builds sibling pairs using full IDs and canonical pair keys (no substring truncation).
}
