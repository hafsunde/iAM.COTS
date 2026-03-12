% Manual patch (Rscript/roxygen generation unavailable in this environment)
\name{sex.invariant.fun}
\alias{sex.invariant.fun}
\title{Recode family columns to sex-invariant sibling/spouse slots}
\usage{
sex.invariant.fun(dta)
}
\arguments{
\item{dta}{Data.frame/data.table with \code{Gruppe} and family phenotype columns.}
}
\value{
Data.table with sex-invariant sibling/spouse slots.
}
\description{
Maps sex-specific family columns to invariant \code{FAM*_Sib} and
\code{FAM*_Spouse} columns.
}
