% Manual patch (Rscript/roxygen generation unavailable in this environment)
\name{z.pheno}
\alias{z.pheno}
\title{Pooled z-standardization across selected phenotype columns}
\usage{
z.pheno(dt, vars)
}
\arguments{
\item{dt}{Data.frame/data.table.}
\item{vars}{Character vector of columns to standardize.}
}
\value{
Data.table with standardized columns.
}
\description{
Computes one pooled mean/SD across selected variables, then standardizes each
selected column using the pooled values.
}
