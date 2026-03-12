% Manual patch (Rscript/roxygen generation unavailable in this environment)
\name{cots.spine.fun}
\alias{cots.spine.fun}
\title{Build COTS family spine from cohort, ancestry, and twin registry data}
\usage{
cots.spine.fun(
  cohort,
  ancestry,
  twins,
  opposite.sex = FALSE,
  include.halfsibs = TRUE,
  include.UZ = FALSE,
  deterministic = TRUE,
  seed = NULL
)
}
\arguments{
\item{cohort}{Data.frame/data.table with \code{w19_1011_lnr_k2_}, \code{mor_lnr_k2_}, \code{far_lnr_k2_}, \code{nuc_fam_id}.}
\item{ancestry}{Pedigree-like table used by \code{find.sibs()}.}
\item{twins}{Twin register table with \code{proband}, \code{pair_id}, \code{sex}, \code{tw_size}, \code{SameSex}, \code{Zygo}.}
\item{opposite.sex}{Logical; if \code{FALSE}, only same-sex sibling pairs are used.}
\item{include.halfsibs}{Logical; include HS families and prioritize HS before FS.}
\item{include.UZ}{Logical; include unknown-zygosity twins (\code{Zygo = NA}).}
\item{deterministic}{Logical; deterministic ordering/tie-breaks.}
\item{seed}{Optional seed used when \code{deterministic = FALSE}.}
}
\value{
Data.table with columns \code{FAM1_M}, \code{FAM1_F}, \code{FAM2_M},
\code{FAM2_F}, \code{FAM1_CH1}, \code{FAM1_CH2}, \code{FAM2_CH1},
\code{FAM2_CH2}, and \code{Gruppe}.
}
\description{
Constructs extended-family units for COTS analyses with priority
\code{MZ -> DZ -> UZ} (UZ optional), then HS before FS.
}
