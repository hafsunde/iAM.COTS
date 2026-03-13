#' Recode family columns to sex-invariant sibling/spouse slots
#'
#' Maps sex-specific parent columns to invariant slots:
#' `FAM*_Sib` and `FAM*_Spouse`.
#'
#' @param dta Data.frame/data.table containing `Gruppe` and family phenotype
#'   columns.
#'
#' @return Data.table with recoded sibling/spouse columns.
#'
#' @export
sex.invariant.fun <- function(dta) {
  dt <- data.table::as.data.table(data.table::copy(dta))

  make_part <- function(subset_dt, ren_old, ren_new) {
    if (!nrow(subset_dt)) {
      return(subset_dt)
    }
    data.table::setnames(subset_dt, old = ren_old, new = ren_new, skip_absent = TRUE)
    subset_dt
  }

  dta_ff <- make_part(
    dt[grepl("ff", Gruppe)],
    ren_old = c("FAM1_F", "FAM2_F", "FAM1_M", "FAM2_M"),
    ren_new = c("FAM1_Sib", "FAM2_Sib", "FAM1_Spouse", "FAM2_Spouse")
  )
  dta_mm <- make_part(
    dt[grepl("mm", Gruppe)],
    ren_old = c("FAM1_M", "FAM2_M", "FAM1_F", "FAM2_F"),
    ren_new = c("FAM1_Sib", "FAM2_Sib", "FAM1_Spouse", "FAM2_Spouse")
  )
  dta_fm <- make_part(
    dt[grepl("fm", Gruppe)],
    ren_old = c("FAM1_F", "FAM2_M", "FAM1_M", "FAM2_F"),
    ren_new = c("FAM1_Sib", "FAM2_Sib", "FAM1_Spouse", "FAM2_Spouse")
  )
  dta_mf <- make_part(
    dt[grepl("mf", Gruppe)],
    ren_old = c("FAM1_M", "FAM2_F", "FAM1_F", "FAM2_M"),
    ren_new = c("FAM1_Sib", "FAM2_Sib", "FAM1_Spouse", "FAM2_Spouse")
  )

  parts <- list(dta_ff, dta_mm, dta_fm, dta_mf)
  parts <- parts[vapply(parts, nrow, integer(1L)) > 0L]
  if (!length(parts)) {
    return(data.table::as.data.table(data.frame()))
  }

  out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  ref_cols <- names(parts[[1L]])
  out[, ..ref_cols]
}
