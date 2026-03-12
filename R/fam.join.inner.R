#' Join Nuclear Families to a Single Relatedness Structure
#'
#' Internal helper used by `fam.join()`. The input `extended` must contain
#' exactly one `fam_link` level.
#'
#' @param extended Data.frame/data.table containing related pairs and a single
#'   `fam_link`.
#' @param nuclear Data.frame/data.table with nuclear family identifiers:
#'   `mor_lnr_k2_`, `far_lnr_k2_`, `CH1`, `CH2`.
#'
#' @return A data.table with columns `FAM1_*`, `FAM2_*`, and `Gruppe`.
#'
#' @keywords internal
#' @noRd
fam.join.inner <- function(extended, nuclear) {
  ext <- data.table::as.data.table(data.table::copy(extended))
  nuc <- data.table::as.data.table(data.table::copy(nuclear))

  if (!"fam_link" %in% names(ext)) {
    warning("fam_link column does not exist in `extended`.")
    return(NULL)
  }

  fam_link <- unique(ext[["fam_link"]])
  if (length(fam_link) != 1L) {
    warning("Input not separated by fam_link, use split(x, f = 'fam_link').")
    return(NULL)
  }
  fam_link <- fam_link[[1L]]

  drop_cols <- intersect(c("nuc_fam_id", "sib_size"), names(nuc))
  if (length(drop_cols)) {
    nuc[, (drop_cols) := NULL]
  }

  nuc_fam1 <- data.table::copy(nuc)
  data.table::setnames(
    nuc_fam1,
    old = c("mor_lnr_k2_", "far_lnr_k2_", "CH1", "CH2"),
    new = c("FAM1_M", "FAM1_F", "FAM1_CH1", "FAM1_CH2"),
    skip_absent = TRUE
  )

  nuc_fam2 <- data.table::copy(nuc)
  data.table::setnames(
    nuc_fam2,
    old = c("mor_lnr_k2_", "far_lnr_k2_", "CH1", "CH2"),
    new = c("FAM2_M", "FAM2_F", "FAM2_CH1", "FAM2_CH2"),
    skip_absent = TRUE
  )

  if (identical(fam_link, "mm")) {
    out <- ext[, .(Gruppe, FAM1_M = proband_FAM1, FAM2_M = proband_FAM2)]
    out <- merge(out, nuc_fam1, by = "FAM1_M", all.x = TRUE, sort = FALSE)
    out <- merge(out, nuc_fam2, by = "FAM2_M", all.x = TRUE, sort = FALSE)
  } else if (identical(fam_link, "ff")) {
    out <- ext[, .(Gruppe, FAM1_F = proband_FAM1, FAM2_F = proband_FAM2)]
    out <- merge(out, nuc_fam1, by = "FAM1_F", all.x = TRUE, sort = FALSE)
    out <- merge(out, nuc_fam2, by = "FAM2_F", all.x = TRUE, sort = FALSE)
  } else if (identical(fam_link, "fm")) {
    out <- ext[, .(Gruppe, FAM1_F = proband_FAM1, FAM2_M = proband_FAM2)]
    out <- merge(out, nuc_fam1, by = "FAM1_F", all.x = TRUE, sort = FALSE)
    out <- merge(out, nuc_fam2, by = "FAM2_M", all.x = TRUE, sort = FALSE)
  } else if (identical(fam_link, "mf")) {
    warning("Opposite-sex siblings are in wrong order. Is `mf` but should be `fm`.")
    return(NULL)
  } else if (identical(fam_link, "f_")) {
    out <- ext[, .(Gruppe, FAM1_F = proband_FAM1)]
    out <- merge(out, nuc_fam1, by = "FAM1_F", all.x = TRUE, sort = FALSE)
    out[, `:=`(FAM2_F = NA, FAM2_M = NA, FAM2_CH1 = NA, FAM2_CH2 = NA)]
  } else if (identical(fam_link, "m_")) {
    out <- ext[, .(Gruppe, FAM1_M = proband_FAM1)]
    out <- merge(out, nuc_fam1, by = "FAM1_M", all.x = TRUE, sort = FALSE)
    out[, `:=`(FAM2_F = NA, FAM2_M = NA, FAM2_CH1 = NA, FAM2_CH2 = NA)]
  } else {
    warning("fam_link not recognized. Use mm/ff/fm or f_/m_.")
    return(NULL)
  }

  needed <- c(
    "FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F",
    "FAM1_CH1", "FAM1_CH2", "FAM2_CH1", "FAM2_CH2", "Gruppe"
  )
  for (nm in setdiff(needed, names(out))) {
    out[, (nm) := NA]
  }

  out[, ..needed]
}
