#' Attach phenotypes to a COTS spine table
#'
#' Joins child and parent phenotypes to a spine produced by `cots.spine.fun()`.
#'
#' @param spine_dt Spine table with family IDs (`FAM1_*`, `FAM2_*`) and `Gruppe`.
#' @param pheno_dt Phenotype source table.
#' @param ph.ch,ph.ma,ph.pa Column names in `pheno_dt` for child, mother, and
#'   father phenotypes.
#' @param proband Column name in `pheno_dt` with person IDs.
#' @param extend Logical; if `TRUE`, attempts to add single-twin extensions
#'   using `pheno_ext`.
#' @param pheno_ext Table with phenotypes for twins without children.
#' @param ph.ext Column name in `pheno_ext` containing the phenotype.
#' @param twins Twin register required when `extend = TRUE`.
#' @param NTR Logical. If `TRUE`, uses NTR-specific twin columns in
#'   `single.twin.fun()`.
#' @param drop_na Deprecated compatibility argument. If provided, it overrides
#'   `missing_policy`.
#' @param missing_policy Missing-data policy for `pheno_dt`:
#'   `"available"` (default, keeps partial data) or `"complete"` (complete-case).
#'
#' @return Data.table with columns:
#'   `FAM1_F`, `FAM2_F`, `FAM1_M`, `FAM2_M`,
#'   `FAM1_CH1`, `FAM1_CH2`, `FAM2_CH1`, `FAM2_CH2`, `Gruppe`.
#'
#' @export
join.phenotypes <- function(spine_dt,
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
                            missing_policy = c("available", "complete")) {
  missing_policy <- match.arg(missing_policy)
  if (!is.null(drop_na)) {
    warning("`drop_na` is deprecated; use `missing_policy`.", call. = FALSE)
    missing_policy <- if (isTRUE(drop_na)) "complete" else "available"
  }

  sp <- data.table::as.data.table(data.table::copy(spine_dt))
  id_candidates <- names(sp)[seq_len(min(8L, ncol(sp)))]
  for (nm in id_candidates) {
    if (!grepl("_id$", nm)) {
      data.table::setnames(sp, nm, paste0(nm, "_id"))
    }
  }

  ph <- data.table::as.data.table(data.table::copy(pheno_dt))
  vars <- c(proband, ph.ch, ph.ma, ph.pa)
  missing_cols <- setdiff(vars, names(ph))
  if (length(missing_cols)) {
    stop("Missing phenotype columns: ", paste(missing_cols, collapse = ", "))
  }
  ph <- ph[, .(
    proband = get(proband),
    CH = get(ph.ch),
    Ma = get(ph.ma),
    Fa = get(ph.pa)
  )]

  if (identical(missing_policy, "complete")) {
    ph <- ph[stats::complete.cases(ph)]
  }

  if (extend) {
    if (is.null(pheno_ext)) {
      stop("Missing `pheno_ext` while `extend = TRUE`.")
    }
    if (is.null(ph.ext)) {
      stop("Missing `ph.ext` while `extend = TRUE`.")
    }
    if (is.null(twins)) {
      stop("Missing `twins` while `extend = TRUE`.")
    }
    single_phenotypes <- single.twin.fun(
      spine_dt = sp,
      pheno = ph,
      pheno_ext = pheno_ext,
      ph.ext = ph.ext,
      twins = twins,
      proband = proband,
      NTR = NTR
    )
  }

  dt <- sp[!grepl("Single", Gruppe)]

  ph_fam1_ch1 <- data.table::copy(ph)
  data.table::setnames(
    ph_fam1_ch1,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM1_CH1_id", "FAM1_CH1", "FAM1_M1", "FAM1_F1")
  )
  ph_fam1_ch2 <- data.table::copy(ph)
  data.table::setnames(
    ph_fam1_ch2,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM1_CH2_id", "FAM1_CH2", "FAM1_M2", "FAM1_F2")
  )
  ph_fam2_ch1 <- data.table::copy(ph)
  data.table::setnames(
    ph_fam2_ch1,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM2_CH1_id", "FAM2_CH1", "FAM2_M1", "FAM2_F1")
  )
  ph_fam2_ch2 <- data.table::copy(ph)
  data.table::setnames(
    ph_fam2_ch2,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM2_CH2_id", "FAM2_CH2", "FAM2_M2", "FAM2_F2")
  )

  dt <- merge(dt, ph_fam1_ch1, by = "FAM1_CH1_id", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, ph_fam1_ch2, by = "FAM1_CH2_id", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, ph_fam2_ch1, by = "FAM2_CH1_id", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, ph_fam2_ch2, by = "FAM2_CH2_id", all.x = TRUE, sort = FALSE)

  dt[, `:=`(
    FAM1_M = data.table::fifelse(!is.na(FAM1_M1), FAM1_M1, FAM1_M2),
    FAM1_F = data.table::fifelse(!is.na(FAM1_F1), FAM1_F1, FAM1_F2),
    FAM2_M = data.table::fifelse(!is.na(FAM2_M1), FAM2_M1, FAM2_M2),
    FAM2_F = data.table::fifelse(!is.na(FAM2_F1), FAM2_F1, FAM2_F2)
  )]

  drop_cols <- grep("_id$", names(dt), value = TRUE)
  aux_cols <- c("FAM1_M1", "FAM1_M2", "FAM1_F1", "FAM1_F2", "FAM2_M1", "FAM2_M2", "FAM2_F1", "FAM2_F2")
  dt[, c(intersect(drop_cols, names(dt)), intersect(aux_cols, names(dt))) := NULL]

  out <- dt[, .(
    FAM1_F, FAM2_F, FAM1_M, FAM2_M,
    FAM1_CH1, FAM1_CH2, FAM2_CH1, FAM2_CH2, Gruppe
  )]

  if (extend) {
    out <- data.table::rbindlist(list(out, single_phenotypes), use.names = TRUE, fill = TRUE)
  }

  out
}
