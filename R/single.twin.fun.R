#' Attach phenotypes for single-twin family structures
#'
#' Internal helper for `join.phenotypes()`. It links childless co-twins and
#' enriches "Single_*" groups with available phenotypes.
#'
#' @param spine_dt Spine table containing family ids (with `_id` suffixes).
#' @param pheno Phenotype table with columns `proband`, `CH`, `Ma`, `Fa`.
#' @param pheno_ext Phenotype table for twins without children.
#' @param ph.ext Column name in `pheno_ext` containing the phenotype.
#' @param twins Twin register table.
#' @param ph.ch,ph.ma,ph.pa Unused legacy placeholders (kept for compatibility).
#' @param proband Proband id column in `pheno_ext`.
#' @param NTR Logical. If `TRUE`, uses NTR-style twin column names.
#'
#' @return Data.table with the standard eight family phenotype columns plus
#'   `Gruppe`.
#'
#' @keywords internal
#' @noRd
single.twin.fun <- function(spine_dt,
                            pheno,
                            pheno_ext,
                            ph.ext,
                            twins,
                            ph.ch = "ph.ch",
                            ph.ma = "ph.ma",
                            ph.pa = "ph.pa",
                            proband = "w19_1011_lnr_k2_",
                            NTR = FALSE) {
  sp <- data.table::as.data.table(data.table::copy(spine_dt))
  ph <- data.table::as.data.table(data.table::copy(pheno))
  tw <- data.table::as.data.table(data.table::copy(twins))
  ext <- data.table::as.data.table(data.table::copy(pheno_ext))
  ext <- ext[, .(proband = get(proband), phenotype = get(ph.ext))]

  mothers <- sp[grepl("Single", Gruppe) & grepl("m", Gruppe)]
  fathers <- sp[grepl("Single", Gruppe) & grepl("f", Gruppe)]

  parent_ids <- unique(c(mothers$FAM1_M_id, fathers$FAM1_F_id))
  parent_ids <- parent_ids[!is.na(parent_ids)]

  if (!length(parent_ids)) {
    return(data.table::data.table(
      FAM1_F = numeric(),
      FAM2_F = numeric(),
      FAM1_M = numeric(),
      FAM2_M = numeric(),
      FAM1_CH1 = numeric(),
      FAM1_CH2 = numeric(),
      FAM2_CH1 = numeric(),
      FAM2_CH2 = numeric(),
      Gruppe = character()
    ))
  }

  if (!NTR) {
    partial_pairs <- unique(
      tw[SameSex == 1 & proband %in% parent_ids, pair_id]
    )
    missing_twins <- tw[pair_id %in% partial_pairs]
    missing_twins[, no_kid := data.table::fifelse(proband %in% parent_ids, 0L, 1L)]
    data.table::setorderv(missing_twins, c("pair_id", "no_kid", "proband"))
    missing_twins <- missing_twins[, .(
      Twin_0 = proband[1L],
      Twin_1 = proband[2L]
    ), by = pair_id]
  } else {
    partial_pairs <- unique(
      tw[LNr_SSB_k2_ %in% parent_ids, ParNr]
    )
    missing_twins <- tw[ParNr %in% partial_pairs]
    missing_twins[, no_kid := data.table::fifelse(LNr_SSB_k2_ %in% parent_ids, 0L, 1L)]
    data.table::setorderv(missing_twins, c("ParNr", "no_kid", "LNr_SSB_k2_"))
    missing_twins <- missing_twins[, .(
      Twin_0 = LNr_SSB_k2_[1L],
      Twin_1 = LNr_SSB_k2_[2L]
    ), by = ParNr]
    data.table::setnames(missing_twins, "ParNr", "pair_id")
  }
  missing_twins[, pair_id := NULL]

  child1 <- data.table::copy(ph)
  data.table::setnames(
    child1,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM1_CH1_id", "FAM1_CH1", "FAM1_M1", "FAM1_F1")
  )
  child2 <- data.table::copy(ph)
  data.table::setnames(
    child2,
    c("proband", "CH", "Ma", "Fa"),
    c("FAM1_CH2_id", "FAM1_CH2", "FAM1_M2", "FAM1_F2")
  )

  build_single <- function(dta, fam1_id, fam2_id, role = c("m", "f")) {
    role <- match.arg(role)
    if (!nrow(dta)) {
      return(data.table::data.table(
        FAM1_F = numeric(),
        FAM2_F = numeric(),
        FAM1_M = numeric(),
        FAM2_M = numeric(),
        FAM1_CH1 = numeric(),
        FAM1_CH2 = numeric(),
        FAM2_CH1 = numeric(),
        FAM2_CH2 = numeric(),
        Gruppe = character()
      ))
    }

    x <- merge(dta, missing_twins, by.x = fam1_id, by.y = "Twin_0", all.x = TRUE, sort = FALSE)
    if (fam2_id %in% names(x)) {
      x[, (fam2_id) := NULL]
    }
    x[, (fam2_id) := Twin_1]
    x[, Twin_1 := NULL]

    x <- merge(x, child1, by = "FAM1_CH1_id", all.x = TRUE, sort = FALSE)
    x <- merge(x, child2, by = "FAM1_CH2_id", all.x = TRUE, sort = FALSE)

    x[, FAM1_M := rowMeans(cbind(FAM1_M1, FAM1_M2), na.rm = TRUE)]
    x[, FAM1_F := rowMeans(cbind(FAM1_F1, FAM1_F2), na.rm = TRUE)]
    x[is.nan(FAM1_M), FAM1_M := NA_real_]
    x[is.nan(FAM1_F), FAM1_F := NA_real_]

    x[, c("FAM1_M1", "FAM1_M2", "FAM1_F1", "FAM1_F2") := NULL]

    x <- merge(x, ext, by.x = fam2_id, by.y = "proband", all.x = TRUE, sort = FALSE)
    x <- x[!is.na(phenotype)]

    if (role == "m") {
      x[, `:=`(
        FAM2_M = phenotype,
        FAM2_F = NA_real_,
        FAM2_CH1 = NA_real_,
        FAM2_CH2 = NA_real_
      )]
      x[, Gruppe := paste0(substr(Gruppe, 8, 11), "m_", substr(Gruppe, 12, 12), "0")]
    } else {
      x[, `:=`(
        FAM2_F = phenotype,
        FAM2_M = NA_real_,
        FAM2_CH1 = NA_real_,
        FAM2_CH2 = NA_real_
      )]
      x[, Gruppe := paste0(substr(Gruppe, 8, 11), "f_", substr(Gruppe, 12, 12), "0")]
    }

    out <- x[, .(
      FAM1_F, FAM2_F, FAM1_M, FAM2_M,
      FAM1_CH1, FAM1_CH2, FAM2_CH1, FAM2_CH2, Gruppe
    )]

    n_flip <- floor(nrow(out) / 2)
    if (n_flip > 0) {
      idx <- seq_len(n_flip)
      tmp <- data.table::copy(out[idx])
      out[idx, `:=`(
        FAM1_F = tmp$FAM2_F,
        FAM2_F = tmp$FAM1_F,
        FAM1_M = tmp$FAM2_M,
        FAM2_M = tmp$FAM1_M,
        FAM1_CH1 = tmp$FAM2_CH1,
        FAM1_CH2 = tmp$FAM2_CH2,
        FAM2_CH1 = tmp$FAM1_CH1,
        FAM2_CH2 = tmp$FAM1_CH2
      )]
    }

    out
  }

  mothers_pheno <- build_single(mothers, fam1_id = "FAM1_M_id", fam2_id = "FAM2_M_id", role = "m")
  fathers_pheno <- build_single(fathers, fam1_id = "FAM1_F_id", fam2_id = "FAM2_F_id", role = "f")

  data.table::rbindlist(list(mothers_pheno, fathers_pheno), use.names = TRUE, fill = TRUE)
}
