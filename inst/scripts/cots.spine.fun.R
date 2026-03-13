#' Build COTS family spine from cohort, ancestry, and twin registry data
#'
#' Constructs extended-family units for COTS analyses, prioritizing twins
#' (`MZ -> DZ -> UZ`), then half siblings (HS), then full siblings (FS), while
#' enforcing one-use-per-parent/spouse across selected families.
#'
#' @param cohort Data.frame/data.table containing child and nuclear family IDs.
#'   Required columns: `w19_1011_lnr_k2_`, `mor_lnr_k2_`, `far_lnr_k2_`,
#'   `nuc_fam_id`.
#' @param ancestry Data.frame/data.table used for sibling pair identification.
#'   Required by `find.sibs()`.
#' @param twins Data.frame/data.table twin registry with at least:
#'   `proband`, `pair_id`, `sex`, `tw_size`, `SameSex`, `Zygo`.
#' @param opposite.sex Logical. If `FALSE`, only same-sex sibling pairs are used.
#' @param include.halfsibs Logical. If `TRUE`, include HS families and prioritize
#'   HS before FS in allocation.
#' @param include.UZ Logical. Include unknown-zygosity twins (`Zygo = NA`).
#'   Defaults to `FALSE`.
#' @param deterministic Logical. If `TRUE` (default), deterministic ordering is
#'   used throughout. If `FALSE`, tie-breaking can vary.
#' @param seed Optional integer seed used when `deterministic = FALSE`.
#'
#' @return Data.table with columns:
#'   `FAM1_M`, `FAM1_F`, `FAM2_M`, `FAM2_F`,
#'   `FAM1_CH1`, `FAM1_CH2`, `FAM2_CH1`, `FAM2_CH2`, `Gruppe`.
#'
#' @export
cots.spine.fun <- function(cohort,
                           ancestry,
                           twins,
                           opposite.sex = FALSE,
                           include.halfsibs = TRUE,
                           include.UZ = FALSE,
                           deterministic = TRUE,
                           seed = NULL) {
  if (!deterministic && !is.null(seed)) {
    set.seed(seed)
  }

  cohort_dt <- data.table::as.data.table(data.table::copy(cohort))
  ancestry_dt <- data.table::as.data.table(data.table::copy(ancestry))
  twins_dt <- data.table::as.data.table(data.table::copy(twins))

  req_cohort <- c("w19_1011_lnr_k2_", "mor_lnr_k2_", "far_lnr_k2_", "nuc_fam_id")
  req_ancestry <- c("proband", "mother", "father", "sex")
  req_twins <- c("proband", "pair_id", "sex", "tw_size", "SameSex", "Zygo")
  miss_cohort <- setdiff(req_cohort, names(cohort_dt))
  miss_ancestry <- setdiff(req_ancestry, names(ancestry_dt))
  miss_twins <- setdiff(req_twins, names(twins_dt))
  if (length(miss_cohort)) {
    stop("Missing columns in `cohort`: ", paste(miss_cohort, collapse = ", "))
  }
  if (length(miss_ancestry)) {
    stop("Missing columns in `ancestry`: ", paste(miss_ancestry, collapse = ", "))
  }
  if (length(miss_twins)) {
    stop("Missing columns in `twins`: ", paste(miss_twins, collapse = ", "))
  }

  message("Creating nuclear family units")
  cohort_subset <- cohort_dt[, ..req_cohort]
  if (!deterministic) {
    cohort_subset <- cohort_subset[sample.int(nrow(cohort_subset))]
  }
  data.table::setorderv(cohort_subset, c("nuc_fam_id", "w19_1011_lnr_k2_"))
  cohort_subset <- cohort_subset[, head(.SD, 2L), by = nuc_fam_id]
  cohort_subset[, sib_size := .N, by = nuc_fam_id]
  cohort_subset[, sib_n := seq_len(.N), by = nuc_fam_id]

  cohort_subset <- data.table::dcast(
    cohort_subset,
    mor_lnr_k2_ + far_lnr_k2_ + nuc_fam_id + sib_size ~ sib_n,
    value.var = "w19_1011_lnr_k2_"
  )
  if ("1" %in% names(cohort_subset)) data.table::setnames(cohort_subset, "1", "CH1")
  if ("2" %in% names(cohort_subset)) data.table::setnames(cohort_subset, "2", "CH2")
  if (!"CH1" %in% names(cohort_subset)) cohort_subset[, CH1 := NA]
  if (!"CH2" %in% names(cohort_subset)) cohort_subset[, CH2 := NA]

  data.table::setorderv(cohort_subset, c("sib_size", "nuc_fam_id"), order = c(-1L, 1L))
  cohort_subset <- cohort_subset[!duplicated(mor_lnr_k2_)]
  cohort_subset <- cohort_subset[!duplicated(far_lnr_k2_)]

  children_count <- data.table::rbindlist(
    list(
      cohort_subset[, .(person = mor_lnr_k2_, n_children = sib_size, implied_sex = 2L)],
      cohort_subset[, .(person = far_lnr_k2_, n_children = sib_size, implied_sex = 1L)]
    ),
    use.names = TRUE
  )

  parent_cols <- c("FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F")

  select_family_units <- function(candidates, nuclear, block_probands = character()) {
    cand <- data.table::as.data.table(data.table::copy(candidates))
    if (!nrow(cand)) {
      return(data.table::as.data.table(data.frame()))
    }
    if (length(block_probands)) {
      cand <- cand[
        !(proband_FAM1 %in% block_probands) &
          (!(proband_FAM2 %in% block_probands) | is.na(proband_FAM2))
      ]
    }
    if (!nrow(cand)) {
      return(data.table::as.data.table(data.frame()))
    }

    joined <- fam.join(cand, nuclear)
    if (is.null(joined) || !nrow(joined)) {
      return(data.table::as.data.table(data.frame()))
    }

    joined <- joined[
      !(is.na(FAM1_M) & is.na(FAM1_F)) &
        (grepl("^Single", Gruppe) | !(is.na(FAM2_M) & is.na(FAM2_F)))
    ]
    if (!nrow(joined)) {
      return(joined)
    }

    cots_select_nonoverlap(
      joined,
      id_cols = parent_cols,
      order_cols = c("Gruppe", "FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F"),
      deterministic = deterministic
    )
  }

  make_twin_group <- function(dt, label) {
    x <- data.table::as.data.table(data.table::copy(dt))
    if (!nrow(x)) {
      return(data.table::data.table(
        pair_id = character(),
        proband_FAM1 = character(),
        proband_FAM2 = character(),
        Gruppe = character(),
        fam_link = character()
      ))
    }
    x[, Gruppe := ifelse(
      !is.na(proband_FAM2),
      paste0(label, "_", role_FAM1, role_FAM2, "_", n_children_FAM1, n_children_FAM2),
      paste0("Single_", label, "_", role_FAM1, n_children_FAM1)
    )]
    x[, fam_link := ifelse(!is.na(proband_FAM2), paste0(role_FAM1, role_FAM2), paste0(role_FAM1, "_"))]
    x[, .(pair_id, proband_FAM1, proband_FAM2, Gruppe, fam_link)]
  }

  message("Creating extended family units with twins")
  tw <- twins_dt[proband %in% children_count$person & tw_size == 2]
  if (!opposite.sex) {
    tw <- tw[SameSex == 1 & (Zygo != 3 | is.na(Zygo))]
  }
  tw <- merge(tw, children_count, by.x = "proband", by.y = "person", all.x = TRUE, sort = FALSE)
  tw <- tw[sex == implied_sex]
  tw[, implied_sex := NULL]
  if (!deterministic && nrow(tw)) {
    tw <- tw[sample.int(nrow(tw))]
  }
  data.table::setorderv(tw, c("pair_id", "sex", "n_children", "proband"), order = c(1L, 1L, -1L, 1L))
  tw[, sib_nr := seq_len(.N), by = pair_id]
  tw <- tw[sib_nr <= 2L]
  tw[, role := data.table::fifelse(sex == 1, "f", data.table::fifelse(sex == 2, "m", NA_character_))]

  tw_wide <- tw[, .(
    proband_FAM1 = proband[1L],
    proband_FAM2 = proband[2L],
    role_FAM1 = role[1L],
    role_FAM2 = role[2L],
    n_children_FAM1 = n_children[1L],
    n_children_FAM2 = n_children[2L]
  ), by = .(pair_id, Zygo)]

  ancestry_subset_MZ <- make_twin_group(tw_wide[Zygo == 1], "MZ")
  ancestry_subset_DZ <- make_twin_group(tw_wide[Zygo == 2], "DZ")
  ancestry_subset_UZ <- if (isTRUE(include.UZ)) make_twin_group(tw_wide[is.na(Zygo)], "UZ") else data.table::as.data.table(data.frame())

  MZtwins_final <- select_family_units(ancestry_subset_MZ, cohort_subset)
  used_twins <- cots_collect_ids(MZtwins_final, parent_cols)

  DZtwins_final <- select_family_units(ancestry_subset_DZ, cohort_subset, block_probands = used_twins)
  used_twins <- unique(c(used_twins, cots_collect_ids(DZtwins_final, parent_cols)))

  UZtwins_final <- if (isTRUE(include.UZ)) {
    select_family_units(ancestry_subset_UZ, cohort_subset, block_probands = used_twins)
  } else {
    data.table::as.data.table(data.frame())
  }

  twins_final <- data.table::rbindlist(list(MZtwins_final, DZtwins_final, UZtwins_final), use.names = TRUE, fill = TRUE)
  twins_final_sample <- cots_collect_ids(twins_final, parent_cols)

  cohort_remaining <- cohort_subset[
    !(mor_lnr_k2_ %in% twins_final_sample) &
      !(far_lnr_k2_ %in% twins_final_sample)
  ]

  build_sibling_candidates <- function(type, block_probands = character()) {
    anc <- ancestry_dt[proband %in% children_count$person]
    if (length(twins_final_sample)) {
      anc <- anc[!(proband %in% twins_final_sample)]
    }
    if (length(block_probands)) {
      anc <- anc[!(proband %in% block_probands)]
    }
    if (!nrow(anc)) {
      return(data.table::as.data.table(data.frame()))
    }

    sib_long <- find.sibs(
      pedigree = anc,
      type = type,
      opposite.sex = opposite.sex
    )
    if (!nrow(sib_long)) {
      return(data.table::as.data.table(data.frame()))
    }

    if (identical(type, "FS") && "pair_id" %in% names(twins_dt)) {
      sib_long <- sib_long[!(pair_id %in% twins_dt$pair_id)]
    }

    sib_long <- merge(sib_long, children_count, by.x = "proband", by.y = "person", all.x = TRUE, sort = FALSE)
    sib_long <- sib_long[sex == implied_sex]
    sib_long[, implied_sex := NULL]
    if (!nrow(sib_long)) {
      return(data.table::as.data.table(data.frame()))
    }

    if (!deterministic) {
      sib_long <- sib_long[sample.int(nrow(sib_long))]
    }
    data.table::setorderv(sib_long, c("pair_id", "sex", "n_children", "proband"), order = c(1L, 1L, -1L, 1L))
    sib_long[, sib_nr := seq_len(.N), by = pair_id]
    sib_long <- sib_long[sib_nr <= 2L]
    sib_long[, role := data.table::fifelse(sex == 1, "f", data.table::fifelse(sex == 2, "m", NA_character_))]

    sib_wide <- sib_long[, .(
      proband_FAM1 = proband[1L],
      proband_FAM2 = proband[2L],
      role_FAM1 = role[1L],
      role_FAM2 = role[2L],
      n_children_FAM1 = n_children[1L],
      n_children_FAM2 = n_children[2L]
    ), by = pair_id]
    sib_wide <- sib_wide[!is.na(proband_FAM1) & !is.na(proband_FAM2)]
    if (!nrow(sib_wide)) {
      return(data.table::as.data.table(data.frame()))
    }

    prefix <- if (identical(type, "HS")) "HS" else "FS"
    sib_wide[, `:=`(
      Gruppe = paste0(prefix, "_", role_FAM1, role_FAM2, "_", n_children_FAM1, n_children_FAM2),
      fam_link = paste0(role_FAM1, role_FAM2)
    )]
    sib_wide[, .(pair_id, proband_FAM1, proband_FAM2, Gruppe, fam_link)]
  }

  message("Creating extended family units with other siblings")
  halfsibs_final <- data.table::as.data.table(data.frame())
  halfsibs_sample <- character()
  cohort_remaining2 <- cohort_remaining

  if (include.halfsibs) {
    hs_candidates <- build_sibling_candidates("HS")
    halfsibs_final <- select_family_units(hs_candidates, cohort_remaining)
    halfsibs_sample <- cots_collect_ids(halfsibs_final, parent_cols)
    cohort_remaining2 <- cohort_remaining[
      !(mor_lnr_k2_ %in% halfsibs_sample) &
        !(far_lnr_k2_ %in% halfsibs_sample)
    ]
  }

  fs_candidates <- build_sibling_candidates("FS", block_probands = halfsibs_sample)
  fullsibs_final <- select_family_units(fs_candidates, cohort_remaining2)

  message("Finalizing")
  if (include.halfsibs) {
    cots_spine <- data.table::rbindlist(list(twins_final, halfsibs_final, fullsibs_final), use.names = TRUE, fill = TRUE)
  } else {
    cots_spine <- data.table::rbindlist(list(twins_final, fullsibs_final), use.names = TRUE, fill = TRUE)
  }

  needed <- c(
    "FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F",
    "FAM1_CH1", "FAM1_CH2", "FAM2_CH1", "FAM2_CH2", "Gruppe"
  )
  for (nm in setdiff(needed, names(cots_spine))) {
    cots_spine[, (nm) := NA]
  }
  cots_spine <- cots_spine[, ..needed]
  cots_spine
}
