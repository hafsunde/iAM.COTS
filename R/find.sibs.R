#' Identify Full- and Half-Sibling Pairs in a Pedigree
#'
#' Builds sibling pairs from parent identifiers and returns one row per person
#' per pair (`pair_id`, `relation`, `proband`, `sex`).
#'
#' @param pedigree Data.frame/data.table containing pedigree information.
#' @param type Character scalar: `"FS"`, `"HS"`, or `"both"`.
#' @param proband Column name for proband id.
#' @param mother Column name for mother id.
#' @param father Column name for father id.
#' @param sex Column name for sex.
#' @param exclusive Logical. If `TRUE`, each proband can appear in at most one
#'   pair (deterministic greedy reduction).
#' @param opposite.sex Logical. If `FALSE`, only same-sex sibling pairs are
#'   returned.
#'
#' @return Data.table with columns `pair_id`, `relation`, `proband`, and `sex`.
#'
#' @details
#' Pair ids are based on full identifiers (canonicalized with `pmin/pmax` logic),
#' not truncated substrings.
#'
#' @export
find.sibs <- function(pedigree,
                      type,
                      proband = "proband",
                      mother = "mother",
                      father = "father",
                      sex = "sex",
                      exclusive = FALSE,
                      opposite.sex = FALSE) {
  if (missing(type)) {
    stop("Type not provided. Use 'FS', 'HS', or 'both'.")
  }
  if (!type %in% c("HS", "FS", "both")) {
    stop("Type must be 'FS', 'HS', or 'both'.")
  }

  ped <- data.table::as.data.table(data.table::copy(pedigree))
  needed <- c(proband, mother, father, sex)
  missing_cols <- setdiff(needed, names(ped))
  if (length(missing_cols)) {
    stop("Missing pedigree columns: ", paste(missing_cols, collapse = ", "))
  }

  proband_col <- proband
  mother_col <- mother
  father_col <- father
  sex_col <- sex

  ped <- ped[, .(
    proband = .SD[[1L]],
    mother = .SD[[2L]],
    father = .SD[[3L]],
    sex = .SD[[4L]]
  ), .SDcols = c(proband_col, mother_col, father_col, sex_col)]
  ped <- ped[!is.na(proband)]
  ped[, proband_chr := as.character(proband)]

  make_pair_table <- function(dt, rel = "FS") {
    if (!nrow(dt)) {
      return(data.table::data.table(
        proband_FAM1 = character(),
        proband_FAM2 = character(),
        relation = character(),
        sex_FAM1 = character(),
        sex_FAM2 = character()
      ))
    }
    out <- dt[, .(
      proband_FAM1 = as.character(proband),
      proband_FAM2 = as.character(i.proband),
      relation = rel,
      sex_FAM1 = as.character(sex),
      sex_FAM2 = as.character(i.sex)
    )]
    out[, pair_key := cots_pair_id(proband_FAM1, proband_FAM2)]
    data.table::setorderv(out, c("pair_key", "proband_FAM1", "proband_FAM2"))
    out <- out[!duplicated(pair_key)]
    out[, pair_key := NULL]
    out
  }

  if (type != "FS") {
    mat <- ped[!is.na(mother)]
    data.table::setkey(mat, mother)
    mat_pairs <- mat[mat, allow.cartesian = TRUE, nomatch = 0L][proband_chr < i.proband_chr]
    mat_pairs[, relation_tmp := ifelse(father == i.father, "FS", "HS")]
    mat_pairs <- mat_pairs[relation_tmp == "HS"]
    mat_hs <- make_pair_table(mat_pairs, rel = "HS")

    pat <- ped[!is.na(father)]
    data.table::setkey(pat, father)
    pat_pairs <- pat[pat, allow.cartesian = TRUE, nomatch = 0L][proband_chr < i.proband_chr]
    pat_pairs[, relation_tmp := ifelse(mother == i.mother, "FS", "HS")]
    pat_pairs <- pat_pairs[relation_tmp == "HS"]
    pat_hs <- make_pair_table(pat_pairs, rel = "HS")

    halfs <- data.table::rbindlist(list(mat_hs, pat_hs), use.names = TRUE, fill = TRUE)
    if (nrow(halfs)) {
      halfs[, pair_key := cots_pair_id(proband_FAM1, proband_FAM2)]
      data.table::setorderv(halfs, c("pair_key", "proband_FAM1", "proband_FAM2"))
      halfs <- halfs[!duplicated(pair_key)]
      halfs[, pair_key := NULL]
    }
  }

  if (type != "HS") {
    fs <- ped[!is.na(mother) & !is.na(father)]
    data.table::setkey(fs, mother, father)
    fs_pairs <- fs[fs, allow.cartesian = TRUE, nomatch = 0L][proband_chr < i.proband_chr]
    fulls <- make_pair_table(fs_pairs, rel = "FS")
  }

  combined <- switch(
    type,
    both = data.table::rbindlist(list(fulls, halfs), use.names = TRUE, fill = TRUE),
    FS = fulls,
    HS = halfs
  )

  if (!nrow(combined)) {
    return(data.table::data.table(
      pair_id = character(),
      relation = character(),
      proband = character(),
      sex = character()
    ))
  }

  if (!opposite.sex) {
    combined <- combined[sex_FAM1 == sex_FAM2]
  }

  if (exclusive) {
    combined <- pair.reducer(combined, deterministic = TRUE)
  }

  combined[, pair_id := cots_pair_id(proband_FAM1, proband_FAM2)]
  out <- data.table::rbindlist(
    list(
      combined[, .(pair_id, relation, proband = proband_FAM1, sex = as.character(sex_FAM1))],
      combined[, .(pair_id, relation, proband = proband_FAM2, sex = as.character(sex_FAM2))]
    ),
    use.names = TRUE
  )
  data.table::setorderv(out, c("pair_id", "proband"))
  out
}
