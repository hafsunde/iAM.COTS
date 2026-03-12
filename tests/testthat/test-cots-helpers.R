make_spine_inputs <- function() {
  cohort <- data.frame(
    w19_1011_lnr_k2_ = paste0("ch", 1:12),
    mor_lnr_k2_ = c("m1", "m1", "m2", "m2", "m3", "m3", "m4", "m4", "m6", "m6", "m7", "m7"),
    far_lnr_k2_ = c("f1", "f1", "f2", "f2", "f3", "f3", "f4", "f4", "f5", "f5", "f6", "f6"),
    nuc_fam_id = c("nf1", "nf1", "nf2", "nf2", "nf3", "nf3", "nf4", "nf4", "nf5", "nf5", "nf6", "nf6"),
    stringsAsFactors = FALSE
  )

  ancestry <- data.frame(
    proband = c("m1", "m2", "m3", "m4", "m6", "m7", "f1", "f2", "f3", "f4", "f5", "f6"),
    mother = c("am1", "am1", "am2", "am2", "am6", "am7", "af1", "af1", "af2", "af2", "af5", "af6"),
    father = c("bm1", "bm2", "bm3", "bm3", "bf6", "bf7", "bf1", "bf1", "bf2", "bf3", "bf5", "bf6"),
    sex = c(2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  twins <- data.frame(
    proband = c("m1", "m5", "f3", "f4", "m6", "m7"),
    tw_size = c(2, 2, 2, 2, 2, 2),
    SameSex = c(1, 1, 1, 1, 1, 1),
    Zygo = c(1, 1, 2, 2, NA, NA),
    pair_id = c("T1", "T1", "T2", "T2", "T3", "T3"),
    sex = c(2, 2, 1, 1, 2, 2),
    stringsAsFactors = FALSE
  )

  list(cohort = cohort, ancestry = ancestry, twins = twins)
}

test_that("cots.spine.fun is deterministic by default", {
  inp <- make_spine_inputs()
  s1 <- cots.spine.fun(inp$cohort, inp$ancestry, inp$twins, include.UZ = TRUE)
  s2 <- cots.spine.fun(inp$cohort, inp$ancestry, inp$twins, include.UZ = TRUE)
  expect_identical(as.data.frame(s1), as.data.frame(s2))
})

test_that("cots.spine.fun excludes UZ by default and includes UZ when requested", {
  inp <- make_spine_inputs()
  no_uz <- cots.spine.fun(inp$cohort, inp$ancestry, inp$twins, include.UZ = FALSE, include.halfsibs = FALSE)
  with_uz <- cots.spine.fun(inp$cohort, inp$ancestry, inp$twins, include.UZ = TRUE, include.halfsibs = FALSE)

  expect_false(any(grepl("^UZ_|^Single_UZ_", no_uz$Gruppe)))
  expect_true(any(grepl("^UZ_|^Single_UZ_", with_uz$Gruppe)))
})

test_that("selected spine has no parent/spouse overlap across families", {
  inp <- make_spine_inputs()
  out <- cots.spine.fun(inp$cohort, inp$ancestry, inp$twins, include.UZ = TRUE)
  out_df <- as.data.frame(out)
  ids <- unlist(out_df[, c("FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F")], use.names = FALSE)
  ids <- ids[!is.na(ids)]
  expect_equal(length(ids), length(unique(ids)))
})

test_that("HS is prioritized before FS when include.halfsibs is TRUE", {
  cohort <- data.frame(
    w19_1011_lnr_k2_ = paste0("ch", 1:8),
    mor_lnr_k2_ = c("m1", "m1", "m2", "m2", "m3", "m3", "m4", "m4"),
    far_lnr_k2_ = c("f1", "f1", "f2", "f2", "f3", "f3", "f4", "f4"),
    nuc_fam_id = c("n1", "n1", "n2", "n2", "n3", "n3", "n4", "n4"),
    stringsAsFactors = FALSE
  )

  ancestry <- data.frame(
    proband = c("f1", "f2", "f3", "f4"),
    mother = c("mx", "mx", "mx", "my"),
    father = c("fa", "fa", "fb", "fc"),
    sex = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  twins <- data.frame(
    proband = character(),
    tw_size = integer(),
    SameSex = integer(),
    Zygo = integer(),
    pair_id = character(),
    sex = integer(),
    stringsAsFactors = FALSE
  )

  out <- cots.spine.fun(
    cohort = cohort,
    ancestry = ancestry,
    twins = twins,
    include.halfsibs = TRUE,
    include.UZ = FALSE
  )

  out_df <- as.data.frame(out)
  hs <- out_df[grepl("^HS_", out_df$Gruppe), , drop = FALSE]
  fs <- out_df[grepl("^FS_", out_df$Gruppe), , drop = FALSE]
  expect_true(nrow(hs) >= 1)

  if (nrow(fs) > 0) {
    hs_ids <- unique(unlist(as.data.frame(hs)[, c("FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F")], use.names = FALSE))
    fs_ids <- unique(unlist(as.data.frame(fs)[, c("FAM1_M", "FAM1_F", "FAM2_M", "FAM2_F")], use.names = FALSE))
    hs_ids <- hs_ids[!is.na(hs_ids)]
    fs_ids <- fs_ids[!is.na(fs_ids)]
    expect_false(any(fs_ids %in% hs_ids))
  }
})

test_that("find.sibs uses full-ID canonical pair keys (no substring collision)", {
  ped <- data.frame(
    proband = c("AAA0000001X", "CCC0000002X", "BBB0000001Y", "DDD0000002Y"),
    mother = c("M1", "M1", "M2", "M2"),
    father = c("F1", "F2", "F3", "F4"),
    sex = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  out <- find.sibs(ped, type = "HS", opposite.sex = TRUE)
  expect_equal(length(unique(out$pair_id)), 2L)
})

test_that("join.phenotypes available-case retains information vs complete-case", {
  spine <- data.frame(
    FAM1_M = "m1",
    FAM1_F = "f1",
    FAM2_M = "m2",
    FAM2_F = "f2",
    FAM1_CH1 = "c1",
    FAM1_CH2 = "c2",
    FAM2_CH1 = "c3",
    FAM2_CH2 = "c4",
    Gruppe = "FS_ff_22",
    stringsAsFactors = FALSE
  )

  pheno <- data.frame(
    id = c("c1", "c2", "c3", "c4"),
    ch = c(1, 2, 3, 4),
    ma = c(NA, 11, 12, 13),
    fa = c(21, 22, 23, 24),
    stringsAsFactors = FALSE
  )

  avail <- join.phenotypes(
    spine_dt = spine,
    pheno_dt = pheno,
    ph.ch = "ch",
    ph.ma = "ma",
    ph.pa = "fa",
    proband = "id",
    extend = FALSE,
    missing_policy = "available"
  )

  comp <- join.phenotypes(
    spine_dt = spine,
    pheno_dt = pheno,
    ph.ch = "ch",
    ph.ma = "ma",
    ph.pa = "fa",
    proband = "id",
    extend = FALSE,
    missing_policy = "complete"
  )

  expect_false(is.na(avail$FAM1_CH1))
  expect_true(is.na(comp$FAM1_CH1))
})
