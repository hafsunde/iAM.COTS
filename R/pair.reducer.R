#' Reduce pair table so each proband appears at most once
#'
#' Internal deterministic greedy reducer.
#'
#' @param dta Data.frame/data.table with `proband_FAM1` and `proband_FAM2`.
#' @param deterministic Logical; if `FALSE`, row order is randomized before
#'   selection.
#'
#' @return Reduced data.table.
#'
#' @keywords internal
#' @noRd
pair.reducer <- function(dta, deterministic = TRUE) {
  dt <- data.table::as.data.table(data.table::copy(dta))
  order_cols <- intersect(c("pair_id", "proband_FAM1", "proband_FAM2"), names(dt))
  cots_select_nonoverlap(
    dt = dt,
    id_cols = c("proband_FAM1", "proband_FAM2"),
    order_cols = order_cols,
    deterministic = deterministic
  )
}
