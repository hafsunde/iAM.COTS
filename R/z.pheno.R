#' Pooled z-standardization across selected phenotype columns
#'
#' Computes one pooled mean and SD across all selected variables, then applies
#' the same standardization to each listed column.
#'
#' @param dt Data.frame/data.table.
#' @param vars Character vector of column names to standardize.
#'
#' @return A data.table with standardized `vars`.
#'
#' @export
z.pheno <- function(dt, vars) {
  x <- data.table::as.data.table(data.table::copy(dt))
  missing_cols <- setdiff(vars, names(x))
  if (length(missing_cols)) {
    stop("Missing columns in `dt`: ", paste(missing_cols, collapse = ", "))
  }

  pooled <- unlist(x[, ..vars], use.names = FALSE)
  pooled_mean <- mean(pooled, na.rm = TRUE)
  pooled_sd <- stats::sd(pooled, na.rm = TRUE)

  if (is.na(pooled_sd) || pooled_sd == 0) {
    warning("Pooled SD is zero or NA; returning NA-standardized values.")
  }

  x[, (vars) := lapply(.SD, function(v) (v - pooled_mean) / pooled_sd), .SDcols = vars]
  x
}
