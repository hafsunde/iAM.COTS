#' Extract Wald SE and confidence intervals from an OpenMx model object
#'
#' Convenience wrapper to extract estimates and SEs from an OpenMx model slot
#' and compute Wald confidence intervals via `HFSutils::SEtoCI()`.
#'
#' @param fit An OpenMx model fit object.
#' @param object Character. Name of the slot inside `fit` to extract, e.g. "A", "S", etc.
#' @param is_correlation Logical. If TRUE, compute CI using Fisher z transform in
#'   `HFSutils::SEtoCI()`.
#'
#' @return A data.frame with columns type, estimate, parameter, SE, lbound, ubound.
#'
#' @details
#' This function assumes `OpenMx::mxSE()` works for the requested `object`.
#'
#' @export
mx.WaldCI <- function(fit, object, is_correlation = FALSE) {

  se_mat <- mxSE(object, fit)

  results <- combine_SE(fit[[object]]$result, se_mat)
  results <- as.data.frame(results)

  names(results)[names(results) == "Row"] <- "type"
  names(results)[names(results) == "Estimate"] <- "estimate"
  names(results)[names(results) == "Col"] <- "parameter"

  results$lbound <- NA_real_
  results$ubound <- NA_real_

  results[, c("lbound", "ubound")] <- HFSutils::SEtoCI(
    estimates = results$estimate,
    std_errors = results$SE,
    is_correlation = is_correlation
  )

  results
}

#' Combine an estimate matrix and an SE matrix into long format
#'
#' Converts a matrix of estimates and a matrix of standard errors of the same
#' dimension into a long data.frame with one row per cell.
#'
#' @param M A numeric matrix of estimates.
#' @param SE A numeric matrix of standard errors. Must have the same dimensions
#'   as M.
#'
#' @return A data.frame with columns Row, Col, Estimate, and SE.
#'
#' @export
combine_SE <- function(M, SE) {
  row_names <- rownames(M)
  col_names <- colnames(M)

  if (is.null(row_names)) {
    row_names <- as.character(1:nrow(M))
  }
  if (is.null(col_names)) {
    col_names <- as.character(1:ncol(M))
  }

  data.frame(
    Row = rep(row_names, times = length(col_names)),
    Col = rep(col_names, each = length(row_names)),
    Estimate = as.vector(M),
    SE = as.vector(SE)
  )
}

#' Combine an estimate matrix with an OpenMx confidence interval table
#'
#' Reshapes a matrix of estimates into long format and merges it with confidence
#' interval bounds extracted from an OpenMx confidence interval object.
#'
#' @param M A numeric matrix of estimates (for example fit$est_PO$result).
#' @param CI A confidence interval table, typically a subset of
#'   fit$output$confidenceIntervals, with rownames containing matrix indices in
#'   the form "[i,j]" and columns named "lbound" and "ubound".
#'
#' @return A data.frame with columns Row, Col, Estimate, lbound, and ubound.
#'
#' @export
combine_CI <- function(M, CI) {
  row_names <- rownames(M)
  col_names <- colnames(M)

  if (is.null(row_names)) {
    row_names <- as.character(1:nrow(M))
  }
  if (is.null(col_names)) {
    col_names <- as.character(1:ncol(M))
  }

  df_long <- data.frame(
    Row = rep(row_names, times = length(col_names)),
    Col = rep(col_names, each = length(row_names)),
    Estimate = as.vector(M)
  )

  rn <- rownames(CI)
  row_indices <- as.numeric(sub(".*\\[(\\d+),(\\d+)\\]", "\\1", rn))
  col_indices <- as.numeric(sub(".*\\[(\\d+),(\\d+)\\]", "\\2", rn))

  CI_df <- data.frame(
    Row = row_names[row_indices],
    Col = col_names[col_indices],
    lbound = CI[, "lbound"],
    ubound = CI[, "ubound"]
  )

  merge(df_long, CI_df, by = c("Row", "Col"), all.x = TRUE, sort = FALSE)
}
