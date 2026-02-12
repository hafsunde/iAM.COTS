#' Extract variance component estimates from an iAM-COTS fit
#'
#' Extracts the variance-component table (est_var) from a fitted OpenMx iAM-COTS
#' model object and returns it as a data frame. The function uses model-based
#' confidence intervals if available. Otherwise, it computes Wald-type confidence
#' intervals from standard errors, or optionally uses bootstrap quantiles if the
#' model was fitted with MxComputeBootstrap.
#'
#' @param fit A fitted OpenMx model (typically returned by iam.cots.fun)
#'   containing an est_var algebra and optionally confidence intervals in
#'   fit$output$confidenceIntervals.
#' @param sim_vars Optional object containing simulated variance-component results
#'   to merge onto the returned table. It is expected to be coercible to a vector
#'   via t(sim_vars) with names corresponding to est_var column names.
#' @param boot.CI Logical. If TRUE and fit$compute is an MxComputeBootstrap
#'   object, bootstrap quantile bounds (0.025 and 0.975) are used for
#'   lbound and ubound. Otherwise Wald-type confidence intervals are used.
#'
#' @return A data frame with columns:
#'   - parameter
#'   - estimate
#'   - lbound
#'   - ubound
#'   If sim_vars is provided, a column named "simulated" is added.
#'
#' @details
#' The function proceeds in the following order:
#' 1. If fit$output$confidenceIntervals is present, bounds are taken from that
#'    object (assumed to align with the first 29 est_var parameters as used in
#'    the iAM-COTS model).
#' 2. If the model used MxComputeBootstrap and boot.CI = TRUE, bootstrap
#'    quantiles are used via mxBootstrapEvalByName("est_var", ...).
#' 3. Otherwise, standard errors are computed using mxSE() and converted to
#'    confidence intervals using SEtoCI().
#'
#' This function assumes the helper function SEtoCI() is available in the
#' package namespace.
#'
#' @seealso \code{\link[OpenMx]{mxSE}}, \code{\link[OpenMx]{mxBootstrapEvalByName}}
#'
#' @export


get.varcomps <- function(fit, sim_vars = NULL, boot.CI = FALSE) {


  results <- data.frame(
    parameter = colnames(fit$est_var$result),
    estimate = unname(t(fit$est_var$result)),
    lbound = NA_real_,
    ubound = NA_real_

  )
  if (!is.null(fit$output$confidenceIntervals)) {
    results[,3:4] <- fit$output$confidenceIntervals[1:29,c(1,3)]
  } else if (is(fit$compute, "MxComputeBootstrap") & boot.CI) {
    message("Calculating CIs from bootstrap samples")
    results[,3:4] <- mxBootstrapEvalByName("est_var", fit, bq=c(.025, .975))[,2:3]
  } else {
    results[,3:4] <- SEtoCI(results$estimate, as.vector(mxSE(est_var, fit)))
  }

  if (!is.null(sim_vars)) {
    simulated <- t(sim_vars) %>%
      as.data.frame() %>%
      mutate(parameter = rownames(.)) %>%
      rename(simulated = V1)
    results <- results %>%
      left_join(simulated, by = "parameter")
  }


  return(results)
}

