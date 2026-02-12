#' Extract parent-offspring correlation decomposition from an iAM-COTS fit
#'
#' Extracts the parent-offspring correlation decomposition (est_PO) from a fitted
#' OpenMx iAM-COTS model object and returns it as a data frame. The function
#' uses model-based confidence intervals if available. Otherwise, it computes
#' Wald-type confidence intervals from standard errors, or optionally uses
#' bootstrap quantiles if the model was fitted with MxComputeBootstrap.
#'
#' @param fit A fitted OpenMx model (typically returned by iam.cots.fun)
#'   containing an est_PO algebra and optionally confidence intervals in
#'   fit$output$confidenceIntervals.
#' @param sim_PO Optional data frame of simulated parent-offspring results.
#'   Must contain columns "parameter", "type", and "estimate". The estimate
#'   column will be renamed to "simulated".
#' @param boot.CI Logical. If TRUE and fit$compute is an MxComputeBootstrap
#'   object, bootstrap quantile bounds (0.025 and 0.975) are used for
#'   lbound and ubound. Otherwise Wald-type confidence intervals are used.
#'
#' @return A data frame with columns:
#'   - type
#'   - parameter
#'   - estimate
#'   - SE (if applicable)
#'   - lbound
#'   - ubound
#'   If sim_PO is provided, a column named "simulated" is added.
#'
#' @details
#' The function proceeds in the following order:
#' 1. If fit$output$confidenceIntervals contains rows matching "est_PO",
#'    those confidence intervals are returned.
#' 2. If the model used MxComputeBootstrap and boot.CI = TRUE, bootstrap
#'    quantiles are used.
#' 3. Otherwise, standard errors are computed using mxSE() and converted
#'    to confidence intervals using SEtoCI().
#'
#' This function assumes the helper functions combine_CI(), combine_SE(),
#' and SEtoCI() are available in the package namespace.
#'
#' @export


get.POcor <- function(fit, sim_PO = NULL, boot.CI = FALSE) {


  if (!is.null(fit$output$confidenceIntervals)) {
    results <- combine_CI(fit$est_PO$result, fit$output$confidenceIntervals[grepl("est_PO", rownames(fit$output$confidenceIntervals)),])  %>%
      as.data.frame() %>%
      rename(type = Row, estimate = Estimate, parameter = Col)
  } else if (is(fit$compute, "MxComputeBootstrap") & boot.CI) {
    results <-   combine_SE(fit$est_PO$result, mxSE(est_PO, fit)) %>%
      as.data.frame() %>%
      rename(type = Row, estimate = Estimate, parameter = Col) %>%
      mutate(lbound=NA_real_,ubound=NA_real_)
    message("Calculating CIs from bootstrap samples")
    results_PO.boot  <- mxBootstrapEvalByName("est_PO", fit, bq=c(.025, .975))
    firstindex <-  c(1,4,7,10)
    results$lbound <- results_PO.boot[c(firstindex,firstindex+1,firstindex+2), 2 ]
    results$ubound <- results_PO.boot[c(firstindex,firstindex+1,firstindex+2), 3 ]
  } else {
    results <-   combine_SE(fit$est_PO$result, mxSE(est_PO, fit)) %>%
      as.data.frame() %>%
      rename(type = Row, estimate = Estimate, parameter = Col) %>%
      mutate(lbound=NA_real_,ubound=NA_real_)
    results[,c("lbound", "ubound")] <- SEtoCI(results$estimate, results$SE)
  }

  if (!is.null(sim_PO)) {
    results <- results %>%
      left_join(rename(sim_PO, simulated = estimate), by = c("parameter", "type"))
  }


  return(results)
}



