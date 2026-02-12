#' Convert standard errors to confidence intervals
#'
#' Computes confidence intervals from point estimates and standard errors.
#' Optionally applies Fisher's z-transformation for correlations.
#'
#' @param estimates Numeric vector of point estimates.
#' @param std_errors Numeric vector of standard errors.
#' @param conf_level Confidence level for the interval. Default is 0.95.
#' @param is_correlation Logical. Set to TRUE if estimates are correlations.
#'
#' @return A data.frame with columns Lower and Upper.
#'
#' @details
#' For correlations, confidence intervals are computed on the Fisher z scale
#' and then back-transformed to the correlation scale.
#'
#' @examples
#' # Standard normal-theory confidence intervals
#' SEtoCI(
#'   estimates = c(0.5, 1.2),
#'   std_errors = c(0.1, 0.2)
#' )
#'
#' # Confidence intervals for correlations using Fisher z
#' SEtoCI(
#'   estimates = c(0.2, 0.4),
#'   std_errors = c(0.05, 0.05),
#'   is_correlation = TRUE
#' )

#' @export
SEtoCI <- function(estimates, std_errors, conf_level = 0.95, is_correlation = FALSE) {

  # Defensive checks to fail early and clearly
  stopifnot(
    is.numeric(estimates),
    is.numeric(std_errors),
    length(estimates) == length(std_errors),
    conf_level > 0 && conf_level < 1
  )

  # Compute the normal-theory critical value
  z_value <- stats::qnorm((1 + conf_level) / 2)

  if (is_correlation) {

    # Fisher z-transformation stabilizes variance of correlations
    z_estimates <- atanh(estimates)

    # Compute bounds on the z scale
    lower_bound_z <- z_estimates - z_value * std_errors
    upper_bound_z <- z_estimates + z_value * std_errors

    # Back-transform to correlation scale
    lower_bound <- tanh(lower_bound_z)
    upper_bound <- tanh(upper_bound_z)

  } else {

    # Standard normal-theory confidence interval
    lower_bound <- estimates - z_value * std_errors
    upper_bound <- estimates + z_value * std_errors
  }

  # Return as a data frame for downstream compatibility
  data.frame(
    Lower = lower_bound,
    Upper = upper_bound
  )
}
