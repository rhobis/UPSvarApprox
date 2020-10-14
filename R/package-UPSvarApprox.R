#' UPSvarApprox: Approximate the variance of the Horvitz-Thompson estimator
#'
#'
#' @description
#' Variance approximations for the
#' Horvitz-Thompson total estimator in Unequal Probability Sampling
#' using only first-order inclusion probabilities.
#' See Matei and Tillé (2005) and Haziza, Mecatti and Rao (2008) for details.
#'
#' @section Variance approximation:
#' The package provides function \code{\link{Var_approx}} for the approximation of the
#' Horvitz-Thompson variance, and function \code{\link{approx_var_est}} for the computation
#' of approximate variance estimators.
#' For both functions, different estimators are implemented,
#' see their documentation for details.
#'
#'
#'
#' @docType package
#'
#' @name UPSvarApprox
#'
#' @references
#'
#' Matei, A.; Tillé, Y., 2005. Evaluation of variance approximations and estimators
#' in maximum entropy sampling with unequal probability and fixed sample size.
#' Journal of Official Statistics 21 (4), 543-570.
#'
#' Haziza, D.; Mecatti, F.; Rao, J.N.K. 2008.
#' Evaluation of some approximate variance estimators under the Rao-Sampford
#' unequal probability sampling design. Metron LXVI (1), 91-108.
#'
#'
NULL
