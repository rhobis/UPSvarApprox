#' UPSvarApprox: Approximate the variance of the Horvitz-thompson estimator
#'
#'
#' @section Variance approximation:
#' The package provides function \code{Var_approx} for the approximation of the
#' Horvitz-Thompson variance and function \code{approx_var_est}, for the computation
#' of an approximate estimate of the variance.
#' Different methods are available for both functions, see their documentation for
#' more details.
#'
#' @section Approximation of Joint-inclusion probabilities:
#' Function \code{jip_approx} provides a number of approximations of the
#' second-order inclusion probabilities that require only the first-order inclusion
#' probabilities.
#' The variance of the Horvitz-Thompson total estimator may be then estimated by
#' plugging the approximated joint probabilities into the Horvitz-Thompson or
#' Sen-Yates-Grundy variance estimator using function \code{HTvar}.
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
