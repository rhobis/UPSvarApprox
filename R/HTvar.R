#' Variance of the Horvitz-Thompson estimator
#'
#' Compute or estimate the variance of the Horvitz-Thompson total estimator
#' by the Horvitz-Thompson or Sen-Yates-Grundy variance estimators.
#'
#' @param y numeric vector representing the variable object of inference
#' @param pikl matrix of second-order (joint) inclusion probabilities; the diagonal
#' must contain the first-order inclusion probabilities.
#' @param sample boolean indicating if sample values are provided and, thus, an
#' estimate of the variance is required
#' @param method string, indicating if the Horvitz-Thompson (\code{"HT"}) or the
#' Sen-Yates-Grundy (\code{"SYG"}) estimator should be computed.
#'
#'
#'
#'
#' @export
#'

HTvar <- function(y, pikl, sample=TRUE, method) {
    method <- match.arg(method,
                        c("HT", "SYG"))
}

