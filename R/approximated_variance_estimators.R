#' Approximated Variance Estimators
#'
#' Approximated variance estimators which use only first-order inclusion probabilities
#'
#' @param y numeric vector of sample observations
#' @param pik numeric vector of first-order inclusion probabilities, of length *N* or *n*,
#' depending on the chosen method (see Details for more information)
#' @param method string indicating the desired approximate variance estimator
#' @param sample Either a numeric vector of length equal to the sample size, containing
#' the indeces of sample units, or a boolean vector of the same length of \code{pik}, indicating which
#' units belong to the sample (\code{TRUE} if the unit is in the sample,
#' \code{FALSE} otherwise.
#' Only used with estimators of the third class (see Details for more information).
#' @param ... two optional parameters can be modified to control the iterative
#' procedures used \code{method} is \code{"MateiTille5"}, \code{"Tille"} or
#' \code{"FixedPoint"}: \code{maxIter} sets the maximum number
#' of iteration to perform and \code{eps} controls the convergence error
#'
#' @details
#'
#' Matei and Tillé (2005) define three classes:
#'
#' 1) first and second-order inclusion probabilities
#' The first class is composed of the Horvitz-Thompson estimator (Horvitz and Thompson 1952) and the Sen-Yates-Grundy estimator
#' (Yates and Grundy 1953; Sen 1953) and is not available in this package.
#' For this estimators, the use may refer to package \code{sampling} ...
#'
#' 2) only first-order inclusion probabilities, only for sample units
#'
#' >>> TBD: add estimators names and formulas
#'
#' 3) only first-order inclusion probabilities, for the entire population
#'
#' >>> TBD: add estimators names and formulas
#'
#'
#'
#' @return a scalar, the estimated variance
#'
#' @examples
#'
#' ### Generate population data ---
#' N <- 500; n <- 50
#'
#' set.seed(0)
#' x <- rgamma(500, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik <- n * x/sum(x)
#' s   <- sample(N, n)
#'
#' ys <- y[s]
#' piks <- pik[s]
#'
#' ### Estimators of class 2 ---
#' approx_var_est(ys, piks, method="Deville1")
#' approx_var_est(ys, piks, method="Deville2")
#' approx_var_est(ys, piks, method="Deville3")
#' approx_var_est(ys, piks, method="Rosen")
#' approx_var_est(ys, piks, method="FixedPoint")
#' approx_var_est(ys, piks, method="Brewer2")
#'
#' ### Estimators of class 3 ---
#' approx_var_est(ys, pik, method="Berger", sample=s)
#' approx_var_est(ys, pik, method="Tille", sample=s)
#' approx_var_est(ys, pik, method="MateiTille1", sample=s)
#' approx_var_est(ys, pik, method="MateiTille2", sample=s)
#' approx_var_est(ys, pik, method="MateiTille3", sample=s)
#' approx_var_est(ys, pik, method="MateiTille4", sample=s)
#' approx_var_est(ys, pik, method="MateiTille5", sample=s)
#' approx_var_est(ys, pik, method="Brewer1", sample=s)
#' approx_var_est(ys, pik, method="Brewer3", sample=s)
#' approx_var_est(ys, pik, method="Brewer4", sample=s)
#'
#'
#' @references
#' Matei, A.; Tillé, Y., 2005. Evaluation of variance approximations and estimators
#' in maximum entropy sampling with unequal probability and fixed sample size.
#' Journal of Official Statistics 21 (4), 543-570.
#'
#' Haziza, D.; Mecatti, F.; Rao, J.N.K. 2008.
#' Evaluation of some approximate variance estimators under the Rao-Sampford
#' unequal probability sampling design. Metron LXVI (1), 91-108.
#'
#'
#' @export
#'
#'

approx_var_est <- function(y, pik, method, sample=NULL, ...){

    ### Check input ---
    method <- match.arg( method,
                         c(  #second class
                             "Deville1",
                             "Deville2",
                             "Deville3",
                             "Rosen",
                             "FixedPoint",
                             "Brewer2",
                             #third class
                             "Berger",
                             "Tille",
                             "MateiTille1",
                             "MateiTille2",
                             "MateiTille3",
                             "MateiTille4",
                             "MateiTille5",
                             "Brewer1",
                             "Brewer3",
                             "Brewer4" )
    )

    argList <- list(...)
    ifelse( is.null(argList$eps), eps <- 1e-05, eps <- argList$eps )
    ifelse( is.null(argList$maxIter),  maxIter <- 1000, maxIter <- argList$maxIter )


    ly <- length(y)
    lp <- length(pik)
    ls <- length(sample)

    if( !identical( class(pik), "numeric" ) ){
        stop( "The argument 'pik' should be a numeric vector!")
    }else if( lp < 2 ){
        stop( "The 'pik' vector is too short!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some 'pik' values are outside the interval [0, 1]")
    }else if( any(pik %in% c(NA, NaN, Inf)) ){
        stop( "The 'pik' vector contains invalid values (NA, NaN, Inf)" )
    }

    if( !identical( class(y), "numeric" ) ){
        stop( "The argument 'y' should be a numeric vector!")
    }else if( ly < 2 ){
        stop( "The 'y' vector is too short!" )
    }else if( any(y %in% c(NA, NaN, Inf)) ){
        stop( "The 'y' vector contains invalid values (NA, NaN, Inf)" )
    }


    if( any(y<0) ){
        message( "Some 'y' values are negative, continuing anyway...")
    }


    #second class of estimators
    if( method %in% c( "Deville1", "Deville2", "Deville3", "Rosen",
                       "FixedPoint", "Brewer2") ){
        if( !identical( length(y), length(pik) ) )
            stop( "Method ", "'", method, "' ",
                  "requires argument 'y' and 'pik' to have the same length, ",
                  "which should be equal to the sample size!"
            )

    }

    #third class of estimators
    if( method %in% c( "Berger", "Tille", "MateiTille1", "MateiTille2",
                       "MateiTille3", "MateiTille4", "MateiTille5",
                       "Brewer1", "Brewer3", "Brewer4" ) ){

        if( sum(pik) != as.integer(sum(pik)) )
            stop( "The sum of pik values is not an integer!")

        if( ly >= lp )
            stop("With method ", "'", method, "'",
                 ", pik values for all population units should be provided. ",
                 "However, the length of pik is <= the length of y, ",
                 "which should be equal to the sample size. ",
                 "Please, check again your input or check the help page for more details!"
            )

        if( missing(sample) )
            stop("The argument 'sample' is required for method=", "'", method, "'")
        if( !is.numeric(sample) & !is.logical(sample) )
            stop( "The argument 'sample' should be either a numeric or logical vector")
        if( is.logical(sample) & !identical( length(pik), length(sample) ) )
            stop( "The arguments 'pik' and 'sample' should have the same length when 'sample' is a logical vector!" )
        if( is.numeric(sample) & !identical( length(y), length(sample) ) )
            stop( "The argument 'sample' should have length equal to the sample size!" )


    }



    ### Call method ---
    v <- switch(method,
                #second class
                "Deville1"    = var_Deville(y, pik, method),
                "Deville2"    = var_Deville(y, pik, method),
                "Deville3"    = var_Deville(y, pik, method),
                "Rosen"       = var_Rosen(y, pik),
                "FixedPoint"  = var_FixedPoint(y, pik, eps=1e-5, maxIter=1000),
                "Brewer2"     = var_Brewer_class2(y, pik),
                #third class
                "Berger"      = var_Berger(y, pik, method, sample),
                "Tille"       = var_Tille(y, pik, sample, eps=1e-5, maxIter=1000),
                "MateiTille1" = var_MateiTille(y, pik, method, sample) ,
                "MateiTille2" = var_MateiTille(y, pik, method, sample),
                "MateiTille3" = var_MateiTille(y, pik, method, sample),
                "MateiTille4" = var_MateiTille(y, pik, method, sample),
                "MateiTille5" = var_MateiTille(y, pik, method, sample),
                "Brewer1"     = var_Brewer_class3(y, pik, method, sample),
                "Brewer3"     = var_Brewer_class3(y, pik, method, sample),
                "Brewer4"     = var_Brewer_class3(y, pik, method, sample)
    )

    # do.call( FUN, list(y, pik))

    ### Return result ---
    return( v )
}



