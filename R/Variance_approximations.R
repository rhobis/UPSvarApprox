#' Approximate the Variance of the Horvitz-Thompson estimator
#'
#' Approximations of the Horvitz-Thompson variance for High-Entropy sampling designs.
#' Such methods use only first-order inclusion probabilities.
#'
#'
#' @param y numeric vector containing the values of the variable of interest, of
#' length equal to population size
#' @param pik numeric vector of first-order inclusion probabilities, of length
#' equal to population size
#' @param n a scalar indicating the sample size
#' @param method string indicating the approximation that should be used
#' @param ... two optional parameters can be modified to control the iterative
#' procedure in \code{method="FixedPoint"}: \code{maxIter} sets the maximum number
#' of iteration to perform and \code{eps} controls the convergence error
#'
#'
#' @details
#' TBD
#'
#'
#' @return
#' a scalar, the approximated variance.
#'
#'
#' @references
#' Matei, A.; Tillé, Y., 2005. Evaluation of variance approximations and estimators
#' in maximum entropy sampling with unequal probability and fixed sample size.
#' Journal of Official Statistics 21 (4), 543-570.
#'
#'
#'
#' @examples
#'
#' N <- 500; n <- 50
#'
#' set.seed(0)
#' x <- rgamma(n=N, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik  <- n * x/sum(x)
#' pikl <- outer(pik, pik, '*'); diag(pikl) <- pik
#'
#' ### Variance approximations ---
#' Var_approx(y, pik, n, method = "Hajek1")
#' Var_approx(y, pik, n, method = "Hajek1")
#' Var_approx(y, pik, n, method = "HartleyRao1")
#' Var_approx(y, pik, n, method = "HartleyRao2")
#' Var_approx(y, pik, n, method = "FixedPoint")
#'
#'
#' @export
#'


### Main function --------------------------------------------------------------
Var_approx <- function(y, pik, n, method, ...){

    method <- match.arg(method,
                        c("Hajek1",
                          "Hajek2",
                          "HartleyRao1",
                          "HartleyRao2",
                          "FixedPoint")
    )


    ### Input check and initialisation ---
    N <- length(y)

    if( length(n)>1 ){
        stop("Argument 'n' should be a scalar!")
    }else if( !is.numeric(n) | n != as.integer(n) ){
        stop( "Argument 'n' must be an integer number")
    }

    if( !(class(y) %in% c("numeric", "integer")) ){
        stop( "The argument 'y' should be a numeric vector!")
    }else if( N < 2 ){
        stop( "The 'y' vector is too short!" )
    }else if( any(y %in% c(NA, NaN, Inf)) ){
        stop( "The 'y' vector contains invalid values (NA, NaN, Inf)" )
    }
    if( any(y<0) ){
        message( "Some 'y' values are negative, continuing anyway...")
    }

    if( !identical( class(pik), "numeric" ) ){
        stop( "The argument 'pik' should be a numeric vector!")
    }else if( !identical(N, length(pik) ) ){
        stop( "The 'pik' vector must have same length as 'y'!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some values of the 'pik' vector are outside the interval [0, 1]")
    }else if( any(pik %in% c(NA, NaN, Inf)) ){
        stop( "The 'pik' vector contains invalid values (NA, NaN, Inf)" )
    }else if( !all.equal( sum(pik), n) ) stop("the sum of 'pik' values should be equal to 'n' ")


    ### Compute Variance ---
    if( identical(method, 'Hajek1') ){

        bk <- pik * (1-pik) * N / (N-1)
        ys <- pik * sum( bk*y/pik ) / sum(bk)
        V  <- sum( bk*(y-ys)^2 / (pik**2) )

    }else if( identical(method, 'Hajek2') ){
        d  <- pik*(1-pik)
        ak <- n * (1-pik) / sum( d )
        yt <- sum( y*ak )
        V  <- sum( d * (y/pik - yt/n)^2 )

    }else if( identical(method, 'HartleyRao1') ){

        p2   <- pik**2
        sp2  <- sum(p2)
        Y    <- sum(y)
        sdif <- (y/pik - Y/n)**2
        V    <- sum( pik * ( 1 - (n-1)*pik/n ) * sdif )
        V    <- V - (n-1)/(n**2) * sum( (2*pik^3  - p2*sp2/2) * sdif )
        V    <- V + 2*(n-1)/(n**3) * ( sum(pik*y) - Y/n * sp2 )^2

    }else if( identical(method, 'HartleyRao2') ){

        Y    <- sum(y)
        sdif <- (y/pik - Y/n)**2
        V    <- sum( pik * (1 - (n-1)/n * pik) * sdif )

    }else if( identical(method, 'FixedPoint') ){

        argList <- list(...)
        ifelse( is.null(argList$eps), eps <- 1e-05, eps <- argList$eps )
        ifelse( is.null(argList$maxIter),  maxIter <- 1000, maxIter <- argList$maxIter )

        ### Compute b_k values ---
        d <- pik*(1-pik)

        necessaryCondition <- all( d/sum(d) < 0.5 )

        if(necessaryCondition){
            b0 <- d * N/(N-1)
            iter <- 1
            err <- Inf
            while( iter<maxIter & err>eps ){
                bk  <- b0**2 / sum(b0) + d
                err <- max( abs(bk - b0) )
                b0  <- bk
                iter <- iter + 1
            }
            if(err > eps) stop("Did not reach convergence")
        }else{
            bk <- N*d / ( (N-1) * sum(d) )
            bk <- d * (bk + 1)
        }

        ys <- pik * sum( bk*y/pik ) / sum(bk)
        V  <- sum( bk*(y-ys)^2 / (pik**2) )

    }

    ### Return result ---
    return(V)
}
