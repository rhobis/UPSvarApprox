### ----------------------------------------------------------------------------
### Approximate Variance Estimators of class 2        --------------------------
### ( pik-values are required only for sample units ) --------------------------
### ----------------------------------------------------------------------------


#' Deville's approximate variance estimators
#'
#' Compute Deville's approximate variance estimators of class 2, which require only
#' first-order inclusion probabilities and only for sample units
#'
#' @param pik numeric vector of first-order inclusion probabilities for sample units
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance
#'

var_Deville <- function(y, pik, method) {

    ### Check input ---
    method <- match.arg( method, c("Deville1", "Deville2", "Deville3") )

    ### Compute c_k values ---
    if( identical(method, "Deville1") ){

        #Compute c_k values ---
        n <- length(y)
        ck <- (1-pik) * n/(n-1)

        ### Estimate variance ---
        ys <- pik * ( sum(ck*y/pik) / sum(ck) )
        delta2 <- (y-ys)**2
        v <- sum( ck*delta2 / (pik**2) )

    }else if( identical(method, "Deville2") ){

        #Compute c_k values ---
        ck <- (1-pik) / sum(1-pik)
        ck <-  1 - sum(ck**2)
        ck <- (1-pik) / ck

        ### Estimate variance ---
        ys <- pik * ( sum(ck*y/pik) / sum(ck) )
        delta2 <- (y-ys)**2
        v <- sum( ck*delta2 / (pik**2) )

    }else if( identical( method, "Deville3" ) ){

        #Compute a_k values ---
        d <- 1-pik
        ak <- d / sum(d)

        ### Estimate variance ---
        n   <- length(y)
        Tht <- sum( y/pik )

        v <- ( y/pik - Tht/n )^2
        v <- sum( d * v )
        v <- v / ( 1 - sum( ak**2 ) )
    }


    ### Return result ---
    return( v )

}


#' Fixed-Point approximate variance estimator
#'
#' Fixed-Point estimator for approximate variance estimation by Deville and Tillé (2005).
#' Estimator of class 2, it requires only first-order inclusion probabilities and
#' only for sample units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for sample units
#' @param maxIter a scalar indicating the maximum number of iterations for the
#' fixed-point procedure
#' @param eps tolerance value for the convergence of the fixed-point procedure
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance
#'

var_FixedPoint <- function(y, pik, maxIter=1000, eps=1e-05) {

    ### Input values---
    n <- length(y)
    d <- (1-pik)

    ### Compute c_k values ---
    necessaryCondition <- all( d/sum(d) < 0.5 )

    if(necessaryCondition){
        c0 <- d * n/(n-1)
        iter <- 0
        err <- Inf
        while(iter<maxIter & err>eps){
            ck   <- (c0**2 / sum(c0)) + d
            err  <- max( abs(ck - c0) )
            c0   <- ck
            iter <- iter + 1
        }
        if(err > eps) stop("Did not reach convergence")
    }else{
        ck <- n*d / ( (n-1)*sum(d) )
        ck <- d * (ck + 1)
    }


    ### Estimate variance ---
    yc <- y*ck
    yc <- outer( yc, yc, '*');    diag(yc) <- 0
    pp <- outer( pik, pik, '*');    diag(pp) <- 1
    v  <- ck - ( ck**2 / sum(ck) )
    v  <- sum( (y/pik)**2 * v )
    v  <- v - sum( yc/pp )/sum(ck)

    ### Return result ---
    return( v )

}




#' Approximate Variance Estimators by Brewer (2002)
#'
#' Computes an approximate variance estimate according to one  of
#' the estimators proposed by Brewer (2002).
#' This estimator belongs to class 2, it requires only first-order inclusion
#' probabilities and only for sample units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for sample units
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance

var_Brewer_class2 <- function(y, pik) {

    ### Compute ck values ---
    n <- length( y )
    ck <- (n-1) / (n-pik)

    ### Estimate variance ---
    v <- ( y/pik - sum(y/pik)/n )^2
    v <- sum( (1/ck - pik) * v)

    ### Return result ---
    return(v)
}


#' Hájek Approximate Variance Estimator
#'
#' Compute an approximate variance estimate using the approximation of
#' joint-inclusion probabilities proposed by Hájek (1964).
#' Estimator of class 2, it requires only first-order inclusion probabilities and
#' only for sample units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for sample units
#' @inheritParams approx_var_est
#'
#'
#' @return a scalar, the estimated variance
#'

var_Hajek <- function(y, pik) {

    ### Estimate variance ---
    n  <- length(y)
    d  <- 1-pik

    ck <- n * d / (n-1)

    B  <- sum(ck*y/pik) / sum(ck)
    ek <- (y/pik - B )**2

    v <- sum(ck*ek)

    ### Return result ---
    return( v )
}



#' Rosén Approximate Variance Estimator
#'
#' Compute an approximate variance estimate using the estimator
#' proposed by Rosén (1991).
#' Estimator of class 2, it requires only first-order inclusion probabilities and
#' only for sample units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for sample units
#' @inheritParams approx_var_est
#'
#'
#' @return a scalar, the estimated variance
#'

var_Rosen <- function(y, pik) {

    ### Estimate variance ---
    n <- length(y)
    d <- 1-pik

    A <- sum( y*d*log(d) / (pik**2) ) / sum( d*log(d) / pik )

    v <- sum( d * ( y/pik - A )^2 )
    v <- n/(n-1) * v

    ### Return result ---
    return( v )
}
