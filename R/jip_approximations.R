#' Approximate Joint-Inclusion Probabilities
#'
#' Approximations of joint-inclusion probabilities by means of first-order
#' inclusion probabilities.
#'
#' @param pik numeric, vector of first-order inclusion probabilities for all
#' population units.
#' @param method string representing one of the available approximation methods.
#'
#'
#' @details
#' Available methods are "Hajek", "HartleyRao", "Tille" and "Brewer18"
#'
#'
#' @return A symmetric matrix of inclusion probabilities, which diagonal is the
#' vector of first-order inclusion probabilities.
#'
#'
#' @references
#'
#' Hartley, H.O.; Rao, J.N.K., 1962. Sampling With Unequal Probability and Without Replacement.
#' The Annals of Mathematical Statistics 33 (2), 350-374.
#'
#' Hájek, J., 1964. Asymptotic Theory of Rejective Sampling with Varying Probabilities from a Finite Population.
#' The Annals of Mathematical Statistics 35 (4), 1491-1523.
#'
#' Tillé, Y., 1996. Some Remarks on Unequal Probability Sampling Designs Without Replacement.
#' Annals of Economics and Statistics 44, 177-189.
#'
#' Brewer, K.R.W.; Donadio, M.E., 2003. The High Entropy Variance of the Horvitz-Thompson Estimator.
#' Survey Methodology 29 (2), 189-196.
#'
#' @export



jip_approx <- function( pik, method ){

    ### Check input ---
    method <- match.arg(method, c( "Hajek", "HartleyRao", "Tille", "Brewer18") )

    if( !identical( class(pik), "numeric" ) ){
        stop( "Argument 'pik' should be a numeric vector!")
    }else if( length(pik) < 2 ){
        stop( "The 'pik' vector is too short!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some values of the 'pik' vector are outside the interval [0, 1]")
    }else if( !all.equal(sum(pik), round(sum(pik)) ) ){
        stop( "The sum of 'pik' values is not an integer!")
    }

    ### Call method ---
    jips <- switch(method,
           "Hajek"       = jip_Hajek(pik),
           "HartleyRao"  = jip_HartleyRao(pik),
           "Tille"       = jip_Tille(pik),
           "Brewer18"    = jip_Brewer(pik, method)
           )

    ### Return result ---
    return( jips )
}



#' Brewer's joint-inclusion probability approximations
#'
#' Approximation of joint inclusion probabilities by one of the estimators
#' proposed by Brewer and Donadio (2003)
#'
#' @inheritParams jip_approx
#'
#' @details
#' \code{"Brewer18"} is the approximation showed in equation (18) of Brewer and Donadio (2003)
#'

jip_Brewer <- function(pik, method){

    method <- match.arg(method, c("Brewer18"))
    n <- sum(pik)

    ### Estimate jips ---
    ci <- (n-1) / (n-(2*n-1)/(n-1) * pik + sum(pik**2)/(n-1))
    cc <- outer(ci,ci, '+')
    pp <- outer(pik,pik, '*')
    out <- cc*pp / 2
    diag(out) <- pik

    ### Return result ---
    return(out)
}




#' Hájek's joint-inclusion probability approximation
#'
#' Estimate joint-inclusion probabilities using Hájek (1964) equation
#'
#'
#' @inheritParams jip_approx
#'
#'

jip_Hajek <- function(pik){

    ### Estimate jips ---
    d   <- sum(pik*(1-pik))
    out <- outer(pik,pik,'*') * (1 - outer(1-pik, 1-pik, '*') / d)
    diag(out) <- pik

    ### Return result ---
    return(out)
}




#' Hartley-Rao approximation of joint-inclusion probabilities
#'
#' Approximation of joint-inclusion probabilities with precision of order
#' \eqn{O(N^{-4})} for the random systematic sampling design, by Hartley and Rao (1962).
#'
#' @inheritParams jip_approx
#'

jip_HartleyRao <- function(pik){

    ### Estimate jips ---
    pp <- outer(pik, pik, '*')
    n  <- sum(pik)
    nn <- (n-1)/n

    p2 <- pik**2
    p3 <- pik**3
    sp2 <- sum(p2)
    sp3 <- sum(p3)

    out <-  nn*pp +
        nn/n * (outer(p2,pik) + outer(pik, p2)) -
        nn/(n**2) * pp * sp2 +
        (2*nn)/(n**2) * (outer(p3,pik) + outer(pik,p3) + outer(p2,p2)) -
        3*nn/n**3 * (outer(p2,pik) + outer(pik,p2)) * sp2 +
        3*nn/n**4 * pp * sp2**2 -
        2*nn/n**3 * pp * sp3

    diag(out) <- pik

    ### Return result ---
    return(out)
}



#' Tillé's approximation of joint-inclusion probabilities
#'
#' Compute the approximation of joint-inclusion probabilities by means of the
#' Iterative Proportional Fitting Procedure (IPFP) proposed by Tillé (1996)
#'
#' @param maxIter a scalar indicating the maximum number of iterations for the
#' fixed-point procedure
#' @param eps tolerance value for the convergence of the fixed-point procedure
#' @inheritParams jip_approx
#'

jip_Tille <- function(pik,eps=1e-06, maxIter=1000){

    n <- sum(pik)

    ### Compute b_k values ---
    b0   <- pik
    iter <- 0
    err  <- Inf
    while( iter<maxIter & err>eps ){
        b1 <- (n-1)*pik / ( sum(b0)-b0 )
        b2 <- b1 * sqrt( n*(n-1) / ( sum(b1)**2 - sum(b1**2) ) )
        b0 <- b2
        tab <- outer(b2,b2,'*');    diag(tab) <- 0
        err <- max( abs( rowSums(tab) - pik*(n-1) ) )
        iter <- iter + 1
    }

    ### Estimate jips ---
    out <- outer(b2,b2,'*')
    diag(out) <- pik

    ### Return result ---
    return(out)
}


