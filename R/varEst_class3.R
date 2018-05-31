### ----------------------------------------------------------------------------
### Approximate Variance Estimators of class 3                  ----------------
### ( pik-values are required for all units in the population ) ----------------
### ----------------------------------------------------------------------------




#' Berger approximate variance estimator
#'
#' Compute the Berger approximate variance estimator (Berger, 1998)
#' Estimator of class 3, it requires only first-order inclusion probabilities but
#' for all population units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for all population units
#' @inheritParams approx_var_est
#'
#'
#' @return a scalar, the estimated variance
#'

var_Berger <- function(y, pik, sample) {

    ### Input values ---
    n <- length(y)

    ### Compute c_k values ---
    d  <- 1-pik
    ck <- d * n / (n-1)
    ck <- ck * sum(d[sample]) / sum(pik*d)
    ck <- ck[sample]

    ### Estimate variance ---
    piks <- pik[sample]
    ys <- piks * sum(ck*y/piks) / sum(ck)
    delta2 <- (y-ys)**2
    v <- sum( ck * delta2 / (piks**2) )

    ### Return result ---
    return( v )

}

#' Hartley and Rao approximate variance estimator
#'
#' Compute the an approximate variance estimator obtained by the approximation
#' of joint-inclusion probabilities proposed by Hartley and Rao (1962).
#' Estimator of class 3, it requires only first-order inclusion probabilities but
#' for all population units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for all population units
#' @inheritParams approx_var_est
#'
#'
#' @return a scalar, the estimated variance
#'

var_HartleyRao <- function(y, pik, sample) {

    ### Input values ---
    n <- length(y)
    piks <- pik[sample]

    ### Compute c_k values ---
    ck <- n / (n-1)
    ck <- ck * ( 1 - piks - sum(piks)/n + sum(pik**2)/n )

    ### Estimate variance ---
    B  <- sum(y/piks) / n
    ek <- y/piks - B

    v <- sum( ck*(ek**2) )

    ### Return result ---
    return( v )

}


#' Approximate Variance Estimators by Matei and Tillé (2005)
#'
#' Computes an approximate variance estimate according to
#' one of the estimators proposed by Matei and Tillé (2005)
#' Estimators of class 3, they require only first-order inclusion probabilities but
#' for all population units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for all population units
#' @param maxIter a scalar indicating the maximum number of iterations for the
#' fixed-point procedure
#' @param eps tolerance value for the convergence of the fixed-point procedure
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance

var_MateiTille <- function(y, pik, method, sample, maxIter=1000, eps=1e-05) {

    ### Check input ---
    method <- match.arg( method,
                         c( "MateiTille1",
                            "MateiTille2",
                            "MateiTille3",
                            "MateiTille4",
                            "MateiTille5" )
    )

    piks <- pik[sample]
    n <- length(y)
    N <- length(pik)


    if( identical( method, "MateiTille1") ){

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


        ### Estimate variance ---
        bks  <- bk[sample]
        ys   <- piks * sum( bks*y/(piks**2) ) / sum( bks/piks )

        v <- sum( bks * (y-ys)^2 / (piks**3) )
        v <- n*(N-1)  * v / ( N*(n-1) )

        ### Return result ---
        return( v )

    }else if( identical( method, "MateiTille2") ){

        Yht <- sum( y/piks )

        ### compute variance ---
        d  <- pik*(1-pik)
        dk <- d/sum(d)

        v <- sum( ( 1-piks ) * ( y/piks - Yht/n)^2 )
        v <- v / ( 1 - sum( (dk**2)/pik ) )

        ### Return result ---
        return( v )

    }else if( identical( method, "MateiTille3") ){

        ### compute variance ---
        d  <- (1-pik)
        dk <- pik*d / sum(pik*d)

        ds <- d[sample]

        v <- sum( ds*(y/piks) ) / sum(ds)
        v <- sum( ds * ( y/piks -  v )^2 )
        v <- v / ( 1 - sum( (dk**2)/pik ) )

        ### Return result ---
        return( v )

    }else if( identical( method, "MateiTille4") ){

        ### compute variance ---
        bk <- pik*(1-pik)*N / (N-1)
        bks <- bk[sample]

        ys <- piks * sum( bks*y / (piks**2) ) / sum( bks/piks )

        v <- sum( bks*( (y-ys)^2 ) / (piks**3) )
        v <- v / ( 1 - sum( bk/(n**2) ) )

        ### Return result ---
        return( v )

    }else if( identical( method, "MateiTille5") ){

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


        ### compute variance ---
        bks <- bk[sample]
        ys <- piks * sum( bks*y / (piks**2) ) / sum( bks/piks )

        v <- sum( bks*( (y-ys)^2 ) / (piks**3) )
        v <- v / ( 1 - sum( bk/(n**2) ) )

        ### Return result ---
        return( v )

    }


}

#' Approximate Variance Estimators by Brewer (2002)
#'
#' Compute an approximate variance estimate using the class of estimators
#' proposed by Brewer (2002).
#' Estimators of class 3, they require only first-order inclusion probabilities but
#' for all population units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for all population units
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance

var_Brewer_class3 <- function(y, pik, method, sample) {

    ### Check input ---
    method <- match.arg( method,
                         c(  "Brewer2",
                             "Brewer3",
                             "Brewer4" )
    )

    n <- length(y)

    ### Compute c_k values
    if( identical( method, "Brewer2") ){

        ck <- n - sum(pik**2)/n
        ck <- (n-1)/ck

    }else if( identical( method, "Brewer3") ){

        ck <- 1 - 2*pik/n + sum(pik**2)/(n**2)
        ck <- (n-1)/(n * ck)
        ck <- ck[sample]

    }else if( identical( method, "Brewer4") ){

        ck <- 1 - ((2*n-1)*pik + sum(pik**2)) / (n*(n-1))
        ck <- (n-1)/(n*ck)
        ck <- ck[sample]
    }

    ### Estimate variance ---
    piks <- pik[sample]

    v <- ( y/piks - sum(y/piks)/n )^2
    v <- sum( (1/ck - piks) * v )

    ### Return result ---
    return(v)
}



#' Tillé's Approximate Variance Estimator
#'
#' Compute an approximate variance estimate using the estimator
#' proposed by Tillé (1996).
#' Estimator of class 3, it requires only first-order inclusion probabilities but
#' for all population units.
#'
#' @param pik numeric vector of first-order inclusion probabilities for all population units
#' @param maxIter a scalar indicating the maximum number of iterations for the
#' fixed-point procedure
#' @param eps tolerance value for the convergence of the fixed-point procedure
#' @inheritParams approx_var_est
#'
#' @return a scalar, the estimated variance

var_Tille <- function(y, pik, sample, maxIter=1000, eps=1e-5) {

    ### Input values ---
    n <- length(y)

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

    ### Estimate variance ---
    bs <- b2[sample]
    piks <- pik[sample]

    yt <- y/piks
    wk <- piks/bs
    yw <- sum( wk*yt ) / sum(wk)
    v  <- sum(wk) * sum( wk*(yt-yw)^2 ) - n*sum( (yt - sum(y/piks)/n)^2)

    ### Return result ---
    return( v )

}



