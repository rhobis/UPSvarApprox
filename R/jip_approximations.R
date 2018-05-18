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
#'
#'
#' @return A symmetric matrix of inclusion probabilities, which diagonal is the
#' vector of first-order inclusion probabilities.
#'
#'
#' @export



jip_approx <- function( pik, method ){

    ### Check input ---
    method <- match.arg(method, c( "Hajek", "Hartley-Rao", "IFPF", "Brewer") )

    if( !identical( class(pik), "numeric" ) ){
        stop( "pik should be a vector!")
    }else if( length(pik) < 2 ){
        stop( "The pik vector is too short!" )
    }else if( any(pik)<0  | any(pik)>1 ){
        stop( "Some pik values are outside the interval [0, 1]")
    }else if( !identical( sum(pik), as.integer( sum(pik) ) ) ){
        stop( "The sum of pik values is not an integer!")
    }

    ### Call method ---
    FUN <- switch(method,
           "Hajek"       = "jip_hajek",
           "Hartley-Rao" = "jip_hartleyrao",
           "IPFP"        = "jip_IPFP",
           "Brewer"      = "jip_brewer"
           )
    jips <- do.call( FUN, list( pik ) )


    ### Return result ---
    return( jips )
}

### ADD MORE METHODS ....



#brw_approx18
jip_brewer <- function(pik){
    #Approximation of joint incl. prob with formula 18 from
    #Brewer & Donadio (2003)
    n <- sum(pik)

    ci <- (n-1) / (n-(2*n-1)/(n-1) * pik + sum(pik**2)/(n-1))
    cc <- outer(ci,ci, '+')
    pp <- outer(pik,pik, '*')
    out <- cc*pp / 2
    diag(out) <- pik

    return(out)
}

#Function to compute Hajek approximation of joint probabilities
jip_hajek <- function(pik){
    #pik = vector of inclusion probabilities in the population

    d <- sum(pik*(1-pik))
    out <- outer(pik,pik,'*') * (1 - outer(1-pik, 1-pik, '*')/d)
    diag(out) <- pik
    return(out)
}

#Approximation of joint incl. prob by Hartley-Rao (1962)
#formula with precision of order O(N^{-4})
jip_hartleyrao <- function(pik){
    #
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
    return(out)
}


jip_IPFP <- function(pik,eps=1e-06, maxIter=1000){
    ### approximated second-order inclusion probability by
    ### means of the Iterative Proportional Fitting Procedure (IPFP)
    ### [see Tillé, 1996 - Some Remarks on Unequal Probability Sampling ...]
    if(any(pik>1)) stop("there are inclusion probabilities greater than one")

    ### Input values ---
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

    out <- outer(b2,b2,'*')
    diag(out) <- pik
    return(out)
}


