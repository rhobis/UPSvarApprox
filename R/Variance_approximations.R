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
#' The variance approximations available in this function are described below,
#' the notation used is that of Matei and Tillé (2005).
#'
#' \itemize{
#'     \item Hájek variance approximation (\code{method="Hajek1"}):
#'     \deqn{ \tilde{Var} = \sum_{i \in U} \frac{b_i}{\pi_i^2}(y_i - y_i^*)^2  }{
#'     \sum b (y - y*)^2 / (\pi^2)}
#'     where
#'     \deqn{y_i^* = \pi_i \frac{ \sum_{j\in U} b_j y_j/\pi_j }{ \sum_{j \in U} b_j } }{
#'     y* = \pi (\sum b*y/\pi) / (\sum b)}
#'     and
#'     \deqn{ b_i = \frac{ \pi_i(1-\pi_i)N }{ N-1 } }{ b= (N\pi(\-\pi)) / (N-1)}
#'
#'     \item Starting from Hajék (1964), Brewer (2002) defined the following estimator
#'     (\code{method="Hajek2"}):
#'     \deqn{ \tilde{Var} = \sum_{i \in U} \pi_i(1-\pi_i) \Bigl( \frac{y_i}{\pi_i} -
#'     \frac{\tilde{Y}}{n} \Bigr)^2 }{
#'     \sum \pi(1-\pi) ( y/\pi - Y*/n)^2 }
#'     where \eqn{\tilde{Y} = \sum{i \in U} a_i y_i}{ Y* = \sum a*y }
#'     and \eqn{a_i = n(1-\pi_i)/\sum_{j \in U} \pi_j(1-\pi_j) }{ a = n(1-\pi) / \sum(\pi(1-\pi))}
#'
#'     \item Hartley and Rao (1962) variance approximation (\code{method="HartleyRao1"}):
#'     \deqn{\tilde{Var} = \sum_{i \in U} \pi_i \Bigl( 1 - \frac{n-1}{n}\pi_i \Bigr) \Bigr( \frac{y_i}{\pi_i} - \frac{Y}{n}  \Bigr)^2
#'     - \frac{n-1}{n^2} \sum_{i \in U} \Biggl( 2\pi_i^3 - \frac{\pi_i^2}{2}\sum_{j \in U} \pi_j^2 \Biggr)\Bigr( \frac{y_i}{\pi_i} - \frac{Y}{n}  \Bigr)^2
#'     + \frac{2(n-1)}{n^3} \Biggl( \sum_{i \in U}\pi_i y_i - \frac{Y}{n}\sum_{i\in U} \pi_i^2 \Biggr)^2 }{
#'     *see pdf version of documentation*}
#'
#'     \item Hartley and Rao (1962) provide a simplified version of the
#'     variance above (\code{method="HartleyRao2"}):
#'     \deqn{ \tilde{Var} = \sum_{i \in U} \pi_i \Bigl( 1 - \frac{n-1}{n}\pi_i \Bigr) \Bigr( \frac{y_i}{\pi_i} - \frac{Y}{n}  \Bigr)^2 }{
#'      Var = \sum \pi ( 1 - ( (n-1)/n )\pi )( y/\pi - Y/n )^2 }
#'
#'      \item \code{method="FixedPoint"} compute the Fixed-Point variance approximation
#'      proposed by Deville and Tillé (2005).
#'      The variance can be expressed in the same form as in \code{method="Hajek1"},
#'      and the coefficients \eqn{b_i}{b} are computed iteratively by the algorithm:
#'      \enumerate{
#'          \item \deqn{b_i^{(0)} = \pi_i (1-\pi_i) \frac{N}{N-1}, \,\, \forall i \in U }{ b0 = \pi(1-\pi)(N/(N-1))}
#'          \item \deqn{ b_i^{(k)} = \frac{(b_i^{(k-1)})^2 }{\sum_{j\in U} b_j^{(k-1)} } + \pi_i(1-\pi_i) }{
#'          b(k) = [ b(i-1) ]^2 / [ \sum b(i-1) ] + \pi(1-\pi) }
#'       }
#'       a necessary condition for convergence is checked and, if not satisfied,
#'       the function returns an alternative solution that uses only one iteration:
#'       \deqn{b_i = \pi_i(1-\pi_i)\Bigl( \frac{N\pi_i(1-\pi_i)}{ (N-1)\sum_{j\in U}\pi_j(1-\pi_j) } + 1 \Bigr) }{
#'       b = \pi(1-\pi)( 1 + (N\pi(1-\pi)) / ( (N-1) \sum \pi(1-\pi) ) ) }
#'
#' }
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
#' Var_approx(y, pik, n, method = "Hajek2")
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
