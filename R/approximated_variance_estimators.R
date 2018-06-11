#' Approximated Variance Estimators
#'
#' Approximated variance estimators which use only first-order inclusion probabilities
#'
#' @param y numeric vector of sample observations
#' @param pik numeric vector of first-order inclusion probabilities, of length *N* or *n*,
#' depending on the chosen method (see Details for more information)
#' @param method string indicating the desired approximate variance estimator.
#' One of "Deville1", "Deville2", "Deville3", "Hajek", "Rosen", "FixedPoint",
#' "Brewer1", "HartleyRao", "Berger", "Tille", "MateiTille1", "MateiTille2",
#' "MateiTille3", "MateiTille4", "MateiTille5", "Brewer2", "Brewer3", "Brewer4".
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
#' The choice of estimator to be used is made through the argument \code{method},
#' the list of methods and their respective equations is presented below.
#'
#' Matei and Tillé (2005) divide the approximated variance estimators into
#' three classes, depending on the quantities they require:
#'
#' \enumerate{
#' \item First and second-order inclusion probabilities:
#' The first class is composed of the Horvitz-Thompson estimator (Horvitz and Thompson 1952)
#' and the Sen-Yates-Grundy estimator (Yates and Grundy 1953; Sen 1953),
#' which are available through function \code{\link[sampling]{varHT}};
#'
#' \item Only first-order inclusion probabilities and only for sample units;
#'
#' \item Only first-order inclusion probabilities, for the entire population.
#' }
#'
#'
#' Haziza, Mecatti and Rao (2008) provide a common form to express most of the
#' estimators in class 2 and 3. The common form proposed is
#'
#' \deqn{\hat{var}(\hat{t}_{HT}) = \sum_{i \in s}c_i e_i^2 }{ var = \sum_s c*e^2}
#'
#' where \eqn{ e_i = \frac{y_i}{\pi_i} - \hat{B} }{ e =  y/\pi - B}, with
#'
#' \deqn{ \hat{B} = \frac{\sum_{i\in s} a_i (y_i/\pi_i) }{\sum_{i\in s} a_i} }{
#' \sum a*(y/\pi) / \sum (a)}
#'
#' and \eqn{a_i}{a} and \eqn{c_i}{c} are parameters that define the different
#' estimators:
#'
#' \itemize{
#' \item \code{method="Hajek"} [Class 2]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i)}{c = (n/(n-1)) (1-\pi)}
#'     \deqn{a_i= c_i}{a=c}
#' \item \code{method="Deville2"} [Class 2]
#'     \deqn{c_i = (1-\pi_i)\Biggl\{ 1 - \sum_{j\in s}\Bigl[ \frac{1-\pi_j}{\sum_{k\in s} (1-\pi_k)} \Bigr]^2 \Biggr\}^{-1} }{
#'     c = (1-\pi) [1 - \sum ( (1-pi)/\sum (1-pi))^2)]^(-1) }
#'     \deqn{a_i= c_i}{a=c}
#' \item \code{method="Deville3"} [Class 2]
#'     \deqn{c_i = (1-\pi_i)\Biggl\{ 1 - \sum_{j\in s}\Bigl[ \frac{1-\pi_j}{\sum_{k\in s} (1-\pi_k)} \Bigr]^2 \Biggr\}^{-1} }{
#'     c = (1-\pi)[1 - \sum ( (1-pi)/\sum (1-pi))^2)]^(-1) }
#'     \deqn{a_i = 1}{a=1}
#' \item \code{method="Rosen"} [Class 2]
#'     \deqn{c_i = \frac{n}{n-1} (1-\pi_i)}{ c = (n/(n-1)) (1-\pi) }
#'     \deqn{a_i= (1-\pi_i)log(1-\pi_i) / \pi_i]}{a= (1-\pi)log(1-\pi) / \pi}
#' \item \code{method="Brewer1"} [Class 2]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i)}{ c = (n/(n-1)) (1-\pi)}
#'     \deqn{a_i= 1}{a=1}
#' \item \code{method="Brewer2"} [Class 3]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i+ \frac{\pi_i}{n} - n^{-2}\sum_{j \in U} \pi_j^2)}{
#'     c = (n/(n-1)) (1-\pi+\pi/n - \sum_U \pi/(n^2))}
#'     \deqn{a_i= 1}{a=1}
#' \item \code{method="Brewer3"} [Class 3]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i - \frac{\pi_i}{n} - n^{-2}\sum_{j \in U} \pi_j^2)}{
#'     c = (n/(n-1)) (1-\pi - \pi/n - \sum_U \pi/(n^2))}
#'     \deqn{a_i= 1}{a=1}
#' \item \code{method="Brewer4"} [Class 3]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i - \frac{\pi_i}{n-1} + n^{-1}(n-1)^{-1}\sum_{j \in U} \pi_j^2)}{
#'     c= (n/(n-1)) (1-\pi - \pi/(n-1) + \sum_U \pi^2 / (n*(n-1)))}
#'     \deqn{a_i= 1}{a=1}
#' \item \code{method="Berger"} [Class 3]
#'     \deqn{c_i = \frac{n}{n-1} (1-\pi_i) \Biggl[ \frac{\sum_{j\in s} (1-\pi_j)}{\sum_{j\in U}} (1-\pi_j) \Biggr] }{
#'     c = (n/(n-1)) (1-\pi) (\sum_s (1-\pi) / \sum_U (1-pi)) }
#'     \deqn{a_i = 1}{a=1}
#' \item \code{method="HartleyRao"} [Class 3]
#'     \deqn{c_i = \frac{n}{n-1}(1-\pi_i - n^{-1}\sum_{j \in s}\pi_i + n^{-1}\sum_{j\in U} \pi_j^2)}{
#'     c= (n/(n-1)) (1-\pi - \sum_s \pi/n + \sum_U \pi^2 / n) }
#'     \deqn{a_i= 1}{a=1}
#'}
#'
#'
#' Some additional estimators are defined in Matei and Tillé (2005):
#'
#' \itemize{
#'     \item \code{method="Deville1"} [Class 2]
#'     \deqn{\hat{var}(\hat{t}_{HT}) = \sum_{i \in s} \frac{c_i}{ var = \pi_i^2} (y_i - y_i^*)^2 }{
#'     \sum (c/\pi^2) (y-y*)^2 }
#'
#'     where
#'     \deqn{ y_i^* = \pi_i \frac{\sum_{j \in s} c_j y_j / \pi_j}{\sum{j \in s} c_j} }{
#'     y* = \pi (\sum c*y/\pi) / (\sum c)}
#'     and \eqn{c_i = (1-\pi_i)\frac{n}{n-1} }{c = (n/(n-1)) (1-\pi)}
#'
#'
#'  \item \code{method="Tille"} [Class 3]
#'      \deqn{
#'      \hat{var}(\hat{t}_{HT}) = \biggl( \sum_{i \in s} \omega_i \biggr)
#'      \sum_{i\in s} \omega_i (\tilde{y}_i - \bar{\tilde{y}}_\omega )^2
#'      - n \sum_{i\in s}\biggl( \tilde{y}_i - \frac{\hat{t}_{HT}}{n} \biggr)^2
#'      }{
#'      var = ( \sum \omega ) \sum \omega( y* - \gamma )^2 - n\sum (y* - \sum(y/\pi)/n )^2}
#'
#'      where  \eqn{\tilde{y}_i = y_i / \pi_i }{ y* = y/\pi},
#'      \eqn{\omega_i = \pi_i / \beta_i}{ \omega = \pi/\beta}
#'      and \eqn{\bar{\tilde{y}}_\omega =
#'                   \biggl( \sum_{i \in s} \omega_i \biggr)^{-1} \sum_{i \in s} \omega_i \tilde{y}_i }{
#'               \gamma = \sum(\omega y*)/\sum(\omega) }
#'
#'    The coefficients \eqn{\beta_i}{\beta} are computed iteratively through the
#'    following procedure:
#'     \enumerate{
#'         \item \eqn{\beta_i^{(0)} = \pi_i, \,\, \forall i\in U}{\beta(0) = \pi, i = 1, ..., N}
#'         \item \eqn{ \beta_i^{(2k-1)} = \frac{(n-1)\pi_i}{\beta^{(2k-2)} - \beta_i^{(2k-2)}}  }{
#'                     \beta(2k-1) = ( (n-1)\pi )/(\sum\beta(2k-2) - \beta(2k-2)) }
#'         \item \eqn{\beta_i^{2k} = \beta_i^{(2k-1)}
#'         \Biggl( \frac{n(n-1)}{(\beta^(2k-1))^2 - \sum_{i\in U} (\beta_k^{(2k-1)})^2 } \Biggr)^(1/2) }{
#'         \beta(2k) = \beta(2k-1) ( n(n-1) / ( (\sum\beta(2k-1))^2 - \sum( \beta(2k-1)^2 ) ) )^(1/2) }
#'     }
#'     \eqn{  \text{with} \beta^{(k)} = \sum_{i\in U} \beta_i^{i}, \,\, k=1,2,3, \dots }{}
#'
#'  \item \code{method="MateiTille1"} [Class 3]
#'      \deqn{\hat{var}(\hat{t}_{HT}) = \frac{n(N-1)}{N(n-1)} \sum_{i\in s} \frac{b_i}{\pi_i^3} (y_i - \hat{y}_i^*)^2 }{
#'            var = n(N-1)/(N(n-1)) \sum (b/\pi^3) (y-y*)^2 }
#'      where \deqn{\hat{y}_i^* = \pi_i \frac{\sum_{i\in s} b_i y_i/\pi_i^2}{\sum_{i\in s} b_i/\pi_i} }{
#'      y* = \pi (\sum b*y/(\pi^2)) / (\sum b/\pi )}
#'
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
#'  \item \code{method="MateiTille2"} [Class 3]
#'      \deqn{ \hat{var}(\hat{t}_{HT}) = \frac{1}{1 - \sum_{i\in U} \frac{d_i^2}{\pi_i} }
#'      \sum_{i\in s} (1-\pi_i) \Bigl( \frac{y_i}{\pi_i} - \frac{\hat{t}_{HT}}{n} \Bigr)^2 }{
#'      var = ( \sum (1-\pi)( y/\pi - \sum(y/\pi)/n)^2  ) / (1-\sum (d^2/\pi)) }
#'
#'      where
#'      \deqn{ d_i = \frac{\pi_i(1-\pi_i)}{\sum_{j\in U} \pi_j(1-\pi_j) } }{
#'      d = (\pi(1-\pi)) / (\sum \pi(1-\pi)) }
#'
#'  \item \code{method="MateiTille3"} [Class 3]
#'      \deqn{ \hat{var}(\hat{t}_{HT}) = \frac{1}{1 - \sum_{i\in U} \frac{d_i^2}{\pi_i} }
#'      \sum_{i\in s} (1-\pi_i) \Bigl( \frac{y_i}{\pi_i} -
#'      \frac{ \sum_{j\in s} (1-\pi_j)\frac{y_j}{\pi_j} }{ \sum_{j\in s} (1-\pi_j)  } \Bigr)^2 }{
#'      var = ( \sum (1-\pi)( y/\pi - (\sum(1-\pi)(y/\pi))/(\sum(1-\pi))  )^2  ) / (1-\sum (d^2/\pi)) }
#'
#'      where \eqn{d_i}{d} is defined as in \code{method="MateiTille2"}.
#'
#'  \item \code{method="MateiTille4"} [Class 3]
#'      \deqn{ \hat{var}(\hat{t}_{HT}) = \frac{1}{1 - \sum_{i\in U} b_i/n^2 }
#'      \sum_{i\in s} \frac{b_i}{\pi_i^3} (y_i - \hat{y}_i^* )^2  }{
#'      var =  \sum ( (b/(\pi^3)) (y-y*)^2 )
#'      }
#'      where
#'      \deqn{  \hat{y}_i^* = \pi_i \frac{ \sum_{j\in s} b_j y_j/\pi_j^2 }{  \sum_{j\in s} b_j/\pi_j } }{
#'       y* = \pi (\sum b*y/(\pi^2)) / (\sum b/\pi )}
#'      and
#'      \deqn{ b_i = \frac{ \pi_i(1-\pi_i)N }{ N-1 } }{ b = ( \pi(1-\pi)N ) / (N-1)}
#'
#'  \item \code{method="MateiTille5"} [Class 3]
#'  This estimator is defined as in \code{method="MateiTille4"}, and the \eqn{b_i}{b}
#'  values are defined as in \code{method="MateiTille1"}
#' }
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
#' approx_var_est(ys, piks, method="Hajek")
#' approx_var_est(ys, piks, method="Rosen")
#' approx_var_est(ys, piks, method="FixedPoint")
#' approx_var_est(ys, piks, method="Brewer1")
#'
#' ### Estimators of class 3 ---
#' approx_var_est(ys, pik, method="HartleyRao", sample=s)
#' approx_var_est(ys, pik, method="Berger", sample=s)
#' approx_var_est(ys, pik, method="Tille", sample=s)
#' approx_var_est(ys, pik, method="MateiTille1", sample=s)
#' approx_var_est(ys, pik, method="MateiTille2", sample=s)
#' approx_var_est(ys, pik, method="MateiTille3", sample=s)
#' approx_var_est(ys, pik, method="MateiTille4", sample=s)
#' approx_var_est(ys, pik, method="MateiTille5", sample=s)
#' approx_var_est(ys, pik, method="Brewer2", sample=s)
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
                             "Hajek",
                             "Rosen",
                             "FixedPoint",
                             "Brewer1",
                             #third class
                             "HartleyRao",
                             "Berger",
                             "Tille",
                             "MateiTille1",
                             "MateiTille2",
                             "MateiTille3",
                             "MateiTille4",
                             "MateiTille5",
                             "Brewer2",
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

    if( !(class(y) %in% c("numeric", "integer")) ){
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
    if( method %in% c( "Deville1", "Deville2", "Deville3", "Rosen", "Hajek",
                       "FixedPoint", "Brewer1") ){
        if( !identical( length(y), length(pik) ) )
            stop( "Method ", "'", method, "' ",
                  "requires argument 'y' and 'pik' to have the same length, ",
                  "which should be equal to the sample size!"
            )

    }

    #third class of estimators
    if( method %in% c( "HartleyRao", "Berger", "Tille", "MateiTille1",
                       "MateiTille2", "MateiTille3", "MateiTille4", "MateiTille5",
                       "Brewer2", "Brewer3", "Brewer4" ) ){

        if( !is.wholenumber( sum(pik) ) )
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
                "Hajek"       = var_Hajek(y, pik),
                "Rosen"       = var_Rosen(y, pik),
                "FixedPoint"  = var_FixedPoint(y, pik, eps=1e-5, maxIter=1000),
                "Brewer1"     = var_Brewer_class2(y, pik),
                #third class
                "HartleyRao"  = var_HartleyRao(y, pik, sample),
                "Berger"      = var_Berger(y, pik, sample),
                "Tille"       = var_Tille(y, pik, sample, eps=1e-5, maxIter=1000),
                "MateiTille1" = var_MateiTille(y, pik, method, sample) ,
                "MateiTille2" = var_MateiTille(y, pik, method, sample),
                "MateiTille3" = var_MateiTille(y, pik, method, sample),
                "MateiTille4" = var_MateiTille(y, pik, method, sample),
                "MateiTille5" = var_MateiTille(y, pik, method, sample),
                "Brewer2"     = var_Brewer_class3(y, pik, method, sample),
                "Brewer3"     = var_Brewer_class3(y, pik, method, sample),
                "Brewer4"     = var_Brewer_class3(y, pik, method, sample)
    )

    # do.call( FUN, list(y, pik))

    ### Return result ---
    return( v )
}



