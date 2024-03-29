UPSvarApprox
======================================================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/UPSvarApprox)](https://cran.r-project.org/package=UPSvarApprox)
[![](https://cranlogs.r-pkg.org/badges/grand-total/UPSvarApprox)](https://cran.r-project.org/package=UPSvarApprox)
[![R badge](https://img.shields.io/badge/-Support%20me-brightgreen)](https://www.buymeacoffee.com/rhobis)

Description 
-----------------

UPSvarApprox provides functions for the approximation of the variance of the 
Horvitz-Thompson total estimator in Unequal Probability Sampling
using only first-order inclusion probabilities.

The main functions are:

- `Var_approx()`: computes and approximation of the variance of the HT estimator; 
- `approx_var_est()`: computes an approximate variance estimate for the HT estimator;



Installation
------------

The development version of the package can be installed from GitHub:

``` r
# if not present, install 'devtools' package
install.packages("devtools")
devtools::install_github("rhobis/UPSvarApprox")
```

Usage
-----

``` r
library(UPSvarApprox)

### Generate population data ---
N <- 500; n <- 50

set.seed(0)
x <- rgamma(500, scale=10, shape=5)
y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )

pik <- n * x/sum(x)
s   <- sample(N, n)

ys <- y[s]
piks <- pik[s]

### Variance approximations ---
Var_approx(y, pik, n, method = "Hajek1")
Var_approx(y, pik, n, method = "Hajek1")
Var_approx(y, pik, n, method = "HartleyRao1")
Var_approx(y, pik, n, method = "HartleyRao2")
Var_approx(y, pik, n, method = "FixedPoint")


### Approximate variance estimators ---

## Estimators of class 2
approx_var_est(ys, piks, method="Deville1")
approx_var_est(ys, piks, method="Deville2")
approx_var_est(ys, piks, method="Deville3")
approx_var_est(ys, piks, method="Rosen")
approx_var_est(ys, piks, method="FixedPoint")
approx_var_est(ys, piks, method="Brewer1")

## Estimators of class 3 
approx_var_est(ys, pik, method="Berger", sample=s)
approx_var_est(ys, pik, method="Tille", sample=s)
approx_var_est(ys, pik, method="MateiTille1", sample=s)
approx_var_est(ys, pik, method="MateiTille2", sample=s)
approx_var_est(ys, pik, method="MateiTille3", sample=s)
approx_var_est(ys, pik, method="MateiTille4", sample=s)
approx_var_est(ys, pik, method="MateiTille5", sample=s)
approx_var_est(ys, pik, method="Brewer2", sample=s)
approx_var_est(ys, pik, method="Brewer3", sample=s)
approx_var_est(ys, pik, method="Brewer4", sample=s)

```

More
----

- Please, report any bug or issue [here](https://github.com/rhobis/UPSvarApprox/issues).
- For more information, please contact the maintainer at `rob.sichera@gmail.com`. 

<br/>

<a href="https://www.buymeacoffee.com/rhobis" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png" alt="Buy Me A Coffee" width="217" height="60"></a>

<br/>

