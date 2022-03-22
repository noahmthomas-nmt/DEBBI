
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DEBBI

<!-- badges: start -->
<!-- badges: end -->

The primary goal of DEBBI (short for ‘Differential Evolution-based
Bayesian Inference’) is to provide efficient access to and enable
reproducible use of Bayesian inference algorithms such as Differential
Evolution Markov Chain Monte Carlo, Differential Evolution Variation
Inference, and Differential Evolution maximum a posteriori estimation.
The second goal of this package is to be compatible with likelihood-free
Bayesian inference methodologies that use approximations of the
likelihood function such as probability density approximation,
kernel-based approximate Bayesian computation, and synthetic likelihood.

## Installation

You can install the development version of DEBBI from
[GitHub](https://github.com/bmgaldo/DEBBI) with:

``` r
# install.packages("devtools")
devtools::install_github("bmgaldo/DEBBI")
```

## DEMCMC Example

Estimate mean parameters of two independent normal distributions with
known standard deviations using DEMCMC

``` r
library(DEBBI)

# simulate from model
dataExample=matrix(rnorm(100,c(-1,1),c(1,1)),nrow=50,ncol=2,byrow = T)

# list parameter names
par_names_example=c("mu_1","mu_2")

# log posterior likelihood function = log likelihood + log prior | returns a scalar
LogPostLikeExample=function(x,data,par_names){
  out=0
  
  names(x)<-par_names
  
  # log prior
  out=out+sum(dnorm(x["mu_1"],0,sd=1,log=T))
  out=out+sum(dnorm(x["mu_2"],0,sd=1,log=T))
  
  # log likelihoods
  out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=T))
  out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=T))
  
  return(out)
}

# Sample from posterior
post <- DEMCMC(LogPostLike=LogPostLikeExample,
                  control_pars=AlgoParsDEMCMC(n_pars=length(par_names_example),
                                     n_iter=500, 
                                     n_chains=12,
                                     init_sd=.01,
                                     init_center=0,
                                     n_cores_use=1,
                                     step_size=NULL,
                                     jitter_size=1e-6,
                                     burnin=100,
                                     parallelType = "none"),
                  data=dataExample,
                  par_names = par_names_example)
#> [1] "initalizing chains..."
#> [1] "1 / 12"
#> [1] "2 / 12"
#> [1] "3 / 12"
#> [1] "4 / 12"
#> [1] "5 / 12"
#> [1] "6 / 12"
#> [1] "7 / 12"
#> [1] "8 / 12"
#> [1] "9 / 12"
#> [1] "10 / 12"
#> [1] "11 / 12"
#> [1] "12 / 12"
#> [1] "chain initialization complete  :)"
#> [1] "running DEMCMC"
#> [1] "iter 100/500"
#> [1] "iter 200/500"
#> [1] "iter 300/500"
#> [1] "iter 400/500"
#> [1] "iter 500/500"


par(mfrow=c(2,2))

  hist(post$samples[,,1],main="marginal posterior distribution",
       xlab=par_names_example[1],prob=T)
  matplot(post$samples[,,1],type='l',ylab=par_names_example[1],
          main="chain trace plot",xlab="iteration")

  hist(post$samples[,,2],xlab=par_names_example[2],prob=T,main="")
  matplot(post$samples[,,2],type='l',ylab=par_names_example[2],xlab="iteration")
```

<img src="man/figures/README-DEMCMC example-1.png" width="100%" />
