#' QSample
#' @description samples or quasi-randomly samples theta from Q(theta|lambda)
#' @param use_lambda  parameters governing approximating distribution Q
#' @param control_pars list of algo control parameters
#' @param S number of samples (integer)
#' @noRd

QSample=function(use_lambda,control_pars,S){
  # returns a collapsed vector for a S by n_pars_model matrix sampled from Q(theta|lambda)
  if(control_pars$use_QMC==T){
    if(control_pars$quasi_rand_seq=='sobol')quantileMat=c(t(randtoolbox::sobol(S,control_pars$n_pars_model)))
    if(control_pars$quasi_rand_seq=='halton')quantileMat=c(t(randtoolbox::halton(S,control_pars$n_pars_model)))
  } else {
    quantileMat=stats::runif(S*control_pars$n_pars_model,0,1)
  }

  out=stats::qnorm(quantileMat,rep(use_lambda[1:control_pars$n_pars_model],S),rep(exp(use_lambda[(control_pars$n_pars_model+1):(control_pars$n_pars_dist)]),S))
  return(out)
}
