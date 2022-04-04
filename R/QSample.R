#' QSample
#' @description samples or quasi-randomly samples theta from Q(theta|lambda)
#' @param use_lambda  parameters governing approximating distribution Q
#' @param control_pars list of algo control parameters
#' @param S number of samples (integer)
#' @noRd

QSample=function(use_lambda,control_pars,S){
  #returns a S by n_pars_model matrix sampled from Q(theta|lambda)
  out=matrix(NA,S,control_pars$n_pars_model)
  if(control_pars$use_QMC==T){
    if(control_pars$quasiSeq=='sobol')quantileMat=randtoolbox::sobol(S,control_pars$n_pars_model)
    if(control_pars$quasiSeq=='halton')quantileMat=randtoolbox::halton(S,control_pars$n_pars_model)
  } else {
    quantileMat=matrix(stats::runif(S,0,1),S,control_pars$n_pars_model)
  }
  for(i in 1:n_pars_model){
    qs=quantileMat[,i]
    out[,i]=stats::qnorm(qs,use_lambda[i],exp(use_lambda[control_pars$n_pars_model+i]))
  }
  return(out)
}
