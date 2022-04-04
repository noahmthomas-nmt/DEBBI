#' QLog
#' @description computes log density of Q(theta|lambda)
#' @param use_theta samples of theta from approximating distribution Q
#' @param use_lambda  parameters governing approximating distribution Q
#' @param S number of samples for monte carlo approximation
#' @param control_pars list of algo control parameters
#' @noRd


QLog=function(use_theta,use_lambda,control_pars,S){
  #returns log density Q(theta|lambda)
  out=0
  out=stats::dnorm(x=use_theta[1:(control_pars$n_pars_model)],
                   mean=use_lambda[1:(control_pars$n_pars_model)],
                   sd=exp(use_lambda[((control_pars$n_pars_model)+1):(control_pars$n_pars_dist)]),
                   log=T)
  out[(out==-Inf) | is.na(out)]=control_pars$neg_inf
  return(sum(out));
}
