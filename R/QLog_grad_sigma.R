#' QLog_grad_sigma
#' @description computes gradient log density of the mean field Q(theta|lambda), uses sigma (standard deviation) parameterization
#' @param use_theta samples of theta from approximating distribution Q
#' @param use_lambda  parameters governing approximating distribution Q
#' @param S number of samples for monte carlo approximation
#' @param control_pars list of algo control parameters
#' @noRd


QLog_grad_sigma <- function(use_theta, use_lambda, control_pars, S = 1) {
  x <- use_theta
  mu <- use_lambda[1:(control_pars$n_params_model)]
  sigma <- exp(use_lambda[((control_pars$n_params_model) + 1):(control_pars$n_params_dist)])
  mu_partials <- (x - mu) / (sigma^2)
  sigma_partials <- ((x - mu)^2) / (sigma^3) - 1 / sigma

  return(c(mu_partials, sigma_partials))
}
