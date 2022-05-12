#' QLog_grad_wrt_lambda
#' @description computes gradient log density of the mean field Q(theta|lambda), uses log sigma (unconstrained) parameterization
#' @param use_theta samples of theta from approximating distribution Q
#' @param use_lambda  parameters governing approximating distribution Q
#' @param S number of samples for Monte Carlo approximation
#' @param control_params list of algo control parameters
#' @noRd


QLog_grad_log_sigma <- function(use_theta, use_lambda, control_params, S = 1) {
  # returns lo g density Q(theta|lambda)

  x <- use_theta
  mu <- use_lambda[1:(control_params$n_params_model)]
  sigma <- exp(use_lambda[((control_params$n_params_model) + 1):(control_params$n_params_dist)])
  mu_partials <- (x - mu) / (sigma^2)
  sigma_partials <- ((x - mu)^2) / (sigma^3) - 1 / sigma
  log_sigma_partials <- sigma_partials * (1 / sigma)
  return(c(mu_partials, log_sigma_partials))
}
