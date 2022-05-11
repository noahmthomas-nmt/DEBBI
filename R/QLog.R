#' QLog
#' @description computes log density of Q(theta|lambda)
#' @param use_theta samples of theta from approximating distribution Q
#' @param use_lambda  parameters governing approximating distribution Q
#' @param S number of samples for Monte Carlo approximation
#' @param control_params list of algo control parameters
#' @noRd


QLog <- function(use_theta, use_lambda, control_params, S = 1) {
  # returns log density Q(theta|lambda)
  out <- 0
  out <- stats::dnorm(
    x = use_theta,
    mean = rep(use_lambda[1:(control_params$n_params_model)], S),
    sd = rep(exp(use_lambda[((control_params$n_params_model) + 1):(control_params$n_params_dist)]), S),
    log = TRUE
  )
  out[(out == -Inf) | is.na(out)] <- control_params$neg_inf
  return(out)
}
