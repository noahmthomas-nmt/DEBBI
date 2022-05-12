#' KLHatGradient_LRVB
#' @description returns a Monte Carlo or quasi-Monte Carlo estimate of KL divergence gradient. lambda must contain variance parameterization of mean field
#' @param lambda samples of theta from approximating distribution Q
#' @param LogPostLike  log posterior likelihood function
#' @param S number of samples to use for the approximation
#' @param control_params list of algo control parameters
#' @param ... additional parameters for LogPostLike
#' @noRd
KLHatGradient_LRVB <- function(lambda, LogPostLike, control_params, S, ...) {
  # Monte Carlo approximation KL divergence up to a constant

  out <- numeric(control_params$n_params_dist) # initalize output vector
  lambda[(control_params$n_params_model + 1):control_params$n_params_dist] <-
    log(lambda[(control_params$n_params_model + 1):control_params$n_params_dist]^(1 / 2))
  # sample from q
  theta_mat <- QSample(use_lambda = lambda, control_params, S)

  # calc mean differences in log densities for theta_mat
  q_log_grad <- QLog_grad_wrt_lambda(theta_mat, use_lambda = lambda, control_params, S)
  q_log_density <- QLog(theta_mat, use_lambda = lambda, control_params, S)
  post_log_density <- apply(matrix(theta_mat,
                                   ncol = control_params$n_params_model,
                                   byrow = T),
  MARGIN = 1, FUN = LogPostLike, ...
  )

  for(i in 1:control_params$n_params_dist){
    out[i] <- sum(q_log_grad[i] * (q_log_density - post_log_density))/S
  }

  return(out)
}
