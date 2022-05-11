#' KLHatforLRVB
#' @description returns a Monte Carlo or quasi-Monte Carlo estimate of KL divergence up to a constant (negative ELBO). lambda must contain log variance
#' @param lambda samples of theta from approximating distribution Q
#' @param LogPostLike  log posterior likelihood function
#' @param S number of samples to use for the approximation
#' @param control_params list of algo control parameters
#' @param ... additional parameters for LogPostLike
#' @noRd
KLHatforLRVB <- function(lambda, LogPostLike, control_params, S, ...) {
  # Monte Carlo approximation KL divergence up to a constant

  out <- 0 # initalize output vector
  lambda[(control_params$n_params_model + 1):control_params$n_params_dist] <-
    log(lambda[(control_params$n_params_model + 1):control_params$n_params_dist]^(1 / 2))
  # sample from q
  theta_mat <- QSample(use_lambda = lambda, control_params, S)

  # calc mean differences in log densities for theta_mat
  q_log_density <- sum(QLog(theta_mat, use_lambda = lambda, control_params, S)) / S
  post_log_density <- mean(apply(matrix(theta_mat,
    ncol = control_params$n_params_model, byrow = T
  ),
  MARGIN = 1, FUN = LogPostLike, ...
  ))
  out <- q_log_density - post_log_density


  return(out)
}
