#' QSample
#' @description samples or quasi-randomly samples theta from Q(theta|lambda)
#' @param use_lambda  parameters governing approximating distribution Q
#' @param control_params list of algo control parameters
#' @param S number of samples (integer)
#' @noRd

QSample <- function(use_lambda, control_params, S) {
  # returns a collapsed vector for a S by n_params_model matrix sampled from Q(theta|lambda)
  if (control_params$use_QMC == TRUE) {
    if (control_params$quasi_rand_seq == "sobol") quantileMat <- c(t(randtoolbox::sobol(S, control_params$n_params_model)))
    if (control_params$quasi_rand_seq == "halton") quantileMat <- c(t(randtoolbox::halton(S, control_params$n_params_model)))
  } else {
    quantileMat <- stats::runif(S * control_params$n_params_model, 0, 1)
  }

  out <- stats::qnorm(quantileMat, rep(use_lambda[1:control_params$n_params_model], S), rep(exp(use_lambda[(control_params$n_params_model + 1):(control_params$n_params_dist)]), S))
  return(out)
}
