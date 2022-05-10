#' ELBO
#' @description returns a monte carlo or quasi-monte carlo estimate of the ELBO
#' @param lambda samples of theta from approximating distribution Q
#' @param LogPostLike  log posterior likelihood function
#' @param control_pars list of algo control parameters
#' @param S number of samples to use for the approximation
#' @param ... additional parameters for LogPostLike
#' @noRd
ELBO <- function(lambda, LogPostLike, control_pars, S, ...) {
  # monte carlo approximation ELBO
  out <- KLHat(lambda, LogPostLike, control_pars, S, ...) * -1
  return(out)
}
