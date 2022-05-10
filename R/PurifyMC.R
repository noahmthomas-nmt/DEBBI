#' PurifyMC
#'
#' @param chain_index which chain you are updating
#' @param current_theta current matrix of theta
#' @param current_log_post_like weight for current chain values
#' @param LogPostLike  log posterior density function (log likelihood + log prior)
#' @param n_chains number of chains
#' @param ... additional arguments for LogPostLike function
#' @noRd

PurifyMC <- function(chain_index, current_theta, current_log_post_like, LogPostLike, ...) {
  theta_use <- current_theta[chain_index, ]
  theta_use <- matrix(theta_use, 1, length(theta_use))

  like <- LogPostLike(theta_use, ...)

  if (!is.na(like) & like != -Inf) {
    current_log_post_like[chain_index] <- like
  }

  return(c(current_log_post_like[chain_index], current_theta[chain_index, ]))
}
