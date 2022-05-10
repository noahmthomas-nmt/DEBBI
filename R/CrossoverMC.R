#' CrossoverMC
#'
#' @param chain_index  which chain you are updating
#' @param param_indices  which parameters you are updating (int vector)
#' @param current_theta current parameter values for chain (numeric vector)
#' @param current_log_post_like  log posterior likelihood for current chain values
#' @param LogPostLike  log likelihood + log prior function (returns scalar)
#' @param step_size step size in DE jump
#' @param jitter_size noise
#' @param n_chains number of chains
#' @param ... additional arguments for LogPostLike function
#' @noRd



CrossoverMC <- function(chain_index, # which chain you are updating
                        param_indices, # which parameters you are updating (int vector)
                        current_theta, # current parameter values for chain (numeric vector)
                        current_log_post_like, # corresponding log post like for (numeric vector)
                        LogPostLike, # log likelihood function (returns scalar)
                        step_size = .8,
                        jitter_size = 1e-6,
                        n_chains, ...) {

  # get statistics about chain
  like_use <- current_log_post_like[chain_index]
  theta_use <- current_theta[chain_index, ]

  # sample parent chains
  parent_chain_indices <- sample(c(1:n_chains)[-chain_index], 2, replace = F)

  # mate parents for proposal
  theta_use[param_indices] <- theta_use[param_indices] +
    step_size * (current_theta[parent_chain_indices[1], param_indices] -
      current_theta[parent_chain_indices[2], param_indices]) +
    stats::runif(length(param_indices), -jitter_size, jitter_size)

  theta_use <- matrix(theta_use, 1, length(theta_use))

  # get log like
  like_proposal <- LogPostLike(theta_use, ...)
  if (is.na(like_proposal)) like_proposal <- -Inf

  # metropolis hasting acceptance rule
  if (stats::runif(1) < exp(like_proposal - like_use)) {
    current_theta[chain_index, ] <- theta_use
    current_log_post_like[chain_index] <- like_proposal
  }

  return(c(current_log_post_like[chain_index], current_theta[chain_index, ]))
}
