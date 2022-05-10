#' CrossoverOptimize
#' @param chain_index  which chain you are updating
#' @param current_params current parameter values for chain (numeric vector)
#' @param current_weight  weight for current chain values
#' @param objFun  function you want to minimize (returns scalar)
#' @param step_size step size in DE jump
#' @param jitter_size noise
#' @param n_chains number of chains
#' @param crossover_rate rate of updating parameters on a give crossover step (0=no parameters updated, 1=all parameters updated)
#' @param ... additional arguments for objFun function
#' @noRd
#'


CrossoverOptimize <- function(chain_index,
                              current_params,
                              current_weight,
                              objFun,
                              step_size = .8,
                              jitter_size = 1e-6,
                              n_chains, crossover_rate = 1, ...) {

  # get statistics about chain
  weight_use <- current_weight[chain_index]
  params_use <- current_params[chain_index, ]
  len_param_use <- length(params_use)

  # use binomial to sample which params to update matching crossover rate frequency
  param_idices_bool <- stats::rbinom(len_param_use, prob = crossover_rate, size = 1)

  # if no params selected, randomly sample 1 parameter
  if (all(param_idices_bool == 0)) {
    param_idices_bool[sample(x = 1:len_param_use, size = 1)] <- 1
  }

  # indices of parameters to be updated
  param_indices <- seq(1, len_param_use, by = 1)[as.logical(param_idices_bool)]

  # sample parent chains
  parent_chain_indices <- sample(c(1:n_chains)[-chain_index], 3, replace = F)

  # mate parents for proposal
  params_use[param_indices] <- current_params[parent_chain_indices[3], param_indices] +
    step_size * (current_params[parent_chain_indices[1], param_indices] -
      current_params[parent_chain_indices[2], param_indices]) +
    stats::runif(length(param_indices), -jitter_size, jitter_size)

  params_use <- matrix(params_use, 1, len_param_use)

  # get log weight
  weight_proposal <- objFun(params_use, ...)
  if (is.na(weight_proposal)) weight_proposal <- -Inf

  # greedy acceptance rule
  if (weight_proposal > weight_use) {
    current_params[chain_index, ] <- params_use
    current_weight[chain_index] <- weight_proposal
  }

  return(c(current_weight[chain_index], current_params[chain_index, ]))
}
