#' PurifyVI
#'
#' @param chain_index  which chain you are updating
#' @param current_params current parameter values for chain (numeric vector)
#' @param current_weight  weight for current chain values
#' @param LogPostLike  log posterior density function (log likelihood + log prior)
#' @param n_chains number of chains
#' @param control_params list of algorithm control parameters
#' @param S number of samples for ELBO Monte Carlo estimate
#' @param ... additional arguments for LogPostLike function
#' @noRd
PurifyVI <- function(chain_index,
                     current_params,
                     current_weight,
                     LogPostLike,
                     n_chains,
                     control_params,
                     S, ...) {

  # get statistics about chain
  weight_use <- current_weight[chain_index]
  params_use <- current_params[chain_index, ]
  len_param_use <- length(params_use)

  params_use <- matrix(params_use, 1, len_param_use)

  # get log weight
  weight_proposal <- ELBO(
    params_use,
    LogPostLike,
    control_params,
    S, ...
  )

  if (is.na(weight_proposal)) weight_proposal <- -Inf

  # greedy acceptance rule
  if (is.finite(weight_proposal)) {
    current_params[chain_index, ] <- params_use
    current_weight[chain_index] <- weight_proposal
  }

  return(c(current_weight[chain_index], current_params[chain_index, ]))
}
