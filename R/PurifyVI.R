#' PurifyVI
#'
#' @param chain_index  which chain you are updating
#' @param current_pars current parameter values for chain (numeric vector)
#' @param current_weight  weight for current chain values
#' @param LogPostLike  log posterior density function (log likelihood + log prior)
#' @param n_chains number of chains
#' @param control_pars list of algorithm control parameters
#' @param S number of samples for ELBO monte carlo estimate
#' @param ... additional arguments for LogPostLike function
#' @noRd
PurifyVI=function(chain_index,
                  current_pars,
                  current_weight,
                  LogPostLike,
                  n_chains,
                  control_pars,
                  S, ... ){

  # get statistics about chain
  weight_use = current_weight[chain_index]
  pars_use = current_pars[chain_index,]
  len_par_use = length(pars_use)

  pars_use = matrix(pars_use,1,len_par_use)

  # get log weight
  weight_proposal = ELBO(pars_use,
                         LogPostLike,
                         control_pars,
                         S,...)

  if(is.na(weight_proposal))weight_proposal = -Inf

  # greedy acceptance rule
  if(is.finite(weight_proposal)) {
    current_pars[chain_index,] <- pars_use
    current_weight[chain_index] <- weight_proposal
  }

  return(c(current_weight[chain_index],current_pars[chain_index,]))

}
