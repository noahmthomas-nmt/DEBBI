#' CrossoverOptimize
#'
#' @param chain_index  which chain you are updating
#' @param current_pars current parameter values for chain (numeric vector)
#' @param current_weight  weight for current chain values
#' @param objFun  function you want to minimize (returns scalar)
#' @param step_size step size in DE jump
#' @param jitter_size noise
#' @param n_chains number of chains
#' @param crossover_rate rate of updating parameters on a give crossover step (0=no parameters updated, 1=all parameters updated)
#' @param ... additional arguments for objFun function
#' @noRd



CrossoverOptimize=function(chain_index,
                           current_pars,
                           current_weight,
                           objFun,
                           step_size=.8,
                           jitter_size=1e-6,
                           n_chains,crossover_rate=1, ... ){

  # get statistics about chain
  weight_use = current_weight[chain_index]
  pars_use = current_pars[chain_index,]
  len_par_use = length(pars_use)

  # use binomial to sample which pars to update matching crossover rate frequency
  par_idices_bool = stats::rbinom(len_par_use, prob = crossover_rate, size=1)

  # if no pars selected, randomly sample 1 parameter
  if(all(par_idices_bool==0)){
    par_idices_bool[sample(x=1:len_par_use,size=1)] <- 1
  }

  # indices of parameters to be updated
  par_indices=seq(1,len_par_use,by=1)[as.logical(par_idices_bool)]

  # sample parent chains
  parent_chain_indices = sample(c(1:n_chains)[-chain_index],2,replace=F)

  # mate parents for proposal
  pars_use[par_indices] = pars_use[par_indices] +
    step_size*(current_pars[parent_chain_indices[1],par_indices] -
                 current_pars[parent_chain_indices[2],par_indices]) +
    stats::runif(1,-jitter_size,jitter_size)

  pars_use = matrix(pars_use,1,len_par_use)

  # get log weight
  weight_proposal = objFun(pars_use,...)
  if(is.na(weight_proposal))weight_proposal = -Inf

  # greedy acceptance rule
  if(weight_proposal > weight_use) {
    current_pars[chain_index,] <- pars_use
    current_weight[chain_index] <- weight_proposal
  }

  return(c(current_weight[chain_index],current_pars[chain_index,]))

}
