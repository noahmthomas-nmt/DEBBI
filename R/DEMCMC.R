#' DEMCMC
#'
#' @description Sample from posterior using Differential Evolution Markov Chain Monte Carlo
#' @param LogPostLike function whose first argument is an n_params-dimensional model parameter vector and returns (scalar) sum of log prior density and log likelihood for the parameter vector.
#' @param control_params control parameters for DEMCMC algorithm. see \code{\link{AlgoParamsDEMCMC}} function documentation for more details. You must specify 'n_params' here.
#' @param ... additional arguments to pass LogPostLike
#' @return list contain posterior samples from DEMCMC in a 'n_samples_per_chain' by 'n_chains' by n_params array and the log posterior likelihood of each sample in a 'n_samples_per_chain' by 'n_chains' array.
#' @export
#' @md
#' @examples
#' # simulate from model
#' dataExample <- matrix(stats::rnorm(100, c(-1, 1), c(1, 1)), nrow = 50, ncol = 2, byrow = TRUE)
#' #
#' # list parameter names
#' param_names_example <- c("mu_1", "mu_2")
#'
#' # log posterior likelihood function = log likelihood + log prior | returns a scalar
#' LogPostLikeExample <- function(x, data, param_names) {
#'   out <- 0
#'
#'   names(x) <- param_names
#'
#'   # log prior
#'   out <- out + sum(dnorm(x["mu_1"], 0, sd = 1, log = TRUE))
#'   out <- out + sum(dnorm(x["mu_2"], 0, sd = 1, log = TRUE))
#'
#'   # log likelihoods
#'   out <- out + sum(dnorm(data[, 1], x["mu_1"], sd = 1, log = TRUE))
#'   out <- out + sum(dnorm(data[, 2], x["mu_2"], sd = 1, log = TRUE))
#'
#'   return(out)
#' }
#'
#' # Sample from posterior
#' DEMCMC(
#'   LogPostLike = LogPostLikeExample,
#'   control_params = AlgoParamsDEMCMC(
#'     n_params = length(param_names_example),
#'     n_iter = 1000,
#'     n_chains = 12
#'   ),
#'   data = dataExample,
#'   param_names = param_names_example
#' )
DEMCMC <- function(LogPostLike, control_params = AlgoParamsDEMCMC(), ...) {

  # import values we will reuse throughout process
  # create memory structures for storing posterior samples
  theta <- array(NA, dim = c(control_params$n_samples_per_chain, control_params$n_chains, control_params$n_params))
  log_post_like <- matrix(-Inf, nrow = control_params$n_samples_per_chain, ncol = control_params$n_chains)


  # chain initialization
  message("initalizing chains...")
  for (chain_idx in 1:control_params$n_chains) {
    count <- 0
    while (log_post_like[1, chain_idx] == -Inf) {
      theta[1, chain_idx, ] <- stats::rnorm(control_params$n_params, control_params$init_center, control_params$init_sd)

      log_post_like[1, chain_idx] <- LogPostLike(theta[1, chain_idx, ], ...)
      count <- count + 1
      if (count > 100) {
        stop("chain initialization failed.
        inspect likelihood and prior or change init_center/init_sd to sample more
             likely parameter values")
      }
    }
    message(paste0(chain_idx, " / ", control_params$n_chains))
  }
  message("chain initialization complete  :)")

  # cluster initialization
  if (!control_params$parallel_type == "none") {
    message(paste0(
      "initalizing ",
      control_params$parallel_type, " cluser with ",
      control_params$n_cores_use, " cores"
    ))

    doParallel::registerDoParallel(control_params$n_cores_use)
    cl_use <- parallel::makeCluster(control_params$n_cores_use,
                                    type = control_params$parallel_type
    )
  }

  message("running DEMCMC")
  thetaIdx <- 1
  for (iter in 1:control_params$n_iter) {
    if (control_params$parallel_type == "none") {
      # cross over step
      temp <- matrix(unlist(lapply(1:control_params$n_chains, CrossoverMC,
                                   param_indices = 1:control_params$n_params,
                                   current_theta = theta[thetaIdx, , ], # current parameter values for chain (numeric vector)
                                   current_log_post_like = log_post_like[thetaIdx, ], # corresponding log like for (numeric vector)
                                   LogPostLike = LogPostLike, # log likelihood function (returns scalar)
                                   step_size = control_params$step_size,
                                   jitter_size = control_params$jitter_size,
                                   n_chains = control_params$n_chains, ...
      )), control_params$n_chains, control_params$n_params + 1, byrow = TRUE)
    } else {
      temp <- matrix(unlist(parallel::parLapplyLB(cl_use, 1:control_params$n_chains, CrossoverMC,
                                                  param_indices = 1:control_params$n_params,
                                                  current_theta = theta[thetaIdx, , ], # current parameter values for chain (numeric vector)
                                                  current_log_post_like = log_post_like[thetaIdx, ], # corresponding log like for (numeric vector)
                                                  LogPostLike = LogPostLike, # log likelihood function (returns scalar)
                                                  step_size = control_params$step_size,
                                                  jitter_size = control_params$jitter_size,
                                                  n_chains = control_params$n_chains, ...
      )), control_params$n_chains, control_params$n_params + 1, byrow = TRUE)
    }
    log_post_like[thetaIdx, ] <- temp[, 1]
    theta[thetaIdx, , ] <- temp[, 2:(control_params$n_params + 1)]
    if (iter < control_params$n_iter) {
      log_post_like[thetaIdx + 1, ] <- temp[, 1]
      theta[thetaIdx + 1, , ] <- temp[, 2:(control_params$n_params + 1)]
    }
    ########
    # PURIFICATION
    ########
    if (iter %% control_params$purify == 0) {
      if (control_params$parallel_type == "none") {
        # cross over step
        temp <- matrix(unlist(lapply(1:control_params$n_chains,
                                     PurifyMC,
                                     param_indices = 1:control_params$n_params,
                                     current_theta = theta[thetaIdx, , ], # current parameter values for chain (numeric vector)
                                     current_log_post_like = log_post_like[thetaIdx, ], # corresponding log like for (numeric vector)
                                     LogPostLike = LogPostLike, # log likelihood function (returns scalar)
                                     n_chains = control_params$n_chains, ...
        )), control_params$n_chains,
        control_params$n_params + 1,
        byrow = TRUE)
      } else {
        temp <- matrix(unlist(parallel::parLapplyLB(cl_use,
                                                    1:control_params$n_chains,
                                                    PurifyMC,
                                                    param_indices = 1:control_params$n_params,
                                                    current_theta = theta[thetaIdx, , ], # current parameter values for chain (numeric vector)
                                                    current_log_post_like = log_post_like[thetaIdx, ], # corresponding log like for (numeric vector)
                                                    LogPostLike = LogPostLike, # log likelihood function (returns scalar)
                                                    n_chains = control_params$n_chains, ...
        )),
        control_params$n_chains,
        control_params$n_params + 1,
        byrow = TRUE
        )
      }
      log_post_like[thetaIdx, ] <- temp[, 1]
      theta[thetaIdx, , ] <- temp[, 2:(control_params$n_params + 1)]
      if (iter < control_params$n_iter) {
        log_post_like[thetaIdx + 1, ] <- temp[, 1]
        theta[thetaIdx + 1, , ] <- temp[, 2:(control_params$n_params + 1)]
      }
    }
    if (iter %% 100 == 0) message(paste0("iter ", iter, "/", control_params$n_iter))
    if (iter > control_params$burnin & iter %% control_params$thin == 0) {
      thetaIdx <- thetaIdx + 1
    }
  }

  # cluster stop
  if (!control_params$parallel_type == "none") {
    parallel::stopCluster(cl = cl_use)
  }

  dimnames(theta)[[3]] <- control_params$param_names
  dimnames(theta)[[2]] <- paste0("Chain", 1:control_params$n_chains)
  dimnames(theta)[[1]] <- paste0("Sample", 1:control_params$n_samples_per_chain)

  return(list("samples" = theta, "log_post_like" = log_post_like))
}
