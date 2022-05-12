#' DEVI
#' @description DE optimization for mean-field variational inference. Minimizes the KL divergence (maximizes the ELBO) between $q(theta|lambda)$ and the target posterior $p(theta|data)$ For a tutorial on variational inference check out Galdo, Bahg, & Turner 2020.
#' @param LogPostLike function whose first argument is an n_params-dimensional model parameter vector and returns (scalar) sum of log prior density and log likelihood for the parameter vector.
#' @param control_params control parameters for DE algorithm. see \code{\link{AlgoParamsDEVI}} function documentation for more details.
#' @param ... additional arguments to pass LogPostLike
#' @return list contain mean in a n_iters_per_chain by n_chains by 2*n_params_model array and the ELBO of each sample in a n_iters_per_chain by n_chains array.
#' @export
#' @md
#' @examples
#' # simulate from model
#' dataExample <- matrix(stats::rnorm(100, c(-1, 1), c(1, 1)), nrow = 50, ncol = 2, byrow = TRUE)
#' ## list parameter names
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
#' # Get variational approximation
#' DEVI(
#'   LogPostLike = LogPostLikeExample,
#'   control_params = AlgoParamsDEVI(
#'     n_params = length(param_names_example),
#'     n_iter = 200,
#'     n_chains = 12
#'   ),
#'   data = dataExample,
#'   param_names = param_names_example
#' )
#'
DEVI <- function(LogPostLike, control_params = AlgoParamsDEVI(), ...) {

  # import values we will reuse throughout process
  # create memory structures for storing posterior samples
  lambda <- array(NA, dim = c(
    control_params$n_iters_per_chain,
    control_params$n_chains,
    control_params$n_params_dist
  ))
  ELBO_values <- matrix(-Inf,
    nrow = control_params$n_iters_per_chain,
    ncol = control_params$n_chains
  )


  # chain initialization
  message("initalizing chains...")
  for (chain_idx in 1:control_params$n_chains) {
    count <- 0
    while (ELBO_values[1, chain_idx] == -Inf) {
      lambda[1, chain_idx, ] <- stats::rnorm(
        control_params$n_params_dist,
        control_params$init_center,
        control_params$init_sd
      )

      ELBO_values[1, chain_idx] <- ELBO(lambda[1, chain_idx, ], LogPostLike, control_params, S = control_params$n_samples_ELBO, ...)
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

  message("running DE to find best variational approximation")
  lambdaIdx <- 1
  for (iter in 1:control_params$n_iter) {

    #####################
    ####### crossover
    #####################
    if (control_params$parallel_type == "none") {
      temp <- matrix(unlist(lapply(1:control_params$n_chains, CrossoverVI,
        current_params = lambda[lambdaIdx, , ], # current parameter values for chain (numeric vector)
        current_weight = ELBO_values[lambdaIdx, ], # corresponding log like for (numeric vector)
        LogPostLike = LogPostLike, # log likelihood function (returns scalar)
        step_size = control_params$step_size,
        jitter_size = control_params$jitter_size,
        n_chains = control_params$n_chains,
        crossover_rate = control_params$crossover_rate,
        control_params = control_params,
        S = control_params$n_samples_ELBO, ...
      )),
      nrow = control_params$n_chains,
      ncol = control_params$n_params_dist + 1, byrow = TRUE
      )
    } else {
      temp <- matrix(unlist(parallel::parLapplyLB(cl_use, 1:control_params$n_chains, CrossoverVI,
        current_params = lambda[lambdaIdx, , ], # current parameter values for chain (numeric vector)
        current_weight = ELBO_values[lambdaIdx, ], # corresponding log like for (numeric vector)
        LogPostLike = LogPostLike, # log likelihood function (returns scalar)
        step_size = control_params$step_size,
        jitter_size = control_params$jitter_size,
        n_chains = control_params$n_chains,
        crossover_rate = control_params$crossover_rate,
        control_params = control_params,
        S = control_params$n_samples_ELBO, ...
      )),
      control_params$n_chains,
      control_params$n_params_dist + 1,
      byrow = TRUE
      )
    }
    # update particle chains
    ELBO_values[lambdaIdx, ] <- temp[, 1]
    lambda[lambdaIdx, , ] <- temp[, 2:(control_params$n_params_dist + 1)]
    if (iter < control_params$n_iter) {
      ELBO_values[lambdaIdx + 1, ] <- temp[, 1]
      lambda[lambdaIdx + 1, , ] <- temp[, 2:(control_params$n_params_dist + 1)]
    }

    #####################
    ####### purify
    #####################
    if (iter %% control_params$purify == 0) {
      if (control_params$parallel_type == "none") {
        temp <- matrix(unlist(lapply(1:control_params$n_chains, PurifyVI,
          current_params = lambda[lambdaIdx, , ], # current parameter values for chain (numeric vector)
          current_weight = ELBO_values[lambdaIdx, ], # corresponding log like for (numeric vector)
          LogPostLike = LogPostLike, # log likelihood function (returns scalar)
          n_chains = control_params$n_chains,
          control_params = control_params,
          S = control_params$n_samples_ELBO, ...
        )),
        nrow = control_params$n_chains,
        ncol = control_params$n_params_dist + 1, byrow = TRUE
        )
      } else {
        temp <- matrix(unlist(parallel::parLapplyLB(cl_use, 1:control_params$n_chains, PurifyVI,
          current_params = lambda[lambdaIdx, , ], # current parameter values for chain (numeric vector)
          current_weight = ELBO_values[lambdaIdx, ], # corresponding log like for (numeric vector)
          LogPostLike = LogPostLike, # log likelihood function (returns scalar)
          control_params = control_params,
          S = control_params$n_samples_ELBO, ...
        )),
        control_params$n_chains,
        control_params$n_params_dist + 1,
        byrow = T
        )
      }

      # update particle chains
      ELBO_values[lambdaIdx, ] <- temp[, 1]
      lambda[lambdaIdx, , ] <- temp[, 2:(control_params$n_params_dist + 1)]
      if (iter < control_params$n_iter) {
        ELBO_values[lambdaIdx + 1, ] <- temp[, 1]
        lambda[lambdaIdx + 1, , ] <- temp[, 2:(control_params$n_params_dist + 1)]
      }
    }

    if (iter %% 100 == 0) message(paste0("iter ", iter, "/", control_params$n_iter))
    if (iter %% control_params$thin == 0) {
      lambdaIdx <- lambdaIdx + 1
    }
  }
  # cluster stop
  if (!control_params$parallel_type == "none") {
    parallel::stopCluster(cl = cl_use)
  }
  maxIdx <- which.max(ELBO_values[control_params$n_iters_per_chain, ])

  means <- lambda[
    control_params$n_iters_per_chain,
    maxIdx, 1:control_params$n_params_model
  ]
  names(means) <- paste0(control_params$param_names, "_mean")

  covariance <- diag(exp(2 * lambda[
    control_params$n_iters_per_chain,
    maxIdx, (control_params$n_params_model + 1):control_params$n_params_dist
  ]))

  if (control_params$LRVB_correction == TRUE) {
    message("Attemtping LRVB covariance correction.")
    control_params_LRVB <- control_params
    control_params_LRVB$use_QMC <- TRUE
    hess <- numDeriv::hessian(KLHatforLRVB, c(means, diag(covariance)),
      method = "Richardson", method.args = list(),
      LogPostLike,
      control_params = control_params_LRVB,
      S = control_params$n_samples_LRVB, ...
    )
    if (any(round(eigen(hess)$values, 4) < 0)) {
      warning("LRVB correction failed, optimization did not reach the neighborhood of a local optima. Returning mean field approximation.")
    } else {
      message("LRVB correction was a success!")
      hess.inv <- solve(hess)
      covariance <- hess.inv[1:control_params$n_params_model, 1:control_params$n_params_model]
    }
  }

  if (control_params$return_trace == TRUE) {
    return(list(
      "means" = means,
      "covariance" = covariance,
      "ELBO" = ELBO_values[control_params$n_iters_per_chain, maxIdx],
      "lambda_trace" = lambda,
      "ELBO_trace" = ELBO_values
    ))
  } else {
    return(list(
      "means" = means,
      "covariance" = covariance,
      "ELBO" = ELBO_values[control_params$n_iters_per_chain, maxIdx]
    ))
  }
}
