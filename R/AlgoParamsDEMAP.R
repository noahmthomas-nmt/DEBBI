#' AlgoParamsDEMAP
#' @description get control parameters for DEMAP function
#' @param n_params number of free parameters estimated
#' @param n_chains number of particle chains, 3*n_params is the default value
#' @param n_iter number of iterations to run the sampling algorithm, 1000 is default
#' @param crossover_rate number on the interval (0,1]. Determines the probability a parameter on a chain is updated on a given crossover step, sampled from a Bernoulli distribution.
#' @param init_sd positive scalar or n_params-dimensional numeric vector, determines the standard deviation of the Gaussian initialization distribution
#' @param init_center scalar or n_params-dimensional numeric vector that determines the mean of the Gaussian initialization distribution
#' @param n_cores_use number of cores used when using parallelization.
#' @param step_size positive scalar, jump size in DE crossover step, default is 2.38/sqrt(2*n_params).
#' @param jitter_size positive scalar, noise is added during crossover step from Uniform(-jitter_size,jitter_size) distribution. 1e-6 is the default value.
#' @param parallel_type string specifying parallelization type. 'none','FORK', or 'PSOCK' are valid values. 'none' is default value.
#' @param return_trace logical, if true, function returns particle trajectories. This is helpful for diagnosing convergence or debugging model code. Function will return an iteration/thin $x$ n_chains $x$ n_params array and the estimated ELBO of each particle in a iteration/thin x n_chains array.
#' @param thin positive integer, only every 'thin'-th iteration will be stored. Default value is 1. Increasing thin will reduce the memory required, while running chains for longer.
#' @return list of control parameters for the DEMAP function
#' @export

AlgoParamsDEMAP <- function(n_params,
                            n_chains = NULL,
                            n_iter = 1000,
                            init_sd = 0.01,
                            init_center = 0,
                            n_cores_use = 1,
                            step_size = NULL,
                            jitter_size = 1e-6,
                            crossover_rate = 1,
                            parallel_type = "none",
                            return_trace = FALSE,
                            thin = 1) {
  # n_params
  ### catch errors
  n_params <- as.integer(n_params)
  if (any(!is.finite(n_params))) {
    stop("ERROR: n_params is not finite")
  } else if (n_params < 1 | length(n_params) > 1) {
    stop("ERROR: n_params must be a postitive integer scalar")
  }

  # n_chains
  ### if null assign default value
  if (is.null(n_chains)) {
    n_chains <- max(3 * n_params, 4)
  }
  ### catch errors
  n_chains <- as.integer(n_chains)
  if (any(!is.finite(n_chains))) {
    stop("ERROR: n_chains is not finite")
  } else if (n_chains < 4 | length(n_chains) > 1) {
    stop("ERROR: n_chains must be a postitive integer scalar, and atleast 4")
  }

  # n_iter
  ### if null assign default value
  if (is.null(n_iter)) {
    n_iter <- 1000
  }
  ### catch errors
  n_iter <- as.integer(n_iter)
  if (any(!is.finite(n_iter))) {
    stop("ERROR: n_iter is not finite")
  } else if (n_iter < 4 | length(n_iter) > 1) {
    stop("ERROR: n_iter must be a postitive integer scalar, and atleast 4")
  }

  # init_sd
  init_sd <- as.numeric(init_sd)
  if (any(!is.finite(init_sd))) {
    stop("ERROR: init_sd is not finite")
  } else if (any(init_sd <= 0 | is.complex(init_sd))) {
    stop("ERROR: init_sd must be positive and real-valued")
  } else if (!(length(init_sd) == 1 | length(init_sd) == n_params)) {
    stop("ERROR: init_sd vector length must be 1 or n_params")
  }

  # init_center
  init_center <- as.numeric(init_center)
  if (any(!is.finite(init_center))) {
    stop("ERROR: init_center is not finite")
  } else if (any(is.complex(init_center))) {
    stop("ERROR: init_center must be real valued")
  } else if (!(length(init_center) == 1 | length(init_center) == n_params)) {
    stop("ERROR: init_center vector length must be 1 or n_params")
  }

  # n_cores_use
  ### assign NULL value default
  if (is.null(n_cores_use)) {
    n_cores_use <- 1
  }
  ### catch any errors
  n_cores_use <- as.integer(n_cores_use)
  if (any(!is.finite(n_cores_use))) {
    stop("ERROR: n_cores_use is not finite")
  } else if (n_cores_use < 1 | length(n_cores_use) > 1) {
    stop("ERROR: n_cores_use must be a postitive integer scalar, and atleast 1")
  }


  # step_size
  ### assign NULL value default
  if (is.null(step_size)) {
    step_size <- 2.38 / sqrt(2 * n_params) # step size recommend in ter braak's 2006 paper
  }
  ### catch any errors
  if (any(!is.finite(step_size))) {
    stop("ERROR: step_size is not finite")
  } else if (any(step_size <= 0 | is.complex(step_size))) {
    stop("ERROR: step_size must be positive and real-valued")
  } else if (!(length(step_size) == 1)) {
    stop("ERROR: step_size vector length must be 1 ")
  }

  # jitter_size
  ### assign NULL value default
  if (is.null(jitter_size)) {
    jitter_size <- 1e-6
  }
  ### catch any errors
  if (any(!is.finite(jitter_size))) {
    stop("ERROR: jitter_size is not finite")
  } else if (any(jitter_size <= 0 | is.complex(jitter_size))) {
    stop("ERROR: jitter_size must be positive and real-valued")
  } else if (!(length(jitter_size) == 1)) {
    stop("ERROR: jitter_size vector length must be 1 ")
  }

  # crossover_rate
  ### if null assign default value
  if (any(is.null(crossover_rate))) {
    crossover_rate <- 1
  }
  ### catch errors
  crossover_rate <- as.numeric(crossover_rate)
  if (any(!is.finite(crossover_rate))) {
    stop("ERROR: crossover_rate is not finite")
  } else if (any(crossover_rate > 1) | any(crossover_rate <= 0) | length(crossover_rate) > 1) {
    stop("ERROR: crossover_rate must be a numeric scalar on the interval (0,1]")
  } else if (is.complex(crossover_rate)) {
    stop("ERROR: crossover_rate cannot be complex")
  }

  # parallel_type
  validParType <- c("none", "FORK", "PSOCK")
  ### assign NULL value default
  if (is.null(parallel_type)) {
    parallel_type <- "none"
  }
  ### catch any errors
  if (!parallel_type %in% validParType) {
    stop("ERROR: invalid parallel_type")
  }

  # thin
  ### if null assign default value
  if (is.null(thin)) {
    thin <- 0
  }
  ### catch errors
  thin <- as.integer(thin)
  if (any(!is.finite(thin))) {
    stop("ERROR: thin is not finite")
  } else if (any(thin < 1) | length(thin) > 1) {
    stop("ERROR: thin must be a scalar postive integer")
  }

  # nSamples Per Chains
  n_iters_per_chain <- floor((n_iter) / thin)
  ### catch errors
  if (n_iters_per_chain < 1 | (!is.finite(n_iters_per_chain))) {
    stop("ERROR: number of samples per chain is negative or non finite.
         n_iters_per_chain=floor((n_iter-burnin)/thin)")
  }

  # purify
  n_iters_per_chain <- floor((n_iter) / thin)
  ### catch errors
  if (n_iters_per_chain < 1 | (!is.finite(n_iters_per_chain))) {
    stop("ERROR: number of samples per chain is negative or non finite.
         n_iters_per_chain=floor((n_iter-burnin)/thin)")
  }


  out <- list(
    "n_params" = n_params,
    "n_chains" = n_chains,
    "n_iter" = n_iter,
    "init_sd" = init_sd,
    "init_center" = init_center,
    "n_cores_use" = n_cores_use,
    "step_size" = step_size,
    "crossover_rate" = crossover_rate,
    "jitter_size" = jitter_size,
    "parallel_type" = parallel_type,
    "thin" = thin,
    "purify" = Inf,
    "n_iters_per_chain" = n_iters_per_chain,
    "return_trace" = return_trace
  )

  return(out)
}
