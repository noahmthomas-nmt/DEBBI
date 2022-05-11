#' AlgoParamsDEVI
#' @description get control parameters for DEVI function
#' @param n_params number of free parameters estimated
#' @param param_names optional vector of parameter names
#' @param n_chains number of particle chains used for optimization, 3*n_params is the default value
#' @param n_iter number of iterations to run the sampling algorithm, 1000 is default
#' @param crossover_rate number on the interval (0,1]. Determines the probability a parameter on a chain is updated on a given crossover step, sampled from a Bernoulli distribution.
#' @param init_sd positive scalar or n_params-dimensional numeric vector, determines the standard deviation of the Gaussian initialization distribution
#' @param init_center scalar or n_params-dimensional numeric vector, determines the mean of the Gaussian initialization distribution
#' @param n_cores_use number of cores used when using parallelization.
#' @param step_size positive scalar, jump size in DE crossover step, default is 2.38/sqrt(2*n_params).
#' @param jitter_size positive scalar, noise is added during crossover step from Uniform(-jitter_size,jitter_size) distribution. 1e-6 is the default value.
#' @param parallel_type string specifying parallelization type. 'none','FORK', or 'PSOCK' are valid values. 'none' is default value.
#' @param return_trace logical, if true, function returns particle trajectories. This is helpful for diagnosing convergence or debugging model code. Function will return an iteration/thin $x$ n_chains $x$ n_params array and the estimated ELBO of each particle in a iteration/thin x n_chains array.
#' @param thin positive integer, only every 'thin'-th iteration will be stored. Default value is 1. Increasing thin will reduce the memory required, while running algorithm for longer.
#' @param burnin number of initial iterations to discard. Default value is 0.
#' @param purify an integer, every 'purify'-th iteration, the Monte Carlo estimator of the ELBO is recalculated. This can help deal with noisy and outlier estimates of the ELBO. Default value is 25. If use_QMC is TRUE, purification is disabled as it is redundant.
#' @param n_samples_ELBO number of samples used for the Monte Carlo estimator of the ELBO (the objective function). default is 10.
#' @param use_QMC logical, if true, a quasi-Monte Carlo estimator is used to estimate ELBO during optimization. default is TRUE.
#' @param LRVB_correction logical, if true, LRVB covariance correction (Giordano, Brodderick, & Jordan 2018; Galdo, Bahg, & Turner 2020) is attempted.
#' @param n_samples_LRVB number of samples used for LRVB correction. default is 25.
#' @param quasi_rand_seq type of low discrepancy sequence used for quasi Monte Carlo integration, either 'sobol' or 'halton'. LRVB correction always use QMC. Default is 'sobol'.
#' @param neg_inf if density for a given value of theta is numerically 0 for q, this value is assigned for log density. This helps with numeric stability of algorithm. Default value is -750.
#' @return list of control parameters for the DEVI function
#' @export
#'

AlgoParamsDEVI <- function(n_params,
                           param_names = NULL,
                           n_chains = NULL,
                           n_iter = 1000,
                           init_sd = 0.01,
                           init_center = 0,
                           n_cores_use = 1,
                           step_size = NULL,
                           jitter_size = 1e-6,
                           parallel_type = "none",
                           use_QMC = TRUE,
                           purify = NULL,
                           quasi_rand_seq = "halton",
                           n_samples_ELBO = 10,
                           LRVB_correction = TRUE,
                           n_samples_LRVB = 25,
                           neg_inf = -750,
                           thin = 1,
                           burnin = 0,
                           return_trace = FALSE,
                           crossover_rate = 1) {
  # n_params
  ### catch errors
  n_params <- as.integer(n_params)
  if (any(!is.finite(n_params))) {
    stop("ERROR: n_params is not finite")
  } else if (n_params < 1 | length(n_params) > 1) {
    stop("ERROR: n_params must be a postitive integer scalar")
  }


  # param_names
  ### catch errors
  if (is.null(param_names)) {
    param_names <- paste0("param_", 1:n_params)
  } else if (!(length(param_names) == n_params)) {
    stop("ERROR: param_names does not match size of n_params")
  }

  dist_param_names <- c(paste0(param_names, "_MEAN"), paste0(param_names, "_VAR"))

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
  } else if (any(init_sd < 0 | is.complex(init_sd))) {
    stop("ERROR: init_sd must be positive and real-valued")
  } else if (!(length(init_sd) == 1 | length(init_sd) == (n_params * 2))) {
    stop("ERROR: init_sd vector length must be 1 or n_params")
  }
  if (any(init_sd == 0)) {
    warning("WARNING an init_sd value is 0")
  }

  # init_center
  init_center <- as.numeric(init_center)
  if (any(!is.finite(init_center))) {
    stop("ERROR: init_center is not finite")
  } else if (any(is.complex(init_center))) {
    stop("ERROR: init_center must be real valued")
  } else if (!(length(init_center) == 1 | length(init_center) == (n_params * 2))) {
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
    step_size <- 2.38 / sqrt(4 * n_params) # step size recommend in ter braak's 2006 paper
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


  # use_QMC
  if (any(is.null(use_QMC))) {
    use_QMC <- TRUE
  }
  if (length(use_QMC) > 1) {
    stop("length(use_QMC)>1, please use a scalar logical")
  }
  if (!((use_QMC == 0) | (use_QMC == 1))) {
    stop("ERROR: use_QMC must be a scalar logical")
  }
  use_QMC <- as.logical(use_QMC)

  # LRVB_correction
  ### assign value if null
  if (any(is.null(LRVB_correction))) {
    LRVB_correction <- TRUE
  }
  ### catch errors
  if (length(LRVB_correction) > 1) {
    stop("length(LRVB_correction)>1, please use a scalar logical")
  }
  if (!((LRVB_correction == 0) | (LRVB_correction == 1))) {
    stop("ERROR: LRVB_correction must be a scalar logical")
  }
  LRVB_correction <- as.logical(LRVB_correction)

  if (LRVB_correction) {
    # if using LRVB correction check for valid samples count
    # n_samples_LRVB
    ### assign value if null
    n_samples_LRVB <- as.integer(n_samples_LRVB)
    if (any(is.null(n_samples_LRVB))) {
      n_samples_LRVB <- 25
    }
    ### catch errors
    if (any(!is.finite(n_samples_LRVB))) {
      stop("ERROR: n_samples_LRVB is not finite")
    } else if (any(n_samples_LRVB < 1) | length(n_samples_LRVB) > 1) {
      stop("ERROR: n_samples_LRVB must be a scalar postive integer")
    }
  }

  # n_samples_ELBO
  ### assign value if null
  n_samples_ELBO <- as.integer(n_samples_ELBO)
  if (any(is.null(n_samples_ELBO))) {
    n_samples_ELBO <- 10
  }
  ### catch errors
  if (any(!is.finite(n_samples_ELBO))) {
    stop("ERROR: n_samples_ELBO is not finite")
  } else if (any(n_samples_ELBO < 1) | length(n_samples_ELBO) > 1) {
    stop("ERROR: n_samples_ELBO must be a scalar postive integer")
  }

  # quasi_rand_seq
  quasi_rand_seq <- tolower(as.character(quasi_rand_seq))
  valid_quasi_rand_seqs <- c("sobol", "halton")
  ### assign NULL value default
  if (is.null(quasi_rand_seq)) {
    quasi_rand_seq <- "sobol"
  }
  ### catch any errors
  if (!quasi_rand_seq %in% valid_quasi_rand_seqs) {
    stop(paste0("ERROR: invalid quasi_rand_seq; must be one of: ", paste(valid_quasi_rand_seqs, sep = ",")))
  }

  # neg_inf
  ### assign NULL value default
  if (any(is.null(neg_inf))) {
    neg_inf <- -750
  }
  ### catch any errors
  if (length(neg_inf) > 1) {
    stop("length(neg_inf)>1, please use a scalar numeric")
  }
  if (!is.numeric(neg_inf)) {
    stop("neg_inf must be a numeric")
  }

  # burnin
  ### if null assign default value
  if (is.null(burnin)) {
    burnin <- 0
  }
  ### catch errors
  burnin <- as.integer(burnin)
  if (any(!is.finite(burnin))) {
    stop("ERROR: burnin is not finite")
  } else if (any(burnin < 0) | any(burnin >= n_iter) | length(burnin) > 1) {
    stop("ERROR: burnin must be a scalar integer from the interval [0,n_iter)")
  }

  # nSamples Per Chains
  n_iters_per_chain <- floor((n_iter - burnin) / thin)
  ### catch errors
  if (n_iters_per_chain < 1 | (!is.finite(n_iters_per_chain))) {
    stop("ERROR: number of iters per chain is negative or non finite.
           n_iters_per_chain=floor((n_iter-burnin)/thin)")
  }

  # purify
  ### default values
  purify <- as.integer(purify)
  if (any(is.null(purify))) {
    purify <- 25
  }
  if (use_QMC == TRUE) {
    purify <- n_iter + 1
  }
  ### catch errors
  if (any(!is.finite(purify))) {
    stop("ERROR: purify is not finite")
  } else if (purify < 1 | length(purify) > 1) {
    stop("ERROR: purify must be a postitive integer scalar")
  }

  out <- list(
    "n_params_model" = n_params,
    "param_names" = param_names,
    "n_chains" = n_chains,
    "n_iter" = n_iter,
    "init_sd" = init_sd,
    "init_center" = init_center,
    "n_cores_use" = n_cores_use,
    "step_size" = step_size,
    "jitter_size" = jitter_size,
    "crossover_rate" = crossover_rate,
    "parallel_type" = parallel_type,
    "return_trace" = return_trace,
    "purify" = purify,
    "use_QMC" = use_QMC,
    "quasi_rand_seq" = quasi_rand_seq,
    "n_samples_ELBO" = n_samples_ELBO,
    "n_samples_LRVB" = n_samples_LRVB,
    "LRVB_correction" = LRVB_correction,
    "thin" = thin,
    "neg_inf" = neg_inf,
    "n_params_dist" = 2 * n_params,
    "n_iters_per_chain" = n_iters_per_chain
  )

  return(out)
}
