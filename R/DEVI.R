#' DEVI
#' @description DE optimization for mean-field variational inference. Minimizes the KL divergence (maximizes the ELBO) between $q(theta|lambda)$ and the target posterior $p(theta|data)$ For a tutorial on variational inference check out Galdo, Bahg, & Turner 2020.
#' @param LogPostLike function whose first arguement is an n_pars-dimensional model parameter vector and returns (scalar) sum of log prior density and log likelihood for the parameter vector.
#' @param control_pars control parameters for DE algo. see \code{\link{AlgoParsDEVI}} function documentation for more details.
#' @param ... additional arguments to pass LogPostLike
#' @return list contain mean in a n_iters_per_chain $x$ n_chains $x$ 2*n_pars_model array and the ELBO of each sample in a n_iters_per_chain x n_chains array.
#' @export
#' @md
#' @examples
#' #simulate from model


#' dataExample=matrix(stats::rnorm(200,c(-1,1),c(.1,.1)),nrow=50,ncol=2,byrow = TRUE)
#' ##list parameter names
#' par_names_example=c("mu_1","mu_2")
#'
#' #log posterior likelihood function = log likelihood + log prior | returns a scalar
#' LogPostLikeExample=function(x,data,par_names){
#'   out=0
#'
#'   names(x)<-par_names
#'
#'   # log prior
#'   out=out+sum(dnorm(x["mu_1"],0,sd=1,log=TRUE))
#'   out=out+sum(dnorm(x["mu_2"],0,sd=1,log=TRUE))
#'
#'  # log likelihoods
#'   out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=TRUE))
#'
#'   return(out)
#' }
#'
#'# Get variational approximation
#'DEVI(LogPostLike=LogPostLikeExample,
#'       control_pars=AlgoParsDEVI(n_pars=length(par_names_example),
#'                                   n_iter=200,
#'                                   n_chains=12),
#'                                   data=dataExample,
#'                                   par_names = par_names_example)



DEVI=function(LogPostLike,control_pars=AlgoParsDEVI(),...){

  # import values we will reuse throughout process
  # create memory structures for storing posterior samples
  lambda=array(NA,dim=c(control_pars$n_iters_per_chain,
                        control_pars$n_chains,
                        control_pars$n_pars_dist))
  ELBO_values=matrix(-Inf,
                     nrow=control_pars$n_iters_per_chain,
                     ncol=control_pars$n_chains)


  # chain initialization
  print('initalizing chains...')
  for(chain_idx in 1:control_pars$n_chains){
    count=0
    while (ELBO_values[1,chain_idx]  == -Inf) {
      lambda[1,chain_idx,] <- stats::rnorm(control_pars$n_pars_dist,
                                           control_pars$init_center,
                                           control_pars$init_sd)

      ELBO_values[1,chain_idx] <- ELBO(lambda[1,chain_idx,],LogPostLike,control_pars,S=control_pars$n_samples_ELBO,...)
      count=count+1
      if(count>100){
        stop('chain initialization failed.
        inspect likelihood and prior or change init_center/init_sd to sample more
             likely parameter values')
      }
    }
    print(paste0(chain_idx," / ",control_pars$n_chains))
  }
  print('chain initialization complete  :)')

  # cluster initialization
  if(!control_pars$parallel_type=='none'){

    print(paste0("initalizing ",
                 control_pars$parallel_type," cluser with ",
                 control_pars$n_cores_use," cores"))

    doParallel::registerDoParallel(control_pars$n_cores_use)
    cl_use <- parallel::makeCluster(control_pars$n_cores_use,
                                    type=control_pars$parallel_type)
  }

  print("running DE to find best variational approximation")
  lambdaIdx=1
  for(iter in 1:control_pars$n_iter){

    #####################
    ####### crossover
    #####################
    if(control_pars$parallel_type=='none'){
      temp=matrix(unlist(lapply(1:control_pars$n_chains,CrossoverVI,
                                current_pars=lambda[lambdaIdx,,],  # current parameter values for chain (numeric vector)
                                current_weight=ELBO_values[lambdaIdx,], # corresponding log like for (numeric vector)
                                LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                step_size=control_pars$step_size,
                                jitter_size=control_pars$jitter_size,
                                n_chains=control_pars$n_chains,
                                crossover_rate=control_pars$crossover_rate,
                                control_pars=control_pars,
                                S=control_pars$n_samples_ELBO,...)),
                  nrow=control_pars$n_chains,
                  ncol=control_pars$n_pars_dist+1,byrow=T)
    } else {
      temp=matrix(unlist(parallel::parLapply(cl_use,1:control_pars$n_chains,CrossoverVI,
                                             current_pars=lambda[lambdaIdx,,],  # current parameter values for chain (numeric vector)
                                             current_weight=ELBO_values[lambdaIdx,], # corresponding log like for (numeric vector)
                                             LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                             step_size=control_pars$step_size,
                                             jitter_size=control_pars$jitter_size,
                                             n_chains=control_pars$n_chains,
                                             crossover_rate=control_pars$crossover_rate,
                                             control_pars,
                                             S=control_pars$n_samples_ELBO,...)),
                  control_pars$n_chains,
                  control_pars$n_pars_dist+1,byrow=T)

    }
    # update particle chains
    ELBO_values[lambdaIdx,]=temp[,1]
    lambda[lambdaIdx,,]=temp[,2:(control_pars$n_pars_dist+1)]
    if(iter<control_pars$n_iter){
      ELBO_values[lambdaIdx+1,]=temp[,1]
      lambda[lambdaIdx+1,,]=temp[,2:(control_pars$n_pars_dist+1)]
    }

    #####################
    ####### purify
    #####################
    if(iter%%control_pars$purify==0){

      if(control_pars$parallel_type=='none'){
        temp=matrix(unlist(lapply(1:control_pars$n_chains,PurifyVI,
                                  current_pars=lambda[lambdaIdx,,],  # current parameter values for chain (numeric vector)
                                  current_weight=ELBO_values[lambdaIdx,], # corresponding log like for (numeric vector)
                                  LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                  n_chains=control_pars$n_chains,
                                  control_pars=control_pars,
                                  S=control_pars$n_samples_ELBO,...)),
                    nrow=control_pars$n_chains,
                    ncol=control_pars$n_pars_dist+1,byrow=T)
      } else {
        temp=matrix(unlist(parallel::parLapply(cl_use,1:control_pars$n_chains,PurifyVI,
                                               current_pars=lambda[lambdaIdx,,],  # current parameter values for chain (numeric vector)
                                               current_weight=ELBO_values[lambdaIdx,], # corresponding log like for (numeric vector)
                                               LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                               step_size=control_pars$step_size,
                                               jitter_size=control_pars$jitter_size,
                                               n_chains=control_pars$n_chains,
                                               crossover_rate=control_pars$crossover_rate,
                                               control_pars,
                                               S=control_pars$n_samples_ELBO,...)),
                    control_pars$n_chains,
                    control_pars$n_pars_dist+1,byrow=T)

      }

      # update particle chains
      ELBO_values[lambdaIdx,]=temp[,1]
      lambda[lambdaIdx,,]=temp[,2:(control_pars$n_pars_dist+1)]
      if(iter<control_pars$n_iter){
        ELBO_values[lambdaIdx+1,]=temp[,1]
        lambda[lambdaIdx+1,,]=temp[,2:(control_pars$n_pars_dist+1)]
      }
    }




    if(iter%%100==0)print(paste0('iter ',iter,'/',control_pars$n_iter))
    if(iter%%control_pars$thin==0){
      lambdaIdx=lambdaIdx+1
    }

  }
  # cluster stop
  if(!control_pars$parallel_type=='none'){
    parallel::stopCluster(cl=cl_use)
  }
  maxIdx=which.max(ELBO_values[control_pars$n_iters_per_chain,])

  means=lambda[control_pars$n_iters_per_chain,
               maxIdx,1:control_pars$n_pars_model]
  names(means)<-paste0(control_pars$par_names,"_mean")

  covariance=diag(exp(2*lambda[control_pars$n_iters_per_chain,
                               maxIdx,(control_pars$n_pars_model+1):control_pars$n_pars_dist]))

  if(control_pars$LRVB_correction==T){
    print("Attemtping LRVB covariance correction.")
    control_pars_LRVB=control_pars
    control_pars_LRVB$use_QMC <- TRUE
    hess=numDeriv::hessian(KLHatforLRVB,c(means,diag(covariance)),
                                          method="Richardson", method.args=list(),
                                          control_pars=control_pars_LRVB,
                                          S=control_pars$n_samples_LRVB,...)
    if(any(round(eigen(hess)$values,4))<0){
      warning("LRVB correction failed, optimization did not reach the neighborhood of a local optima. Returning mean field approximation.")
    } else {
      hess.inv=solve(hess)
      covariance <- hess.inv[1:n.pars.model,1:n.pars.model]
    }
  }
  if(control_pars$return_trace==T){
    return(list('means'=means,
                'covariance'=covariance,
                'ELBO'=ELBO_values[control_pars$n_iters_per_chain,maxIdx],
                'lambda_trace'=lambda,
                'ELBO_trace'=ELBO_values,
                'control_pars'=control_pars))
  } else {
    return(list('means'=means,
                'covariance'=covariance,
                'ELBO'=ELBO_values[control_pars$n_iters_per_chain,maxIdx],
                'control_pars'=control_pars))
  }
}
