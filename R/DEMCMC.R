DEMCMC=function(LogPostLike,control_pars=AlgoParsDEMCMC(),...){

  # import values we will reuse throughout process
  # create memory structures for storing posterior samples
  theta=array(NA,dim=c(control_pars$nSamplesPerChain,control_pars$n_chains,control_pars$n_pars))
  log_post_like=matrix(-Inf,nrow=control_pars$nSamplesPerChain,ncol=control_pars$n_chains)


  # chain initialization
  print('initalizing chains...')
  for(chain_idx in 1:control_pars$n_chains){
    count=0
    while (log_post_like[1,chain_idx]  == -Inf) {
      theta[1,chain_idx,] <- stats::rnorm(control_pars$n_pars,control_pars$init_center,control_pars$init_sd)

      log_post_like[1,chain_idx] <- LogPostLike(theta[1,chain_idx,],...)
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
  if(!control_pars$parallelType=='none'){

    print(paste0("initalizing ",
                 control_pars$parallelType," cluser with ",
                 control_pars$n_cores_use," cores"))

    doParallel::registerDoParallel(control_pars$n_cores_use)
    cl_use <- parallel::makeCluster(control_pars$n_cores_use,
                          type=control_pars$parallelType)
  }

  print("running DEMCMC")
  thetaIdx=1
  for(iter in 1:control_pars$n_iter){

    if(control_pars$parallelType=='none'){
      # cross over step
      temp=matrix(unlist(lapply(1:control_pars$n_chains,CrossoverMC,
                                par_indices=1:control_pars$n_pars,
                                current_theta=theta[thetaIdx,,],  # current parameter values for chain (numeric vector)
                                current_log_post_like=log_post_like[thetaIdx,], # corresponding log like for (numeric vector)
                                LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                step_size=control_pars$step_size,
                                jitter_size=control_pars$jitter_size,
                                n_chains=control_pars$n_chains,...)),control_pars$n_chains,control_pars$n_pars+1,byrow=T)
    } else {
      temp=matrix(unlist(parallel::parLapply(cl_use,1:control_pars$n_chains,CrossoverMC,
                                   par_indices=1:control_pars$n_pars,
                                   current_theta=theta[thetaIdx,,],  # current parameter values for chain (numeric vector)
                                   current_log_post_like=log_post_like[thetaIdx,], # corresponding log like for (numeric vector)
                                   LogPostLike=LogPostLike, # log likelihood function (returns scalar)
                                   step_size=control_pars$step_size,
                                   jitter_size=control_pars$jitter_size,
                                   n_chains=control_pars$n_chains,...)),control_pars$n_chains,control_pars$n_pars+1,byrow=T)

    }
    log_post_like[thetaIdx,]=temp[,1]
    theta[thetaIdx,,]=temp[,2:(control_pars$n_pars+1)]

    #   # purification step
    #   if(iter%%purify.rate==0){
    #     temp=unlist(lapply(1:n_chains,PurifyMC,
    #                                 current_theta=theta[iter,],  # current parameter values for chain (numeric vector)
    #                                 current_log_post_like=log_post_like[iter,], # corresponding log like for (numeric vector)
    #                                 LogPostLike=LogPostLike)) # log likelihood function (returns scalar),n.chains,n.pars+1,byrow=T)
    #     log_post_like[iter,]=temp
    #   }

    if(iter>control_pars$burnin & iter%%control_pars$thin==0){
      thetaIdx=thetaIdx+1
    }

    print(paste0('iter ',iter,'/',control_pars$n_iter))
  }
  # cluster initialization
  if(!control_pars$parallelType=='none'){
    parallel::stopCluster(cl=cl_use)
  }
  return(list('samples'=theta,'log_post_like'=log_post_like,'control_pars'=control_pars))
}
