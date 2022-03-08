

Crossover.MC=function(chain_index, # which chain you are updating
                      par_indices, # which parameters you are updating (int vector)
                      current_theta,  # current parameter values for chain (numeric vector)
                      current_log_post_dens, # corresponding log post dens for (numeric vector)
                      Log.Post.Dens, # log likelihood function (returns scalar)
                      step_size=.8,
                      jitter_size=1e-6,
                      n_chains, ... ){
  
  # get statistics about chain
  like_use = current_log_post_dens[chain_index]
  theta_use = current_theta[chain_index,]		
  
  # sample parent chains
  parent_chain_indices = sample(c(1:n_chains)[-chain_index],2,replace=F)
  
  # mate parents for proposal
  theta_use[par_indices] = theta_use[par_indices] +
    step_size*(current_theta[parent_chain_indices[1],par_indices] -
                 current_theta[parent_chain_indices[2],par_indices]) +
    runif(1,-jitter_size,jitter_size)
  
  theta_use = matrix(theta_use,1,length(theta_use))
  
  # get log like
  like_proposal = Log.Post.Dens(theta_use,...)
  if(is.na(like_proposal))like_proposal = -Inf
  
  # metropolis hasting acceptance rule
  if(runif(1) < exp(like_proposal - like_use)) {							
    current_theta[chain_index,] <- theta_use
    current_log_post_dens[chain_index] <- like_proposal
  }
  
  return(c(current_log_post_dens[chain_index],current_theta[chain_index,]))
  
}

# Recalculate the log posterior density for each sample. 
# Helps when likelihood function is probabilistic
Purify.MC=function(chain_index,current_theta,Log.Post.Dens,...){
  
  theta_use = current_theta[chain_index,]		
  theta_use = matrix(theta_use,1,length(theta_use))
  
  like=Log.Post.Dens(theta_use,...)
  
  if(!is.na(like) & like!=-Inf){
    current_log_post_dens[chain_index]=like    
  }
  
  return(c(current_log_post_dens[chain_index],current_theta[chain_index,]))
}

DEMCMC.Algo.Pars=function(n_pars, 
                          n_chains=NULL, 
                          n_iter=1000, 
                          init_sd=0.01, 
                          init_center=0,
                          n_cores_use=1, 
                          step_size=NULL,
                          jitter_size=1e-6,
                          parallelType='none',
                          burnin=0,
                          thin=1){
  # n_pars
  ### catch errors
  n_pars=as.integer(n_pars)
  if(any(!is.finite(n_pars))){
    stop('ERROR: n_pars is not finite')
  }  else if( n_pars<1 | length(n_pars)>1){
    stop('ERROR: n_pars must be a postitive integer scalar')
  }
  
  # n_chains  
  ### if null assign default value
  if(is.null(n_chains)){
    n_chains=max(3*n_pars,4)
  }
  ### catch errors
  n_chains=as.integer(n_chains)
  if(any(!is.finite(n_chains))){
    stop('ERROR: n_chains is not finite')
  } 
  else if( n_chains<4 | length(n_chains)>1){
    stop('ERROR: n_chains must be a postitive integer scalar, and atleast 4')
  }
  
  # n_iter
  ### if null assign default value
  if(is.null(n_iter)){
    n_iter=1000
  }
  ### catch errors
  n_iter=as.integer(n_iter)
  if(any(!is.finite(n_iter))){
    stop('ERROR: n_iter is not finite')
  } 
  else if( n_iter<4 | length(n_iter)>1){
    stop('ERROR: n_iter must be a postitive integer scalar, and atleast 4')
  }
  
  # init_sd
  init_sd=as.numeric(init_sd)
  if(any(!is.finite(init_sd))){
    stop('ERROR: init_sd is not finite')
  } else if(any(init_sd<=0 | is.complex(init_sd))){
    stop('ERROR: init_sd must be positive and real-valued')
  } else if(!(length(init_sd)==1 | length(init_sd)==n_pars)){
    stop('ERROR: init_sd vector length must be 1 or n_pars')
  }
  
  # init_center
  init_center=as.numeric(init_center)
  if(any(!is.finite(init_center))){
    stop('ERROR: init_center is not finite')
  } else if(any(is.complex(init_center))){
    stop('ERROR: init_center must be real valued')
  } else if(!(length(init_center)==1 | length(init_center)==n_pars)){
    stop('ERROR: init_center vector length must be 1 or n_pars')
  }
  
  # n_cores_use
  ### assign NULL value default
  if(is.null(n_cores_use)){
    n_cores_use=1  
  }
  ### catch any errors
  n_cores_use=as.integer(n_cores_use)
  if(any(!is.finite(n_cores_use))){
    stop('ERROR: n_cores_use is not finite')
  } else if( n_cores_use<1 | length(n_cores_use)>1){
    stop('ERROR: n_cores_use must be a postitive integer scalar, and atleast 1')
  }
  
  
  # step_size
  ### assign NULL value default
  if(is.null(step_size)){
    step_size=2.38/sqrt(2*n_pars) # step size recommend in ter braak's 2006 paper
  }
  ### catch any errors
  if(any(!is.finite(step_size))){
    stop('ERROR: step_size is not finite')
  } else if(any(step_size<=0 | is.complex(step_size))){
    stop('ERROR: step_size must be positive and real-valued')
  } else if(!(length(step_size)==1)){
    stop('ERROR: step_size vector length must be 1 ')
  }
  
  #jitter_size
  ### assign NULL value default
  if(is.null(jitter_size)){
    jitter_size=1e-6
  }
  ### catch any errors
  if(any(!is.finite(jitter_size))){
    stop('ERROR: jitter_size is not finite')
  } else if(any(jitter_size<=0 | is.complex(jitter_size))){
    stop('ERROR: jitter_size must be positive and real-valued')
  } else if(!(length(jitter_size)==1)){
    stop('ERROR: jitter_size vector length must be 1 ')
  }
  
  #parallelType
  validParType=c('none','FORK','PSOCK')
  ### assign NULL value default
  if(is.null(parallelType)){
    parallelType='none'
  }
  ### catch any errors
  if(!parallelType %in% validParType){
    stop('ERROR: invalid parallelType')
  } 
  
  # burnin
  ### if null assign default value
  if(is.null(burnin)){
    burnin=0
  }
  ### catch errors
  burnin=as.integer(burnin)
  if(any(!is.finite(burnin))){
    stop('ERROR: burnin is not finite')
  } 
  else if(any(burnin<0) | any(burnin>=n_iter) | length(burnin)>1){
    stop('ERROR: burnin must be a scalar integer from the interval [0,n_iter)')
  }
  
  # thin
  ### if null assign default value
  if(is.null(thin)){
    thin=0
  }
  ### catch errors
  thin=as.integer(thin)
  if(any(!is.finite(thin))){
    stop('ERROR: thin is not finite')
  } 
  else if(any(thin<1) | length(thin)>1){
    stop('ERROR: thin must be a scalar postive integer')
  }
  
  #nSamples Per Chains
  nSamplesPerChain=floor((n_iter-burnin)/thin)
  ### catch errors
  if(nSamplesPerChain<1 | (!is.finite(nSamplesPerChain))){
    stop('ERROR: number of samples per chain is negative or non finite. 
         nSamplesPerChain=floor((n_iter-burnin)/thin)')
  }
  
  # purify
  nSamplesPerChain=floor((n_iter-burnin)/thin)
  ### catch errors
  if(nSamplesPerChain<1 | (!is.finite(nSamplesPerChain))){
    stop('ERROR: number of samples per chain is negative or non finite. 
         nSamplesPerChain=floor((n_iter-burnin)/thin)')
  }
  
  
  out=list('n_pars'=n_pars,
           'n_chains'=n_chains, 
           'n_iter'=n_iter, 
           'init_sd'=init_sd, 
           'init_center'=init_center, 
           'n_cores_use'=n_cores_use,
           'step_size'=step_size,
           'jitter_size'=jitter_size,
           'parallelType'=parallelType,
           'burnin'=burnin,
           'thin'=thin,
           'purify'=Inf,
           'nSamplesPerChain'=nSamplesPerChain)
  
  return(out)
}
########################################## 
Run.DEMCMC=function(Log.Post.Dens,control_pars=DEMCMC.Algo.Pars(),...){
  
  
  # import values we will reuse throughout process
  # create memory structures for storing posterior
  nSamplesKeep=
  theta=array(NA,dim=c(control_pars$nSamplesPerChain,control_pars$n_chains,control_pars$n_pars))
  log_post_dens=matrix(-Inf,nrow=control_pars$nSamplesPerChain,ncol=control_pars$n_chains)
  
  
  # chain initialization
  print('initalizing chains...')
  for(chain_idx in 1:control_pars$n_chains){
    count=0
    while (log_post_dens[1,chain_idx]  == -Inf) {
      theta[1,chain_idx,] <- rnorm(control_pars$n_pars,control_pars$init_center,control_pars$init_sd)
      
      log_post_dens[1,chain_idx] <- Log.Post.Dens(theta[1,chain_idx,],...)  
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
    library(doParallel)
    
    print(paste0("initalizing ",
                 control_pars$parallelType," cluser with ",
                 control_pars$n_cores_use," cores"))
    
    registerDoParallel(control_pars$n_cores_use)
    cl_use <- makeCluster(control_pars$n_cores_use,
                          type=control_pars$parallelType)  
  }
  
  print("running DEMCMC")
  thetaIdx=1
  for(iter in 1:control_pars$n_iter){

    if(control_pars$parallelType=='none'){
      # cross over step
      temp=matrix(unlist(lapply(1:control_pars$n_chains,Crossover.MC,
                                par_indices=1:control_pars$n_pars,
                                current_theta=theta[thetaIdx,,],  # current parameter values for chain (numeric vector)
                                current_log_post_dens=log_post_dens[thetaIdx,], # corresponding log like for (numeric vector)
                                Log.Post.Dens=Log.Post.Dens, # log likelihood function (returns scalar)
                                step_size=control_pars$step_size,
                                jitter_size=control_pars$jitter_size,
                                n_chains=control_pars$n_chains,...)),control_pars$n_chains,control_pars$n_pars+1,byrow=T)
    } else {
      temp=matrix(unlist(parLapply(cl_use,1:control_pars$n_chains,Crossover.MC,
                                   par_indices=1:control_pars$n_pars,
                                   current_theta=theta[thetaIdx,,],  # current parameter values for chain (numeric vector)
                                   current_log_post_dens=log_post_dens[thetaIdx,], # corresponding log like for (numeric vector)
                                   Log.Post.Dens=Log.Post.Dens, # log likelihood function (returns scalar)
                                   step_size=control_pars$step_size,
                                   jitter_size=control_pars$jitter_size,
                                   n_chains=control_pars$n_chains,...)),control_pars$n_chains,control_pars$n_pars+1,byrow=T)
      
    }
    log_post_dens[thetaIdx,]=temp[,1]
    theta[thetaIdx,,]=temp[,2:(control_pars$n_pars+1)]
    
    #   # purification step
    #   if(iter%%purify.rate==0){
    #     temp=unlist(lapply(1:n_chains,Purify.MC,
    #                                 current_theta=theta[iter,],  # current parameter values for chain (numeric vector)
    #                                 current_log_post_dens=log_post_dens[iter,], # corresponding log like for (numeric vector)
    #                                 Log.Post.Dens=Log.Post.Dens)) # log likelihood function (returns scalar),n.chains,n.pars+1,byrow=T)
    #     log_post_dens[iter,]=temp
    #   }
    
    if(iter>control_pars$burnin & iter%%control_pars$thin==0){
      thetaIdx=thetaIdx+1
    }
    
    print(paste0('iter ',iter,'/',control_pars$n_iter))
  }
  # cluster initialization
  if(!control_pars$parallelType=='none'){
    stopCluster(cl=cl_use)
  }
  return(list('samples'=theta,'log_post_density'=log_post_dens,'control_pars'=control_pars))
}


