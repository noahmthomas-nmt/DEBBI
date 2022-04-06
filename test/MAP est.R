

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

AlgoParsDEMAP=function(n_pars,
                       n_chains=NULL,
                       n_iter=1000,
                       init_sd=0.01,
                       init_center=0,
                       n_cores_use=1,
                       step_size=NULL,
                       jitter_size=1e-6,
                       crossover_rate=1,
                       parallel_type='none',
                       return_trace=FALSE,
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

  # crossover_rate
  ### if null assign default value
  if(any(is.null(crossover_rate))){
    crossover_rate=1
  }
  ### catch errors
  crossover_rate=as.numeric(crossover_rate)
  if(any(!is.finite(crossover_rate))){
    stop('ERROR: crossover_rate is not finite')
  }
  else if(any(crossover_rate>1) | any(crossover_rate<=0) | length(crossover_rate)>1){
    stop('ERROR: crossover_rate must be a numeric scalar on the interval (0,1]')
  } else if(is.complex(crossover_rate)){
    stop('ERROR: crossover_rate cannot be complex')
  }

  #parallel_type
  validParType=c('none','FORK','PSOCK')
  ### assign NULL value default
  if(is.null(parallel_type)){
    parallel_type='none'
  }
  ### catch any errors
  if(!parallel_type %in% validParType){
    stop('ERROR: invalid parallel_type')
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
  n_iters_per_chain=floor((n_iter)/thin)
  ### catch errors
  if(n_iters_per_chain<1 | (!is.finite(n_iters_per_chain))){
    stop('ERROR: number of samples per chain is negative or non finite.
         n_iters_per_chain=floor((n_iter-burnin)/thin)')
  }

  # purify
  n_iters_per_chain=floor((n_iter)/thin)
  ### catch errors
  if(n_iters_per_chain<1 | (!is.finite(n_iters_per_chain))){
    stop('ERROR: number of samples per chain is negative or non finite.
         n_iters_per_chain=floor((n_iter-burnin)/thin)')
  }


  out=list('n_pars'=n_pars,
           'n_chains'=n_chains,
           'n_iter'=n_iter,
           'init_sd'=init_sd,
           'init_center'=init_center,
           'n_cores_use'=n_cores_use,
           'step_size'=step_size,
           'crossover_rate'=crossover_rate,
           'jitter_size'=jitter_size,
           'parallel_type'=parallel_type,
           'thin'=thin,
           'purify'=Inf,
           'n_iters_per_chain'=n_iters_per_chain,
           'return_trace'=return_trace)

  return(out)
}


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
  parent_chain_indices = sample(c(1:n_chains)[-chain_index],3,replace=F)

  # mate parents for proposal
  pars_use[par_indices] = current_pars[parent_chain_indices[3],par_indices] +
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

QSample=function(use_lambda,control_pars,S){
  # returns a collapsed vector for a S by n_pars_model matrix sampled from Q(theta|lambda)
  if(control_pars$use_QMC==T){
    if(control_pars$quasi_rand_seq=='sobol')quantileMat=c(t(randtoolbox::sobol(S,control_pars$n_pars_model)))
    if(control_pars$quasi_rand_seq=='halton')quantileMat=c(t(randtoolbox::halton(S,control_pars$n_pars_model)))
  } else {
    quantileMat=stats::runif(S*control_pars$n_pars_model,0,1)
  }

  out=stats::qnorm(quantileMat,rep(use_lambda[1:control_pars$n_pars_model],S),rep(exp(use_lambda[(control_pars$n_pars_model+1):(control_pars$n_pars_dist)]),S))
  return(out)
}


KLHat=function(lambda,LogPostLike,control_pars,S,...){
  # monte carlo approximation KL divergence up to a constant

  out=0 #initalize output vector

  # sample from q
  theta_mat <- QSample(use_lambda=lambda,control_pars,S)

  # calc mean differences in log densities for theta_mat
  q_log_density=QLog(theta_mat,use_lambda = lambda,control_pars,S)
  post_log_density=mean(apply(matrix(theta_mat,ncol=control_pars$n_pars_model,byrow = T),MARGIN = 1,FUN=LogPostLike,...))
  out <- q_log_density-post_log_density


  return(out)
}



ELBO=function(lambda,LogPostLike,control_pars,S,...){
  # monte carlo approximation ELBO
  out=KLHat(lambda,LogPostLike,control_pars,S,...)*-1
  return(out)
}

AlgoParsDEVI=function(n_pars,
                      par_names=NULL,
                      n_chains = NULL,
                      n_iter = 1000,
                      init_sd = 0.01,
                      init_center = 0,
                      n_cores_use = 1,
                      step_size = NULL,
                      jitter_size = 1e-6,
                      parallel_type = 'none',
                      use_QMC = T,
                      quasi_rand_seq = 'halton',
                      n_samples_ELBO = 10,
                      LRVB_correction = TRUE,
                      n_samples_LRVB = 25,
                      neg_inf = -750,
                      thin = 1,
                      burnin = 0,
                      return_trace = FALSE,
                      crossover_rate=1){
  # n_pars
  ### catch errors
  n_pars=as.integer(n_pars)
  if(any(!is.finite(n_pars))){
    stop('ERROR: n_pars is not finite')
  }  else if( n_pars<1 | length(n_pars)>1){
    stop('ERROR: n_pars must be a postitive integer scalar')
  }

  # par_names
  ### catch errors
  if(is.null(par_names)){
    par_names=paste0('par',1:n_pars)
  }  else if(!(length(par_names)==n_pars)){
    stop('ERROR: par_names does not match size of n_pars')
  }

  dist_par_names=c(paste0(par_names,'_MEAN'),paste0(par_names,'_VAR'))

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
  } else if(any(init_sd<0 | is.complex(init_sd))){
    stop('ERROR: init_sd must be positive and real-valued')
  } else if(!(length(init_sd)==1 | length(init_sd)==n_pars)){
    stop('ERROR: init_sd vector length must be 1 or n_pars')
  }
  if(any(init_sd==0)){
    warning('WARNING an init_sd value is 0')
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

  # crossover_rate
  ### if null assign default value
  if(any(is.null(crossover_rate))){
    crossover_rate=1
  }
  ### catch errors
  crossover_rate=as.numeric(crossover_rate)
  if(any(!is.finite(crossover_rate))){
    stop('ERROR: crossover_rate is not finite')
  }
  else if(any(crossover_rate>1) | any(crossover_rate<=0) | length(crossover_rate)>1){
    stop('ERROR: crossover_rate must be a numeric scalar on the interval (0,1]')
  } else if(is.complex(crossover_rate)){
    stop('ERROR: crossover_rate cannot be complex')
  }

  #parallel_type
  validParType=c('none','FORK','PSOCK')
  ### assign NULL value default
  if(is.null(parallel_type)){
    parallel_type='none'
  }
  ### catch any errors
  if(!parallel_type %in% validParType){
    stop('ERROR: invalid parallel_type')
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


  # use_QMC
  if(any(is.null(use_QMC))){
    use_QMC=TRUE
  }
  if(length(use_QMC)>1){
    stop('length(use_QMC)>1, please use a scalar logical')
  }
  if(!((use_QMC==0) | (use_QMC==1))){
    stop('ERROR: use_QMC must be a scalar logical')
  }
  use_QMC=as.logical(use_QMC)


  # LRVB_correction
  ### assign value if null
  if(any(is.null(LRVB_correction))){
    LRVB_correction=TRUE
  }
  ### catch errors
  if(length(LRVB_correction)>1){
    stop('length(LRVB_correction)>1, please use a scalar logical')
  }
  if(!((LRVB_correction==0) | (LRVB_correction==1))){
    stop('ERROR: LRVB_correction must be a scalar logical')
  }
  LRVB_correction=as.logical(LRVB_correction)

  if(LRVB_correction){
    # if using LRVB correction check for valid samples count
    # n_samples_LRVB
    ### assign value if null
    n_samples_LRVB=as.integer(n_samples_LRVB)
    if(any(is.null(n_samples_LRVB))){
      n_samples_LRVB=25
    }
    ### catch errors
    if(any(!is.finite(n_samples_LRVB))){
      stop('ERROR: n_samples_LRVB is not finite')
    }
    else if(any(n_samples_LRVB<1) | length(n_samples_LRVB)>1){
      stop('ERROR: n_samples_LRVB must be a scalar postive integer')
    }
  }

  # n_samples_ELBO
  ### assign value if null
  n_samples_ELBO=as.integer(n_samples_ELBO)
  if(any(is.null(n_samples_ELBO))){
    n_samples_ELBO=10
  }
  ### catch errors
  if(any(!is.finite(n_samples_ELBO))){
    stop('ERROR: n_samples_ELBO is not finite')
  }
  else if(any(n_samples_ELBO<1) | length(n_samples_ELBO)>1){
    stop('ERROR: n_samples_ELBO must be a scalar postive integer')
  }

  # quasi_rand_seq
  quasi_rand_seq=tolower(as.character(quasi_rand_seq))
  valid_quasi_rand_seqs=c("sobol","halton")
  ### assign NULL value default
  if(is.null(quasi_rand_seq)){
    quasi_rand_seq='sobol'
  }
  ### catch any errors
  if(!quasi_rand_seq %in% valid_quasi_rand_seqs){
    stop(paste0('ERROR: invalid quasi_rand_seq; must be one of: ', paste(valid_quasi_rand_seqs,sep = ',')))
  }

  # neg_inf
  ### assign NULL value default
  if(any(is.null(neg_inf))){
    neg_inf=-750
  }
  ### catch any errors
  if(length(neg_inf)>1){
    stop('length(neg_inf)>1, please use a scalar numeric')
  }
  if(!is.numeric(neg_inf)){
    stop('neg_inf must be a numeric')
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

  #nSamples Per Chains
  n_iters_per_chain=floor((n_iter-burnin)/thin)
  ### catch errors
  if(n_iters_per_chain<1 | (!is.finite(n_iters_per_chain))){
    stop('ERROR: number of iters per chain is negative or non finite.
           n_iters_per_chain=floor((n_iter-burnin)/thin)')
  }

  out=list('n_pars_model'=n_pars,
           'par_names'=par_names,
           'n_chains'=n_chains,
           'n_iter'=n_iter,
           'init_sd'=init_sd,
           'init_center'=init_center,
           'n_cores_use'=n_cores_use,
           'step_size'=step_size,
           'jitter_size'=jitter_size,
           'crossover_rate'=crossover_rate,
           'parallel_type'=parallel_type,
           'return_trace'=return_trace,
           'purify'=Inf,
           'use_QMC'=use_QMC,
           'quasi_rand_seq'=quasi_rand_seq,
           'n_samples_ELBO'=n_samples_ELBO,
           'n_samples_LRVB'=n_samples_LRVB,
           'LRVB_correction'=LRVB_correction,
           'thin'=thin,
           'neg_inf'=neg_inf,
           'n_pars_dist'=2*n_pars,
           'n_iters_per_chain'=n_iters_per_chain)

  return(out)
}

dataExample=matrix(stats::rnorm(200,c(-1,1),c(.1,.1)),nrow=50,ncol=2,byrow = TRUE)
##list parameter names
par_names_example=c("mu_1","mu_2")

#log posterior likelihood function = log likelihood + log prior | returns a scalar
LogPostLikeExample=function(x,data,par_names){
  out=0

  names(x)<-par_names

  # log prior
  out=out+sum(dnorm(x["mu_1"],0,sd=1,log=TRUE))
  out=out+sum(dnorm(x["mu_2"],0,sd=1,log=TRUE))

  # log likelihoods
  out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=TRUE))
  out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=TRUE))

  return(out)
}

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

out=DEVI(LogPostLike=LogPostLikeExample,
         control_pars=AlgoParsDEVI(n_pars=length(par_names_example),
                                   n_iter=200,
                                   n_samples_ELBO = 5,
                                   n_chains=25,return_trace=T,crossover_rate = .8,use_QMC = T),
         data=dataExample,
         par_names = par_names_example)
par(mfrow=c(2,2))
matplot(out$lambda_trace[,,2],type="l")
matplot(out$lambda_trace[,,1],type="l")
matplot(out$lambda_trace[,,3],type="l")
matplot(out$lambda_trace[,,4],type="l")
