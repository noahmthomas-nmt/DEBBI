QLog=function(use_theta,use_lambda,control_pars){
  #returns log density Q(theta|lambda)
  out=0
  out=stats::dnorm(x=use_theta[1:(control_pars$n_pars_model)],
            mean=use_lambda[1:(control_pars$n_pars_model)],
            sd=exp(use_lambda[((control_pars$n_pars_model)+1):(control_pars$n_pars_post)]),
            log=T)
  out[(out==-Inf) | is.na(out)]=control_pars$neg_inf
  return(sum(out));
}


QSample=function(use_lambda,control_pars,S){
  #returns a S by n_pars_model matrix sampled from Q(theta|lambda)
  out=matrix(NA,S,control_pars$n_pars_model)
  if(control_pars$use_QMC==T){
    if(quasiSeq=="sobol")quantileMat=randtoolbox::sobol(S,control_pars$n_pars_model)
    if(quasiSeq=="halton")quantileMat=randtoolbox::halton(S,control_pars$n_pars_model)
  } else {
    quantileMat=matrix(stats::runif(S,0,1),S,control_pars$n_pars_model)
  }
  for(i in 1:n_pars_model){
    qs=quantileMat[,i]
    out[,i]=stats::qnorm(qs,use_lambda[i],exp(use_lambda[control_pars$n_pars_model+i]))
  }
  return(out)
}

KLHat=function(lambda,LogPostLike,control_pars,S,...){
  #monte carlo approximation KL divergence up to a constant

  out=0 #initalize output vector

  #sample from q
  theta_mat<-QSample(use_lambda=lambda,control_pars,S)

  #calc mean differences in log densities for theta_mat
  for(s in 1:S){ #can parallelize/vector this
    out=out+QLog(use_theta = theta_mat[s,],use_lambda = lambda,control_pars)/S
    out=out-LogPostLike(theta_mat[s,],...)/S
  }


  return(out)
}

AlgoParsDEVI_meanField=function(n_pars,
                        n_chains=NULL,
                        n_iter=1000,
                        init_sd=0.01,
                        init_center=0,
                        n_cores_use=1,
                        step_size=NULL,
                        jitter_size=1e-6,
                        parallel_type='none',
                        use_QMC=T,
                        qmc_type='sobol',
                        n_samples_ELBO=25,
                        n_samples_LRVB=25,
                        LRVB_correction=T,
                        neg_inf=-750){
  out=list('n_pars_models'=n_pars,
           'n_chains'=n_chains,
           'n_iter'=n_iter,
           'init_sd'=init_sd,
           'init_center'=init_center,
           'n_cores_use'=n_cores_use,
           'step_size'=step_size,
           'jitter_size'=jitter_size,
           'parallel_type'=parallel_type,
           'purify'=Inf)

  return(out)
}
