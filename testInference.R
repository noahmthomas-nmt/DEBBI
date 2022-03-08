
setwd("~/Dropbox/DEBBI")

data=matrix(rnorm(100,c(-1,1),c(1,1)),nrow=50,ncol=2,byrow = T)

par_names_1=c("mu_1","mu_2")

Log.Post.Dens.1=function(x,data,par_names){
  out=0
  
  names(x)<-par_names
  
  # log prior
  out=out+sum(dnorm(x["mu_1"],0,sd=1,log=T))
  out=out+sum(dnorm(x["mu_2"],0,sd=1,log=T))
  
  # log likelihoods
  out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=T))
  out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=T))
  
  return(out)
}

source("DEMCMC.R")

out <- Run.DEMCMC(Log.Post.Dens=Log.Post.Dens.1,
                  control_pars=DEMCMC.Algo.Pars(n_pars=length(par_names_1),
                                     n_iter=1000, n_chains=12,
                                     init_sd=.01,init_center=0,
                                     n_cores_use=1,step_size=NULL,
                                     jitter_size=1e-6,parallelType = "none"),data=data,
                  par_names = par_names_1)


matplot(out$samples[,,],type='l')

matplot(out$log_post_density,type='l')
