# Recalculate the log posterior likelihood for each sample.
# Helps when model likelihood function is probabilistic
PurifyMC=function(chain_index,current_theta,current_log_post_like,LogPostLike,...){

  theta_use = current_theta[chain_index,]
  theta_use = matrix(theta_use,1,length(theta_use))

  like=LogPostLike(theta_use,...)

  if(!is.na(like) & like!=-Inf){
    current_log_post_like[chain_index]=like
  }

  return(c(current_log_post_like[chain_index],current_theta[chain_index,]))
}
