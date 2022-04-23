iterate_lambda <- function(path = NULL) {
  
  #####################
  ### Make lambda_0 ###
  #####################
  #list the prior emissions inventory sectors
  prior.list <- list.files(file.path(path, 'include', 'sectors'),
                           full.names = TRUE)
  
  truth.list <- list.files(file.path(path, 'include'),
                           pattern = 'truth',
                           full.names = TRUE)
  truth <- brick(truth.list)
  
  scaling_factor <- list.files(file.path(path, '..', 'inversion_results'),
                               pattern = 'scaling_factor',
                               full.names = TRUE)
  scaling_factor <- scaling_factor[length(scaling_factor)]
  
  if(length(scaling_factor) == 0) {
    #determine the total prior emissions
    for(j in 1:length(prior.list)) {
      if(j == 1)
        total.prior <- brick(prior.list[j])
      if(j > 1)
        total.prior <- total.prior + brick(prior.list[j])
    }
  } else {
    scaling_factor <- brick(scaling_factor)
    #determine the total prior emissions
    for(j in 1:length(prior.list)) {
      if(j == 1)
        total.prior <- scaling_factor[j]*brick(prior.list[j])
      if(j > 1)
        total.prior <-
          total.prior + scaling_factor[j]*brick(prior.list[j])
    }
  }
  
  #construct lambda_0
  for(j in 1:length(prior.list)) {
    prior.sector <- brick(prior.list[j])
    if(j == 1) {
      lambda_tmp <- calc(prior.sector/total.prior, mean)
    }
    
    if(j > 1) {
      proportion <- calc(prior.sector/total.prior, mean)
      lambda_tmp <- addLayer(lambda_tmp, proportion)
    }
  }
  
  if(i == 1)
    lambda_0 <- lambda_tmp
  if(i > 1)
    lambda_0 <- lambda_0 + lambda_tmp
  
  names(lambda_0) <- basename(prior.list)
  
  return(lambda_0)
  #####################
  #####################
}