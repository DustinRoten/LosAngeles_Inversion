compose.bayesian.equation.v2 <- function(homedir = NULL, path = NULL,
                                         lambda = NULL, V_shat.path = NULL,
                                         xco2.error.path = 'ext/included.errors.csv',
                                         which.sam = 1, obs.error = 0) {
  
  #source all necessary functions
  setwd(homedir); source('r/dependencies.r')
  
  if(!dir.exists(file.path(path, '../inversion_results')))
    dir.create(file.path(path, '../inversion_results'))
    
  ### Make S (Error) ###
  #read and average prior
  prior.list <- list.files(file.path(path, 'include/sectors'),
                           full.names = TRUE)
  for(i in 1:length(prior.list)) {
    
    #generate the new sector
    sector <- lambda[[i]]*brick(prior.list[[i]])
    
    if(i == 1)
      prior <- sector
    if(i > 1)
      prior <- prior + sector
    
  }; avg.prior <- calc(prior, mean)
  
  #read and average truthf
  truth <- list.files(file.path(path, 'include'),
                      pattern = 'truth', full.names = TRUE)
  avg.truth <- calc(brick(truth), mean)
  
  ratio <- avg.truth/avg.prior
  ratio[is.na(ratio)] <- 0
  
  spatial.error <- abs(ratio - lambda)
  spatial.error[is.na(spatial.error)] <- 0
  names(spatial.error) <- names(lambda)
  
  #construct the spatial covariance matrix
  #modified remove the top 5% of values.
  if(length(V_shat.path) == 0) {
    
    #generate S and/or V_shat
    S <- spatial.correlation(spatial.error = spatial.error,
                             plot.title = paste0('SAM ', which.sam),
                             plot.output = file.path(path, '..',
                                                     'inversion_results',
                                                     paste0('var_flux_',
                                                            which.sam)))
    
    #save the file
    saveRDS(S, file = file.path(path, '..', 'inversion_results',
                                paste0('Vshat_', 0)))
    
  } else {S <- readRDS(V_shat.path)}
  
  #obtain the M and Z matrices
  #combine observations
  Z_tot <- readRDS(file = list.files(file.path(path, 'include'),
                                     pattern = paste0('Z_', which.sam),
                                     full.names = TRUE))
  
  #combine the M matrices
  M_tot <- readRDS(file = list.files(file.path(path, 'include'),
                                     pattern = paste0('M_', which.sam),
                                     full.names = TRUE))
  
  ### Make R ###
  xco2.error <- data.frame(error = abs(Z_tot$Z - rowSums(M_tot)),
                           lon = Z_tot$LON,
                           lat = Z_tot$LAT)
  
  R <- spatial.correlation(spatial.error = xco2.error,
                           vgm.binwidth = 4,
                           vgm.cutoff = 30,
                           included.xco2.errors = xco2.error.path,
                           plot.title = paste0(which.sam, ' SAM'),
                           plot.output = file.path(path, '..',
                                                   'inversion_results',
                                                   paste0('var_xco2_',
                                                          which.sam,
                                                          '_.jpg')))
  
  #make a vector of lambda values
  for(i in 1:nlayers(lambda)) {
    if(i == 1)
      lambda_vec <- values(lambda[[i]])
    if(i > 1)
      lambda_vec <- c(lambda_vec, values(lambda[[i]]))
  }

  #Generate the new lambda vector
  #fancy math!
  term_1 <- S %*% t(M_tot)
  term_2 <- solve((M_tot %*% term_1) + R)
  term_3 <- Z_tot[,1] - (M_tot %*% lambda_vec)
  new.lambda <- lambda_vec + (term_1 %*% (term_2 %*% term_3))
  
  #create a new raster for the cellwise corrections
  lambda_raster <- lambda; values(lambda_raster) <- new.lambda[,1]
  
  #save the raster
  writeRaster(lambda_raster, format = 'CDF',
              filename = file.path(path, '../inversion_results',
                                   paste0('lambda_', which.sam,
                                          '_.nc')),
              overwrite = TRUE)
  
  remove('term_3'); remove('lambda')
  remove('lambda_raster'); remove('lambda_vec')
  gc()
  
  #Generate the new covariance matrix
  #V_shat = S - SM^T (MSM^T + R)^-1 MS
  #V_shat = (M^t R^-1 M + S^-1)^-1
  #even more fancy math!
  term_3 <- M_tot %*% S
  V_shat <- S - (term_1 %*% (term_2 %*% term_3))
  
  #save the new covariance matrix
  saveRDS(V_shat, file = file.path(path, '../inversion_results',
                                   paste0('Vshat_', which.sam)))
  
  #clean up
  remove('term_1'); remove('term_2')
  gc()
  
}
