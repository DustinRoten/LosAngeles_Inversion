compose.bayesian.equation <- function(homedir = NULL, path = NULL,
                                      xco2.error.path = 'ext/included.errors.csv',
                                      which.sams = 1, obs.error = 0) {
  
  #source all necessary functions
  setwd(homedir); source('r/dependencies.r')
  
  if(!dir.exists(file.path(path, 'inversion_results')))
    dir.create(file.path(path, 'inversion_results'))
  
  #list the SAMs needed
  SAM.list <-
    list.files(path, full.names = TRUE, pattern = 'out_')[1:which.sams]
  
  #sector list
  sectors <- unique(list.files(file.path(SAM.list, 'include/sectors')))
  
  ### Make Lambda ###
  lambda_0 <- make_lambda(SAM.list)
  lambda_0.vec <- as.vector(values(lambda_0))
  lambda_0.vec[is.na(lambda_0.vec)] <- 0
  lambda_0.vec[is.nan(lambda_0.vec)] <- 0
  
  # Identify the locations of sector-specific flux
  bin.lambda <- lambda_0
  bin.lambda[bin.lambda > 0] <- 1
  bin.lambda[bin.lambda < 0] <- 0 #this shouldn't be needed
  
  ### Make S (Error) ###
  S <- spatial.correlation(spatial.error =
                             abs((1.25*lambda_0)-lambda_0)*bin.lambda,
                           plot.title = paste0(which.sams, ' SAMs'),
                           plot.output = file.path(path,
                                                   'inversion_results',
                                                   paste0('var_flux_',
                                                          which.sams)))
  
  #obtain the M and Z matrices
  Z_tot <- M_tot <- NULL
  for(i in 1:length(SAM.list)) {
    
    #combine observations
    Z_i <- readRDS(file = list.files(file.path(SAM.list[i], 'inversion'),
                                     pattern = paste0('Z_', i),
                                     full.names = TRUE))
    Z_tot <- rbind(Z_tot, Z_i)
    
    #combine the M matrices
    M_i <- readRDS(file = list.files(file.path(SAM.list[i], 'inversion'),
                                     pattern = paste0('M_', i),
                                     full.names = TRUE))
    M_tot <- rbind(M_tot, M_i)
    
  }
  
  ### Make R ###
  xco2.error <- data.frame(error = abs(Z_tot$Z - rowSums(M_tot)),
                           lon = Z_tot$LON,
                           lat = Z_tot$LAT)
  
  R <- spatial.correlation(spatial.error = xco2.error,
                           vgm.binwidth = 2,
                           vgm.cutoff = 25,
                           included.xco2.errors = xco2.error.path,
                           plot.title = paste0(which.sams, ' SAMs'),
                           plot.output = file.path(path,
                                                   'inversion_results',
                                                   paste0('var_xco2_',
                                                          which.sams,
                                                          '_.jpg')))
  
  #fancy math!
  term_1 <- S %*% t(M_tot)
  term_2 <- solve((M_tot %*% term_1) + R)
  term_3 <- Z_tot[,1] - (M_tot %*% lambda_0.vec)
  lambda = lambda_0.vec + (term_1 %*% (term_2 %*% term_3))
  
  saveRDS(lambda, file = paste0('TMP_', which.sams, '.RDS'))
  
  #create a new raster for the cellwise corrections
  lambda_raster <- lambda_0; values(lambda_raster) <- lambda[,1]
  
  #save the raster
  writeRaster(lambda_raster, format = 'CDF',
              filename = file.path(path, 'inversion_results',
                                   paste0('lambda_', which.sams, '_.nc')),
              overwrite = TRUE)
  
}