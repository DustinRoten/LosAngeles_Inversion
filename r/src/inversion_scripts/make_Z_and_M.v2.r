make_Z_and_M.v2 <- function(homedir = NULL, path = NULL, which.sam = 1,
                            lambda = NULL, obs.error = 0) {
  
  #source all necessary functions
  setwd(homedir); source('r/dependencies.r')
  
  #initialize the output vector and matrix
  Z <- LON <- LAT <- NULL
  M <- NULL
  
  #obtain the true emissions
  truth <- list.files(file.path(path, 'include'),
                      pattern = 'truth',
                      full.names = TRUE)
  truth <- brick(truth)
  
  #list the prior emissions inventory sectors
  prior.list <- list.files(file.path(path, 'include',
                                     'sectors'),
                           full.names = TRUE)
  
  #check on the status of lambda. If no lambda exists, make the default.
  if(length(lambda) == 0) {
    for(i in 1:length(prior.list)) {
      tmp <- calc(brick(prior.list[i]), mean)
      tmp[abs(tmp) != 0] <- 1
      
      if(i == 1)
        lambda <- tmp
      if(i > 1)
        lambda <- addLayer(lambda, tmp)
    }
    names(lambda) <- gsub('.nc', '', basename(prior.list))
    writeRaster(lambda,
                filename = file.path(path, '../inversion_results',
                                     paste0('lambda_', 0, '_')),
                format = 'CDF', overwrite = TRUE)

  } else {
    lambda <- brick(lambda)
    names(lambda) <- gsub('.nc', '', basename(prior.list))
  }
  
  #list the footprints
  foot.list <- list.files(file.path(path, 'footprints'),
                          full.names = TRUE)
  foot.list <- list.files(foot.list, full.names = TRUE)
  
  #begin the process of building Z and M
  for(i in 1:length(foot.list)) {
    
    #read in the footprint file
    foot <- brick(foot.list[i])
    
    #calculate the "true" observation and add instrument error
    Z[i] <- sum(cellStats(truth*foot, sum)) + rnorm(1, 0, obs.error)
    LON[i] <- as.numeric(gsub('.nc', '',
                              str_split_fixed(basename(foot.list[i]),
                                              pattern = '_', n = 3)[2]))
    LAT[i] <- as.numeric(gsub('.nc', '',
                              str_split_fixed(basename(foot.list[i]),
                                              pattern = '_', n = 3)[3]))
    
    #create a vector of ppm influence values that includes all sectors
    m <- NULL
    for(j in 1:length(prior.list)) {
      prior.sector <- brick(prior.list[j])
      m_s <- values(calc(lambda[[j]]*prior.sector*foot, sum))
      m <- c(m, m_s)
    } #close the prior loop
    
    #add the generated line to the M matrix
    M <- rbind(M, t(m))
    
  } #close the footprint loop (Z and M)
  
  #save the observation vector
  Z <- data.frame(Z, LON, LAT)
  saveRDS(Z, file = file.path(path, 'include',
                              paste0('Z_', which.sam, '_.RDS')))
  
  #save the observation vector
  saveRDS(M, file = file.path(path, 'include',
                              paste0('M_', which.sam, '_.RDS')))
  
} #closes the function
