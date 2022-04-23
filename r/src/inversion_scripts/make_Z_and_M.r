make_Z_and_M <- function(homedir = NULL, path = NULL, which.sam = 1, obs.error = 0) {
  
  #source all necessary functions
  setwd(homedir); source('r/dependencies.r')
  
  #initialize the output vector and matrix
  Z <- LON <- LAT <- NULL
  M <- NULL
  
  #create the output folder for the inversion
  if(!dir.exists(file.path(path, 'inversion')))
    dir.create(file.path(path, 'inversion'))
  
  #obtain the true emissions
  truth <- list.files(file.path(path, 'include'),
                      pattern = 'truth',
                      full.names = TRUE)
  truth <- brick(truth)
  
  #list the prior emissions inventory sectors
  prior.list <- list.files(file.path(path, 'include',
                                     'sectors'),
                           pattern = 'sector.prior',
                           full.names = TRUE)
  
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
      m_s <- values(calc(prior.sector*foot, sum))
      m <- c(m, m_s)
    } #close the prior loop
    
    #add the generated line to the M matrix
    M <- rbind(M, t(m))
    
  } #close the footprint loop (Z and M)
  
  #save the observation vector
  Z <- data.frame(Z, LON, LAT)
  saveRDS(Z, file = file.path(path, 'inversion',
                              paste0('Z_', which.sam, '_.RDS')))
  
  #save the observation vector
  saveRDS(M, file = file.path(path, 'inversion',
                              paste0('M_', which.sam, '_.RDS')))
  
} #closes the function