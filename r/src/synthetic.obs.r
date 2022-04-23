synthetic.obs <- function(footprints = NULL, start.time = NULL, errors = NULL, obs.error = NA,
                          inner.emissions = NULL, modeled.emissions = NULL, outer.emissions = NULL,
                          biosphere.emissions = NULL, obs.output = NULL, background.output = NULL,
                          std.bkg = T, anthropogenic.emissions.only = F, filter.noise = F,
                          regridded = FALSE) {
  
  #' If modeled and true emissions need to be regridded,
  #' ensure that the original observations are obtained
  #' BEFORE the regridding happens!
  if(regridded)
    original.obs <- read.csv(file.path(dirname(obs.output),
                                       'xco2_observations.csv'))
  
  #the variable has a long name (for description)
  #it is shortened here
  anthro <- anthropogenic.emissions.only
  
  #read in the inner and outer emissions files
  inner.emissions.raster <- brick(inner.emissions)
  modeled.emissions.raster <- brick(modeled.emissions)
  
  #only read if anthro is FALSE!
  if(!anthro) {
    
    #read in outer domain (background) emissions
    outer.emissions.raster <- brick(outer.emissions)
    
    #bio emissions are optional
    if(!is.na(biosphere.emissions))
      bio.emissions.raster <- brick(biosphere.emissions)
  }
  
  site_code.list <- NULL
  XCO2_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 7))
  names(XCO2_df) <- c('date.time', 'lon', 'lat',
                      'modeled.xco2', 'true.xco2', 'xco2',
                      'footprint_path')
  bkg.XCO2_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 2))
  names(bkg.XCO2_df) <- c('time', 'avg_CO2')
  for(i in 1:length(footprints)) {
    
    #get the footprint's directory name
    site_code.list[i] <- basename(dirname(footprints[i]))
    
    #begin with a progress message
    cat(paste0('Generating synthetic observations: ',
               round(100*i/length(footprints), 2),
               '% complete.     '), '\r')
    
    #get date.time here
    tmp.date.time <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[1]
    
    #get lon here
    tmp.lon <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[2]
    tmp.lon <- gsub('.nc', '', tmp.lon)
    
    #get lat here
    tmp.lat <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[3]
    tmp.lat <- gsub('.nc', '', tmp.lat)
    
    tmp.footprint <- brick(footprints[i])
    tmp.layers <- names(tmp.footprint)
    
    hourly.XCO2 <- hourly.XCO2.modeled <- NULL
    for(j in 1:length(tmp.layers)) {
      
      #read in each footprint layer
      eval(parse(text = paste0('footprint.layer <- tmp.footprint$',
                               tmp.layers[j])))
      
      #read in each inner (truth) layer
      eval(parse(text = paste0('inner.layer <- inner.emissions.raster$',
                               tmp.layers[j])))
      
      #read in each modeled (prior) layer
      eval(parse(text = paste0('modeled.layer <- modeled.emissions.raster$',
                               tmp.layers[j])))
      
      #Check whether the bio and background emissions are required
      if(!anthro) {
        
        #obtain the outer layer
        #read in each modeled (prior) layer
        eval(parse(text = paste0('outer.layer <- outer.emissions.raster$',
                                 tmp.layers[j])))
        
        #same as above with bio layers
        if(!is.na(biosphere.emissions)) {
          
          #read in each modeled (prior) layer
          eval(parse(text = paste0('bio.layer <- bio.emissions.raster$',
                                   tmp.layers[j])))
          
        } #closes bio emissions
      } #closes anthro only
      
      if(anthro) {
        
        #inner domain emissions
        hourly.XCO2[j] <-
          suppressWarnings(cellStats(footprint.layer*inner.layer,
                                     sum))
        hourly.XCO2.modeled[j] <-
          suppressWarnings(cellStats(footprint.layer*modeled.layer,
                                     sum))
          
      } else if(!anthro) {
        
        #determine the domain's total emissions
        total.layer <-
          extend(inner.layer, outer.layer, value = 0) +
          outer.layer
        
        #same for the modeled domains total emissions
        total.layer.mod <-
          extend(modeled.layer, outer.layer, value = 0) +
          outer.layer
        
        #add the (optional) biospheric layer
        if(!is.na(biosphere.emissions)) {
          total.layer <- total.layer + bio.layer
          total.layer.mod <- total.layer.mod + bio.layer
        }
        
        #inner domain emissions
        hourly.XCO2[j] <-
          suppressWarnings(cellStats(footprint.layer*total.layer,
                                     sum))
        hourly.XCO2.modeled[j] <-
          suppressWarnings(cellStats(footprint.layer*total.layer.mod,
                                     sum))
        
      } #include biospheric emissions
    } #footprint layers (j) loop
    
    #sum each layer's contribution
    XCO2 <- sum(hourly.XCO2); XCO2.modeled <- sum(hourly.XCO2.modeled)
    
    #include the error here
    #how the error is added depends on how it's quantified
    if(!is.na(obs.error) & !regridded) {
      
      #' If an observation is provided in the main script,
      #' *inversion.analysis.r* then the other error method
      #' will be overridden.
      collect.error <- rnorm(1, 0, obs.error)
      
    } else if(is.na(obs.error) & !regridded) {
      
      #mean and s.d. for a normal distribution
      if(!all(is.na(errors[,2])) & !all(is.na(errors[,3]))) {
        
        #add standard deviations in quadrature
        combined.sd <- sqrt(sum(errors[,3]^2))
        
        # Obtain the error value from a normal distribution
        collect.error <- rnorm(1, 0, combined.sd)
        
        #max percentage of value used in a uniform distribution
      } else if((all(is.na(errors[,2])) & all(is.na(errors[,3]))) &
                all(is.na(errors[,5])) & !all(is.na(errors[j,4]))) {
        
        #add percent errors in quadrature
        combined.percent <- sqrt(sum(errors[,4]^2))
        
        #Obtain the error value from a percent error
        collect.error <- runif(1,
                               -1*combined.percent*XCO2,
                               combined.percent*XCO2)
        
        #RMSE used as s.d. in normal distribution (mu = 0)
      } else if((all(is.na(errors[,2])) & all(is.na(errors[,3])) &
                 all(is.na(errors[,4]))) & !all(is.na(errors[,5]))) {
        
        #add the RMSE values in quadrature
        combined.rmse <- sqrt(sum(errors[,5]^2))
        
        #start adding the error term for a normal distribution
        collect.error <- rnorm(1, 0, combined.rmse)
        
      } else{stop('Check error file!')}
    } #grab obs error from R setup file
    
    ### Inner domain - obs.rds
    #add to the XCO2 dataframe with error
    if(!regridded) {
      
      #add the random error to the "true" value
      OBS <- XCO2 + collect.error
      
    } else if(regridded) {
      
      #use the lon/lat to determine what the original observation
      #value is from the original dataframe.
      orig.idx <- which(original.obs$lon == as.numeric(tmp.lon) &
                          original.obs$lat == as.numeric(tmp.lat))
      
      #assign the observation value
      OBS <- original.obs[orig.idx,'xco2']

    }
    
    add.line <-
      data.frame(tmp.date.time,
                 as.numeric(tmp.lon),
                 as.numeric(tmp.lat),
                 XCO2.modeled,
                 XCO2, OBS,
                 footprints[i])
    names(add.line) <- names(XCO2_df)
    
    #add the line to the dataframe
    XCO2_df <- rbind(XCO2_df, add.line)
    
  } #closes the footprint loop
  
  #write the synthetic footprints to a file
  write.csv(XCO2_df,
            file = file.path(dirname(obs.output),
                             'xco2_observations.csv'),
            row.names = FALSE)
  
  #generate obs.rds
  obs.rds_df <- data.frame(site_code = site_code.list,
                           seconds_since_1970_01_01 =
                             paste(as.numeric(start.time)),
                           co2_ppm = XCO2_df$xco2)
  
  if(filter.noise) {
    sig.idx <- which(as.numeric(obs.rds_df[,3]) > obs.error)
    obs.rds_df <- obs.rds_df[sig.idx,]
  }
  
  #format and save the file
  obs.rds <- as.matrix(obs.rds_df)
  saveRDS(obs.rds, file = obs.output)
  
  if(std.bkg == F) {
    
    stop('Non-standard background not currently supported!')

  } else if(std.bkg == T) {
    
    #sets background to zero
    bkg.XCO2_df <- data.frame(time = start.time,
                              avg_CO2 = 0)
    saveRDS(bkg.XCO2_df, file = background.output)
  }
  
}
