make_prior.outer.truth.v2_ncdf4 <- function(site = NULL, odiac.path = NULL,
                                            odiac.background.path = NULL, vulcan.path = NULL,
                                            sector.list = NULL, times = NULL,
                                            inner.extent = NULL, full.extent = NULL,
                                            downscaling.extension = NULL, carma.path = NULL,
                                            obs.dir = NULL, inner.name = NULL, outer.name = NULL,
                                            truth.name = NULL, flux.threshold = 75,
                                            lps.names = c('_LPS')) {
  
  #remove other indicators beginning with '_' in the sector list
  new.sector.list <- grep(x = sector.list$category,
                          pattern = '_', invert = TRUE)
  new.vulcan.sector_df <- sector.list[new.sector.list,]
  
  #get the list of sectors making up the large point sources
  .LPS <- grep(x = sector.list$category, pattern = '_LPS')
  .LPS <- sector.list[.LPS,]
  
  #check the format of the input times
  if(!any(class(times) == 'POSIXt'))
    stop('Input times must be POSIXt format!')
  
  require(ggmap)
  lon.lat <- 
    suppressWarnings(get.lon.lat(site,
                                 dlon = (inner.extent[2]-inner.extent[1])/2,
                                 dlat = (inner.extent[4]-inner.extent[3])/2))
  
  prev.YYYYMM <- 'NA'
  for(i in 1:length(times)) {
    
    #format the timestep
    YYYYMM <- strftime(times[i], format = '%Y%m')
    day_of_week <- lubridate::wday(times[i], week_start = 1)
    hour_of_day <- hour(times[i])
    
    #get ODIAC layer
    if(YYYYMM != prev.YYYYMM) {
      
      #get the ODIAC file for the prior
      ODIAC.FullDomain <- get.odiac(tiff.path = odiac.path,
                                    nc.extent = full.extent,
                                    YYYYMM, convert.units = TRUE)
      ODIAC.monthly <- crop(ODIAC.FullDomain, inner.extent)
      
      #update the prev.YYYYMM flag
      prev.YYYYMM <- YYYYMM
       
    }
      
    #daily ('weekly') factors
    daily.path <- list.files(downscaling.extension,
                             pattern = 'weekly',
                             full.names = TRUE)
    daily.factors <- brick(daily.path)[[day_of_week]]
    daily.factors <- projectRaster(from = daily.factors,
                                   to = ODIAC.monthly,
                                   res = res(ODIAC.monthly),
                                   crs = crs(ODIAC.monthly),
                                   method = 'ngb')
    
    #hourly downscaling
    hourly.path <- list.files(downscaling.extension,
                              pattern = 'hourly',
                              full.names = TRUE)
    hourly.factors <- brick(hourly.path)[[hour_of_day + 1]]
    hourly.factors <- projectRaster(from = hourly.factors,
                                    to = ODIAC.monthly,
                                    res = res(ODIAC.monthly),
                                    crs = crs(ODIAC.monthly),
                                    method = 'ngb')
    
    #create the ODIAC raster below and make sure there are no NAs or NANs
    ODIAC <-
      daily.factors*hourly.factors*ODIAC.monthly
    ODIAC[is.na(ODIAC)] <- 0
    ODIAC[is.nan(ODIAC)] <- 0
    gc()
    
    #remove the large point sources in a separate raster file
    ODIAC.nolps <- ODIAC
    ODIAC.nolps[ODIAC.nolps > flux.threshold] <- NA
    
    #interpolate the NA values in the modified ODIAC raster
    ODIAC.nolps <- focal(ODIAC.nolps, w = matrix(1, nrow = 3, ncol = 3),
                   na.rm = TRUE, fun = mean, NAonly = TRUE)
    ODIAC.nolps[is.na(ODIAC.nolps)] <- 0
    
    #Get the Vulcan large point sources and add them to a single raster.
    for(lps.idx in 1:length(lps.names)) {
      lps.layer <-
        get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                                ref.raster = ODIAC.FullDomain[[1]],
                                category = lps.names[lps.idx],
                                sector_df = .LPS,
                                date.time = times[i],
                                nc.extent = outer.extent,
                                include.maritime = FALSE,
                                include.elec_prod = TRUE)
      if(lps.idx == 1)
        lps.adds <- lps.layer
      if(lps.idx > 1)
        lps.adds <- lps.adds + lps.layer
    }; lps.adds[lps.adds <= flux.threshold] <- 0
    
    #' Get the total Vulcan flux for the timestep.
    #' This raster layer will be used in the final weighting
    #' calculation for the custom ODIAC inventory.
    total.layer <-
      get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                              ref.raster = ODIAC.FullDomain[[1]],
                              category = NA,
                              sector_df = new.vulcan.sector_df,
                              date.time = times[i],
                              nc.extent = outer.extent,
                              include.maritime = FALSE,
                              include.elec_prod = TRUE)
    
    #trim the vulcan files
    ODIAC.area <- ODIAC
    ODIAC.area[ODIAC.area > 0] <- 1
    total.layer <- ODIAC.area*total.layer
    lps.adds <- ODIAC.area*lps.adds
    
    #Do the math
    ODIAC <- (1-(lps.adds/total.layer))*ODIAC.nolps + lps.adds
    
    #determine total.vulcan for the timestep here
    #units are umol/m2/s
    #a reference raster (ref.raster) is supplied for mapping.
    total.vulcan.bkg <-
      get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                              ref.raster = ODIAC.FullDomain[[1]],
                              category = NA, sector_df,
                              date.time = times[i],
                              nc.extent = outer.extent,
                              include.maritime = FALSE,
                              include.elec_prod = TRUE)
    
    #remove remaining ocean-based emissions
    spatial.ODIAC <- ODIAC.FullDomain
    spatial.ODIAC[ODIAC.FullDomain > 0] <- 1
    total.vulcan.bkg <- spatial.ODIAC*total.vulcan.bkg
    
    #get the inner domain to be optimized
    total.vulcan <- crop(total.vulcan.bkg, inner.extent)
    
    #add inner and outer layers to the output stacks
    if(i == 1) {
      
      #save the "truth" extent
      names(ODIAC) <- as.numeric(times)[i]
      downscaled.ODIAC <- ODIAC
      
      #create the null inner area
      inverse.spatial <- total.vulcan; values(inverse.spatial) <- 0
      inverse.spatial <- extend(inverse.spatial, total.vulcan.bkg,
                                value = 1)
      
      #save the outer extent
      vulcan.background <- total.vulcan.bkg*inverse.spatial
      names(vulcan.background) <- as.numeric(times)[i]
      
      #save the truth raster
      names(total.vulcan) <- as.numeric(times)[i]
      vulcan <- total.vulcan
      
    } else if(i > 1) {
      
      #if the process has already started...
      #add the new ODIAC truth layer to the prior stack
      names(ODIAC) <- as.numeric(times)[i]
      downscaled.ODIAC <- addLayer(downscaled.ODIAC, ODIAC)
      
      #create the null inner area
      inverse.spatial <- total.vulcan; values(inverse.spatial) <- 0
      inverse.spatial <- extend(inverse.spatial, total.vulcan.bkg,
                                value = 1)
      
      #save the outer extent
      vulcan.layer <- total.vulcan.bkg*inverse.spatial
      names(vulcan.layer) <- as.numeric(times)[i]
      vulcan.background <-
        addLayer(vulcan.background, vulcan.layer)
      
      #save the prior raster
      names(total.vulcan) <- as.numeric(times)[i]
      vulcan <- addLayer(vulcan, total.vulcan)
    }
    
    #include a progress message
    message <- paste0('\n', 'Calculated prior emissions for time ',
                      times[i], '. ', round(100*i/length(times), 2),
                      '% complete.')
    cat(message, '\n')
    
  } #closes time loop  
  ##############################
  ##############################
  ##############################
  
  
  
  #########################################
  ##### Save the Prior Emissions File #####
  #########################################
  #once all of the layers have been added, remove NAs
  #save the stack
  downscaled.ODIAC[is.na(downscaled.ODIAC)] <- 0
  writeRaster(downscaled.ODIAC,
              filename = file.path(obs.dir, inner.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, inner.name), write = TRUE)
  
  #' Define these dimensions as prescribed by L. Kunik and D. Mallia
  #' in previous code.
  #convert time
  time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                        as.numeric(times),
                        longname = "seconds since R epoch: 1970-01-01 00:00:00")
  ntime <- length(time_dim$vals)
  
  #convert lats
  lat_dim <- ncdim_def("lat", "degrees_north",
                       ncvar_get(nc, 'latitude'),
                       longname = "latitude (center of cell)")
  nlat <- length(lat_dim$vals)
  
  #convert lons
  lon_dim <- ncdim_def("lon", "degrees_east",
                       ncvar_get(nc, 'longitude'),
                       longname = "longitude (center of cell)")
  nlon <- length(lon_dim$vals)
  
  #assign these new dimensions and variables
  prior_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "prior emissions")
  prior_vars <- list(prior_emiss_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  if(length(times) == 1) prior_values <- ncvar_get(nc, "layer")
  if(length(times) > 1) prior_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper footprint value 
  #' (Kunik & Mallia).
  nc_filename <- file.path(obs.dir, inner.name)
  nc_emiss <- nc_create(nc_filename, prior_vars)
  ncvar_put(nc_emiss, prior_emiss_var, prior_values,
            count = c(nlon, nlat, ntime))
  nc_close(nc_emiss)
  ########################################
  ########################################
  ########################################
  
  
  ##############################################
  ##### Save the Background Emissions File #####
  ##############################################
  #once all of the layers have been added, remove the NAs
  #save the stack
  vulcan.background[is.na(vulcan.background)] <- 0
  writeRaster(vulcan.background,
              filename = file.path(obs.dir, outer.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, outer.name), write = TRUE)
  
  #' Define these dimensions as prescribed by L. Kunik and D. Mallia
  #' in previous code.
  #convert time
  time_dim <-
    ncdim_def("time", "seconds_since_1970_01_01",
              as.numeric(times),
              longname = "seconds since R epoch: 1970-01-01 00:00:00")
  ntime <- length(time_dim$vals)
  
  #convert lats
  lat_dim <- ncdim_def("lat", "degrees_north",
                       ncvar_get(nc, 'latitude'),
                       longname = "latitude (center of cell)")
  nlat <- length(lat_dim$vals)
  
  #convert lons
  lon_dim <- ncdim_def("lon", "degrees_east",
                       ncvar_get(nc, 'longitude'),
                       longname = "longitude (center of cell)")
  nlon <- length(lon_dim$vals)
  
  #assign these new dimensions and variables
  outer_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "outer emissions")
  outer_vars <- list(outer_emiss_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  if(length(times) == 1) outer_values <- ncvar_get(nc, "layer")
  if(length(times) > 1) outer_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper footprint value 
  #' (Kunik & Mallia).
  nc_filename <- file.path(obs.dir, outer.name)
  nc_outer <- nc_create(nc_filename, outer_vars)
  ncvar_put(nc_outer, outer_emiss_var, outer_values,
            count = c(nlon, nlat, ntime))
  nc_close(nc_outer)
  ########################################
  ########################################
  ########################################
  
  
  
  ####################################
  ##### Save True Emissions File #####
  ####################################
  #once all of the layers have been added, remove the NAs
  #save the stack
  vulcan[is.na(vulcan)] <- 0
  writeRaster(vulcan,
              filename = file.path(obs.dir, truth.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, truth.name), write = TRUE)
  
  #' Define these dimensions as prescribed by L. Kunik and D. Mallia
  #' in previous code.
  #convert time
  time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                        as.numeric(times),
                        longname = "seconds since R epoch: 1970-01-01 00:00:00")
  ntime <- length(time_dim$vals)
  
  #convert lats
  lat_dim <- ncdim_def("lat", "degrees_north",
                       ncvar_get(nc, 'latitude'),
                       longname = "latitude (center of cell)")
  nlat <- length(lat_dim$vals)
  
  #convert lons
  lon_dim <- ncdim_def("lon", "degrees_east",
                       ncvar_get(nc, 'longitude'),
                       longname = "longitude (center of cell)")
  nlon <- length(lon_dim$vals)
  
  #assign these new dimensions and variables
  truth_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "true emissions")
  truth_vars <- list(truth_emiss_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  if(length(times) == 1) truth_values <- ncvar_get(nc, "layer")
  if(length(times) > 1) truth_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper footprint value 
  #' (Kunik & Mallia).
  nc_filename <- file.path(obs.dir, truth.name)
  nc_truth <- nc_create(nc_filename, truth_vars)
  ncvar_put(nc_truth, truth_emiss_var, truth_values,
            count = c(nlon, nlat, ntime))
  nc_close(nc_truth)
  #####################################
  #####################################
  #####################################
  
} #closes function