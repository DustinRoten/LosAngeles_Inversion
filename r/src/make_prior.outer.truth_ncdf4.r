make_prior.outer.truth_ncdf4 <- function(site = NULL, odiac.prior.path = NULL,
                                         odiac.background.path = NULL, vulcan.path = NULL,
                                         category = NA, sector.list = NULL, times = NULL,
                                         inner.extent = NULL, full.extent = NULL,
                                         downscaling.extension = NULL, carma.path = NULL,
                                         obs.dir = NULL, inner.name = NULL, outer.name = NULL,
                                         truth.name = NULL) {
  
  
  
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
      ODIAC.prior.monthly <- get.odiac(tiff.path = odiac.prior.path,
                                       nc.extent = inner.extent,
                                       YYYYMM, convert.units = TRUE)
      
      #get the ODIAC file for the background
      ODIAC.background.monthly <- get.odiac(tiff.path =
                                              odiac.background.path,
                                            nc.extent = full.extent,
                                            YYYYMM, convert.units = TRUE)
      
      #remove the prior domain from the background (replace with zeros)
      tmp.raster <- ODIAC.prior.monthly; values(tmp.raster) <- 0
      tmp.raster <- extend(x = tmp.raster, y = ODIAC.background.monthly,
                           value = 1)
      ODIAC.background.monthly <- tmp.raster*ODIAC.background.monthly
     
      #update the prev.YYYYMM flag
      prev.YYYYMM <- YYYYMM
       
    }
      
    #daily ('weekly') factors
    daily.path <- list.files(downscaling.extension,
                             pattern = 'weekly',
                             full.names = TRUE)
    daily.factors <- brick(daily.path)[[day_of_week]]
    daily.factors <- projectRaster(from = daily.factors,
                                   to = ODIAC.background.monthly,
                                   res = res(ODIAC.background.monthly),
                                   crs = crs(ODIAC.background.monthly),
                                   method = 'ngb')
    
    #hourly downscaling
    hourly.path <- list.files(downscaling.extension,
                              pattern = 'hourly',
                              full.names = TRUE)
    hourly.factors <- brick(hourly.path)[[hour_of_day + 1]]
    hourly.factors <- projectRaster(from = hourly.factors,
                                    to = ODIAC.background.monthly,
                                    res = res(ODIAC.background.monthly),
                                    crs = crs(ODIAC.background.monthly),
                                    method = 'ngb')
    
    #create the ODIAC prior raster below
    ODIAC.prior <-
      daily.factors*hourly.factors*ODIAC.prior.monthly
    ODIAC.background <-
      daily.factors*hourly.factors*ODIAC.background.monthly
    gc()
    
    #determine total.vulcan for the timestep here
    #units are umol/m2/s
    #a reference raster (ref.raster) is supplied for mapping.
    total.vulcan <- get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                                            ref.raster = ODIAC.prior[[1]],
                                            category, sector_df,
                                            date.time = times[i],
                                            nc.extent = nc.extent,
                                            include.maritime = FALSE,
                                            include.elec_prod = TRUE)
    
    #remove remaining ocean-based emissions
    spatial.ODIAC <- ODIAC.prior
    spatial.ODIAC[ODIAC.prior > 0] <- 1
    total.vulcan <- spatial.ODIAC*total.vulcan
    
    if(is.na(category)) {
      #add inner and outer layers to the output stacks
      if(i == 1) {
        #save the inner extent
        names(ODIAC.prior) <- as.numeric(times)[i]
        downscaled.ODIAC.prior <- ODIAC.prior
        
        #save the outer extent
        #create the null inner area
        names(ODIAC.background) <- as.numeric(times)[i]
        downscaled.ODIAC.background <- ODIAC.background
        
        #save the truth raster
        names(total.vulcan) <- as.numeric(times)[i]
        vulcan.truth <- total.vulcan
        
      } else if(i > 1) {
        #if the process has already started...
        #add the new ODIAC prior layer to the prior stack
        names(ODIAC.prior) <- as.numeric(times)[i]
        downscaled.ODIAC.prior <-
          addLayer(downscaled.ODIAC.prior, ODIAC.prior)
  
        #save the outer extent
        names(ODIAC.background) <- as.numeric(times)[i]
        downscaled.ODIAC.background <-
          addLayer(downscaled.ODIAC.background, ODIAC.background)
        
        #save the truth raster
        names(total.vulcan) <- as.numeric(times)[i]
        vulcan.truth <- addLayer(vulcan.truth, total.vulcan)
      }
      
      #include a progress message
      message <- paste0('\n', 'Calculated prior emissions for time ',
                        times[i], '. ', round(100*i/length(times), 2),
                        '% complete.')
      cat(message, '\n', '\n')
      
    } else if(!is.na(category)) {
      #CHOOSE SECTORS HERE!!
    }
    
  } #closes time loop  
  ##############################
  ##############################
  ##############################
  
  
  
  ########################################
  ##### Save the prior_emiss.nc File #####
  ########################################
  #once all of the layers have been added, remove NAs
  #save the stack
  downscaled.ODIAC.prior[is.na(downscaled.ODIAC.prior)] <- 0
  writeRaster(downscaled.ODIAC.prior,
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
  ncvar_put(nc_emiss, prior_emiss_var, prior_values, count = c(nlon, nlat, ntime))
  nc_close(nc_emiss)
  ########################################
  ########################################
  ########################################
  
  
  ########################################
  ##### Save the outer_emiss.nc File #####
  ########################################
  #once all of the layers have been added, remove the NAs
  #save the stack
  downscaled.ODIAC.background[is.na(downscaled.ODIAC.background)] <- 0
  writeRaster(downscaled.ODIAC.background,
              filename = file.path(obs.dir, outer.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, outer.name), write = TRUE)
  
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
  ncvar_put(nc_outer, outer_emiss_var, outer_values, count = c(nlon, nlat, ntime))
  nc_close(nc_outer)
  ########################################
  ########################################
  ########################################
  
  
  
  #####################################
  ##### Save Truth Emissions File #####
  #####################################
  #once all of the layers have been added, remove the NAs
  #save the stack
  vulcan.truth[is.na(vulcan.truth)] <- 0
  writeRaster(vulcan.truth,
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
  ncvar_put(nc_truth, truth_emiss_var, truth_values, count = c(nlon, nlat, ntime))
  nc_close(nc_truth)
  #####################################
  #####################################
  #####################################
  
} #closes function