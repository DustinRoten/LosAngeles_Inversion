make_sectors <- function(out.dir, odiac.path, vulcan.path, vulcan.sector.path,
                         nc.extent, odiac.downscaling, vulcan_sectors) {
  
  #create a cateogry directory
  dir.string <- file.path(out.dir, 'include', 'Categories')
  if(!dir.exists(dir.string)) dir.create(dir.string)
  
  #get the prior and truth rasters
  prior <- brick(file.path(out.dir, 'include', 'prior_emiss.nc'))
  truth <- brick(file.path(out.dir, 'include', 'truth_emiss.nc'))
  
  #check that the timesteps are identical for both rasters
  if(any(names(prior) != names(truth)))
    stop('Mismatch in raster names/timesteps!')
  
  #obtain the timesteps in POSIX time
  times <-
    as.POSIXct(as.numeric(gsub('X', '', names(prior))),
               origin = '1970-01-01', tz = 'UTC')
  
  for(cat in vulcan_sectors$category) {
    
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
                                      nc.extent = nc.extent,
                                      YYYYMM, convert.units = TRUE)
        ODIAC.monthly <- crop(ODIAC.FullDomain, nc.extent)
        
        #update the prev.YYYYMM flag
        prev.YYYYMM <- YYYYMM
        
      }
      
      #daily ('weekly') factors
      daily.path <- list.files(odiac.downscaling,
                               pattern = 'weekly',
                               full.names = TRUE)
      daily.factors <- brick(daily.path)[[day_of_week]]
      daily.factors <- projectRaster(from = daily.factors,
                                     to = ODIAC.monthly,
                                     res = res(ODIAC.monthly),
                                     crs = crs(ODIAC.monthly),
                                     method = 'ngb')
      
      #hourly downscaling
      hourly.path <- list.files(odiac.downscaling,
                                pattern = 'hourly',
                                full.names = TRUE)
      hourly.factors <- brick(hourly.path)[[hour_of_day + 1]]
      hourly.factors <- projectRaster(from = hourly.factors,
                                      to = ODIAC.monthly,
                                      res = res(ODIAC.monthly),
                                      crs = crs(ODIAC.monthly),
                                      method = 'ngb')
      
      #create the ODIAC raster below
      ODIAC <-
        daily.factors*hourly.factors*ODIAC.monthly
      gc()
      
      ########################
      ##### Total Vulcan #####
      ########################
      #determine total.vulcan for the timestep here
      #units are umol/m2/s
      #a reference raster (ref.raster) is supplied for mapping.
      total.vulcan <-
        get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                                ref.raster = ODIAC.FullDomain[[1]],
                                category = NA, sector_df,
                                date.time = times[i],
                                nc.extent = nc.extent,
                                include.maritime = FALSE,
                                include.elec_prod = TRUE)
      
      #remove remaining ocean-based emissions
      spatial.ODIAC <- ODIAC.FullDomain
      spatial.ODIAC[ODIAC.FullDomain > 0] <- 1
      total.vulcan <- spatial.ODIAC*total.vulcan
      
      #get the inner domain to be optimized
      total.vulcan <- crop(total.vulcan, nc.extent)
      ########################

      
      
      #########################
      ##### Vulcan Sector #####      
      #########################
      #determine total.vulcan for the timestep here
      #units are umol/m2/s
      #a reference raster (ref.raster) is supplied for mapping.
      cat.vulcan <-
        get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                                ref.raster = ODIAC.FullDomain[[1]],
                                category = cat,
                                sector_df = vulcan_sectors,
                                date.time = times[i],
                                nc.extent = nc.extent,
                                include.maritime = FALSE,
                                include.elec_prod = TRUE)
      
      #remove remaining ocean-based emissions
      spatial.ODIAC <- ODIAC.FullDomain
      spatial.ODIAC[ODIAC.FullDomain > 0] <- 1
      cat.vulcan <- spatial.ODIAC*cat.vulcan
      
      #get the inner domain to be optimized
      cat.vulcan <- crop(cat.vulcan, nc.extent)
      #########################
      
      #add inner and outer layers to the output stacks
      if(i == 1) {
        
        #save the "truth" extent
        names(ODIAC) <- as.numeric(times)[i]
        downscaled.ODIAC <- ODIAC
        
        #save the prior raster
        names(total.vulcan) <- as.numeric(times)[i]
        vulcan <- total.vulcan
        
        #save the category from vulcan
        names(cat.vulcan) <- as.numeric(times)[i]
        category.vulcan <- cat.vulcan
        
        #calculate the category from ODIAC
        cat.odiac <- (cat.vulcan/total.vulcan)*ODIAC
        names(cat.odiac) <- as.numeric(times)[i]
        category.odiac <- cat.odiac
        
      } else if(i > 1) {
        
        #if the process has already started...
        #add the new ODIAC truth layer to the prior stack
        names(ODIAC) <- as.numeric(times)[i]
        downscaled.ODIAC <- addLayer(downscaled.ODIAC, ODIAC)
        
        #save the prior raster
        names(total.vulcan) <- as.numeric(times)[i]
        vulcan <- addLayer(vulcan, total.vulcan)
        
        #save the category from vulcan
        names(cat.vulcan) <- as.numeric(times)[i]
        category.vulcan <- addLayer(category.vulcan, cat.vulcan)
        
        #calculate the category from ODIAC
        cat.odiac <- (cat.vulcan/total.vulcan)*ODIAC
        names(cat.odiac) <- as.numeric(times)[i]
        category.odiac <- addLayer(category.odiac, cat.odiac)
        
      }
      
    } #closes times loop
    
    ##############################################
    ##### Save the Vulcan Cateogry Emissions #####
    ##############################################
    #once all of the layers have been added, remove NAs
    #save the stack
    category.vulcan[is.na(category.vulcan)] <- 0
    cat.name <- paste0(gsub(' ', '', cat), '_vulcan.nc')
    writeRaster(category.vulcan,
                filename = file.path(dir.string, cat.name),
                format = 'CDF', overwrite = TRUE)
    
    #open the raster as a NetCDF file and correct the information
    nc <- nc_open(file.path(dir.string, cat.name), write = TRUE)
    
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
    cat_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                                 list(lon_dim, lat_dim, time_dim),
                                 longname = "cat emissions")
    cat_vars <- list(cat_emiss_var)
    
    #get the variable values from the raster.
    #(defaults to "variable")
    if(length(times) == 1) cat_values <- ncvar_get(nc, "layer")
    if(length(times) > 1) cat_values <- ncvar_get(nc, "variable")
    nc_close(nc)
    
    #' Here, default name after `writeRaster()` is 'layer'.
    #' This must be changed to the proper footprint value 
    #' (Kunik & Mallia).
    nc_filename <- file.path(dir.string, cat.name)
    nc_emiss <- nc_create(nc_filename, cat_vars)
    ncvar_put(nc_emiss, cat_emiss_var, cat_values,
              count = c(nlon, nlat, ntime))
    nc_close(nc_emiss)
    ########################################
    ########################################
    ########################################
    
    
    
    #############################################
    ##### Save the ODIAC Cateogry Emissions #####
    #############################################
    #once all of the layers have been added, remove NAs
    #save the stack
    category.odiac[is.na(category.odiac)] <- 0
    cat.name <- paste0(gsub(' ', '', cat), '_odiac.nc')
    writeRaster(category.odiac,
                filename = file.path(dir.string, cat.name),
                format = 'CDF', overwrite = TRUE)
    
    #open the raster as a NetCDF file and correct the information
    nc <- nc_open(file.path(dir.string, cat.name), write = TRUE)
    
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
    cat_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "cat emissions")
    cat_vars <- list(cat_emiss_var)
    
    #get the variable values from the raster.
    #(defaults to "variable")
    if(length(times) == 1) cat_values <- ncvar_get(nc, "layer")
    if(length(times) > 1) cat_values <- ncvar_get(nc, "variable")
    nc_close(nc)
    
    #' Here, default name after `writeRaster()` is 'layer'.
    #' This must be changed to the proper footprint value 
    #' (Kunik & Mallia).
    nc_filename <- file.path(dir.string, cat.name)
    nc_emiss <- nc_create(nc_filename, cat_vars)
    ncvar_put(nc_emiss, cat_emiss_var, cat_values,
              count = c(nlon, nlat, ntime))
    nc_close(nc_emiss)
    ########################################
    ########################################
    ########################################
    
  } #closes category loop
  
} #closes function