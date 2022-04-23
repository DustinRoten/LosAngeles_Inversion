make_prior.outer.truth.v3_ncdf4 <- function(site = NULL, odiac.prior.path = NULL,
                                         odiac.background.path = NULL, vulcan.path = NULL,
                                         sector.list = NULL, times = NULL,
                                         inner.extent = NULL, full.extent = NULL,
                                         downscaling.extension = NULL, carma.path = NULL,
                                         obs.dir = NULL, inner.name = NULL, outer.name = NULL,
                                         truth.name = NULL, bias = NULL, ignore.lps = FALSE,
                                         lps.threshold = 75, lps.names = c('Power Industry',
                                                                          'Manufacturing')) {
  
  #check the format of the input times
  if(!any(class(times) == 'POSIXt'))
    stop('Input times must be POSIXt format!')
  
  if(!dir.exists(file.path(obs.dir, 'sectors')))
    dir.create(file.path(obs.dir, 'sectors'))
  
  require(ggmap)
  lon.lat <- 
    suppressWarnings(get.lon.lat(site,
                                 dlon = (inner.extent[2]-inner.extent[1])/2,
                                 dlat = (inner.extent[4]-inner.extent[3])/2))
  
  for(i in 1:nrow(sector.list)) {
    
    #get the appropriate category
    category <- sector.list$category[i]
    lps.flag <- category %in% lps.names
    
    
    prev.YYYYMM <- 'NA'
    for(j in 1:length(times)) {
      
      #format the timestep
      YYYYMM <- strftime(times[j], format = '%Y%m')
      day_of_week <- lubridate::wday(times[j], week_start = 1)
      hour_of_day <- hour(times[j])
      
      ######################################
      ### Generate the Hourly ODIAC File ###
      ######################################
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
        
      } #closes the ODIAC update loop
      
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
      
      if(j == 1)
        downscaled.ODIAC.background <- ODIAC.background
      if(j > 1)
        downscaled.ODIAC.background <-
        addLayer(downscaled.ODIAC.background, ODIAC.background)
      
      gc()
      ######################################
      ######################################
      
      #determine total.vulcan for the timestep here
      #units are umol/m2/s
      #a reference raster (ref.raster) is supplied for mapping.
      category.vulcan.layer <- get_vulcan_emiss_hourly(vulcan_emiss = vulcan.path,
                                                       ref.raster = ODIAC.prior[[1]],
                                                       category = category,
                                                       sector_df = sector.list,
                                                       date.time = times[j],
                                                       nc.extent = nc.extent,
                                                       include.maritime = FALSE,
                                                       include.elec_prod = TRUE)
      
      #remove remaining ocean-based emissions
      spatial.ODIAC <- ODIAC.prior
      spatial.ODIAC[ODIAC.prior > 0] <- 1
      category.vulcan.layer <- spatial.ODIAC*category.vulcan.layer
      
      #some QA/QC
      category.vulcan.layer[is.na(category.vulcan.layer)] <- 0
      category.vulcan.layer[is.nan(category.vulcan.layer)] <- 0
      
      if(j == 1) {
        
        ### Create the true cat vulcan ###
        true.cat.vulcan <- category.vulcan.layer
        
        ### Create the bias cat vulcan ###
        #conditionals in cases where large point sources are not ignored
        if(!ignore.lps | (ignore.lps & !lps.flag))
          bias.true.cat.vulcan <- bias*category.vulcan.layer
        
        #in case large point sources are ignored
        if(ignore.lps & lps.flag)
          bias.true.cat.vulcan <- category.vulcan.layer
        
      } else if(j > 1) {
        
        ### Add to the true cat vulcan ###
        true.cat.vulcan <-
          addLayer(true.cat.vulcan, category.vulcan.layer)
        
        ### Add the bias cat vulcan ###
        #large point sources NOT ignored
        if(!ignore.lps | (ignore.lps & !lps.flag))
          bias.true.cat.vulcan <-
            addLayer(bias.true.cat.vulcan, bias*category.vulcan.layer)
        
        #large point sources are ignored
        if(ignore.lps & lps.flag)
          bias.true.cat.vulcan <-
            addLayer(bias.true.cat.vulcan, category.vulcan.layer)
      }
      #####################################
      #####################################
      #####################################
    } #closes the timestep loop (j)
    
    
    
    #######################################################
    ##### Save Prior and True Category Emissions File #####
    #######################################################
    #once all of the layers have been added, remove the NAs
    #save the stack
    for(k in 1:2) {
      
      # k == 1 true value
      if(k == 1) {
        #create the file title
        category.title <- paste0('sector.truth_', gsub(' ', '', category),
                                 '.nc')
        writeRaster(true.cat.vulcan,
                    filename = file.path(obs.dir, 'sectors', category.title),
                    format = 'CDF', overwrite = TRUE)
        #open the raster as a NetCDF file and correct the information
        nc <- nc_open(file.path(obs.dir, 'sectors', category.title),
                      write = TRUE)
      } else if(k == 2) {
        # k == 2 bias value
        #create the file title
        category.title <- paste0('sector.prior_', gsub(' ', '', category),
                                 '.nc')
        writeRaster(bias.true.cat.vulcan,
                    filename = file.path(obs.dir, 'sectors', category.title),
                    format = 'CDF', overwrite = TRUE)
        #open the raster as a NetCDF file and correct the information
        nc <- nc_open(file.path(obs.dir, 'sectors', category.title),
                      write = TRUE)
      }
      
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
                                   longname = "category emissions")
      cat_vars <- list(cat_emiss_var)
      
      #get the variable values from the raster.
      #(defaults to "variable")
      if(length(times) == 1) cat_values <- ncvar_get(nc, "layer")
      if(length(times) > 1) cat_values <- ncvar_get(nc, "variable")
      nc_close(nc)
      
      #' Here, default name after `writeRaster()` is 'layer'.
      #' This must be changed to the proper footprint value 
      #' (Kunik & Mallia).
      nc_filename <- file.path(obs.dir, 'sectors', category.title)
      nc_cat <- nc_create(nc_filename, cat_vars)
      ncvar_put(nc_cat, cat_emiss_var, cat_values, count = c(nlon, nlat, ntime))
      nc_close(nc_cat)
    } #close bias loop (k)
    #####################################
    #####################################
    #####################################
    
    #add the category to the total Vulcan
    if(i == 1) {
      total.vulcan <- true.cat.vulcan
      bias.total.vulcan <- bias.true.cat.vulcan
    } else if(i > 1) {
      total.vulcan <- total.vulcan + true.cat.vulcan
      bias.total.vulcan <- bias.total.vulcan + bias.true.cat.vulcan
    }
    
  } #closes the category loop (i)
  total.vulcan[is.na(total.vulcan)] <- 0
  bias.total.vulcan[is.na(bias.total.vulcan)] <- 0
  ##############################
  ##############################
  ##############################
  
  
  
  ############################################################
  ##### Save the prior_emiss.nc and truth_emiss.nc Files #####
  ############################################################
  for(i in 1:2) {
    
    # i == 1 true value
    if(i == 1) {
      writeRaster(total.vulcan,
                  filename = file.path(obs.dir, truth.name),
                  format = 'CDF', overwrite = TRUE)
      #open the raster as a NetCDF file and correct the information
      nc <- nc_open(file.path(obs.dir, truth.name), write = TRUE)
    } else if(i == 2) {
      # i == 2 bias value
      writeRaster(bias.total.vulcan,
                  filename = file.path(obs.dir, inner.name),
                  format = 'CDF', overwrite = TRUE)
      #open the raster as a NetCDF file and correct the information
      nc <- nc_open(file.path(obs.dir, inner.name), write = TRUE)
    }
    
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
    total_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                                 list(lon_dim, lat_dim, time_dim),
                                 longname = "Vulcan emissions")
    total_vars <- list(total_emiss_var)
    
    #get the variable values from the raster.
    #(defaults to "variable")
    if(length(times) == 1) total_values <- ncvar_get(nc, "layer")
    if(length(times) > 1) total_values <- ncvar_get(nc, "variable")
    nc_close(nc)
    
    #' Here, default name after `writeRaster()` is 'layer'.
    #' This must be changed to the proper footprint value 
    #' (Kunik & Mallia).
    
    if(i == 1)
      nc_filename <- file.path(obs.dir, truth.name)
    if(i == 2)
      nc_filename <- file.path(obs.dir, inner.name)
    
    nc_total <- nc_create(nc_filename, total_vars)
    ncvar_put(nc_total, total_emiss_var, total_values, count = c(nlon, nlat, ntime))
    nc_close(nc_total)
  }
  ############################################################
  ############################################################
  ############################################################
  
  
  
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
  
} #closes function
