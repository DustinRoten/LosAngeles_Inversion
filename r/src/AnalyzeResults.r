#one giant function to analyze the Bayesian inversion's output

if(F) {
  p.table <- p.table[1,]
  home.dir <- p.table$home.dir
  work.dir <- p.table$work.dir
  data.dir <- p.table$data.dir
  site <- p.table$site
  local.tz <- p.table$local.tz
  prior.file.name <- p.table$prior.file.name
  truth.file.name <- p.table$truth.file.name
  posterior.file.name <- p.table$posterior.file.name
  prior.uncert.file.name <- p.table$prior.uncert.file.name
  which.SAM <- p.table$which.SAM
  api.key <- p.table$api.key
  zoom <- p.table$zoom
  analysis.output.dir.name <- p.table$analysis.output.dir.name
}

AnalyzeResults <- function(home.dir, work.dir, data.dir, site, local.tz,
                           prior.file.name, truth.file.name, posterior.file.name,
                           prior.uncert.file.name, posterior.uncert.file.name,
                           observations, which.SAM = NA, api.key = NA, zoom = 9,
                           analysis.output.dir.name) {
  
  #set the working directory and source dependencies
  setwd(work.dir); source('r/dependencies.r')
  
  #get the ggmap for the spatial plots
  register_google(key = api.key)
  map <- get_map(location = site, maptype = 'satellite', zoom = 10)
  
  #list the inversion output directories
  output.dirs <- which.SAM
  
  #check for sector rasters
  if(dir.exists(file.path(output.dirs, 'include', 'sectors'))) {
    sector.list <- read.csv(list.files('ext',
                                       pattern = 'vulcan',
                                       full.names = TRUE))
    n.sectors <- nrow(sector.list)
    
    #prepare dataframes to store results
    xco2_df.sectors <- data.frame(matrix(NA, nrow = 0,
                                         ncol = (10 + n.sectors)))
    names(xco2_df.sectors) <- c('time', 'lon', 'lat', 'obs',
                                'prior', 'prior.err', 'truth',
                                'posterior', 'posterior.err',
                                'sum.foot.density',
                                sector.list$category)
  } else {
    
    #prepare dataframes to store results
    xco2_df <- data.frame(matrix(NA, nrow = 0, ncol = 10))
    names(xco2_df) <- c('time', 'lon', 'lat', 'obs',
                        'prior', 'prior.err', 'truth',
                        'posterior', 'posterior.err',
                        'sum.foot.density')
    
  }
  
  for(i in 1:length(output.dirs)) {
    
    #################################################
    ##### Construct the appropriate directories #####
    #################################################
    #obtain the 'footprints' directory from the output
    foots.dir <- file.path(output.dirs[i], 'footprints')
    
    #obtain the 'include' directory from the output
    include.dir <- file.path(output.dirs[i], 'include')
    
    #obtain the 'out' directory from the output
    inversion.out.dir <- file.path(output.dirs[i], 'inversion/out')
    #################################################
    
    
    
    ##############################################
    ##### Read in the necessary raster files #####
    ##############################################
    #read in the observations
    obs <- read.csv(file.path(include.dir, observations),
                    header = TRUE)
    
    #obtain a list of footprint raster files
    foot.files <- list.files(list.files(foots.dir, full.names = TRUE),
                             full.names = TRUE)
    
    #read in the prior data
    prior <- brick(file.path(include.dir, prior.file.name))
    n.layers <- names(prior)
    
    #read in the truth data
    truth <- brick(file.path(include.dir, truth.file.name))
    
    #read in the posterior data
    posterior <- brick(file.path(inversion.out.dir,
                                 posterior.file.name))
    
    #read in the prior uncertainty data
    #(co-variance matrix)
    prior.uncert <- brick(file.path(inversion.out.dir,
                                    prior.uncert.file.name))
    prior.uncert <- sqrt(prior.uncert)
    
    #read in the posterior uncertainty data
    #(co-variance matrix)
    posterior.uncert <- brick(file.path(inversion.out.dir,
                                        posterior.uncert.file.name))
    posterior.uncert <- sqrt(posterior.uncert)
    ##############################################
    
    #prepare a raster to average footprint information in
    temp.foot <- raster(foot.files[1])
    averaged.footprint <- raster(vals = 0,
                                 crs = crs(temp.foot),
                                 ext = extent(temp.foot),
                                 resolution = res(temp.foot))
    remove('temp.foot')
    raster::origin(averaged.footprint) <- c(0,0)
    
    tmp.layer <- averaged.footprint
    for(x in 1:(nlayers(prior)-1)) {
    
      averaged.footprint <- 
        addLayer(averaged.footprint,
                 tmp.layer)
      
    }; remove('tmp.layer')
    names(averaged.footprint) <- names(prior)
    
    ##########################################
    ##### Generate results in XCO2 space #####
    ##########################################
    for(j in 1:length(foot.files)) {
      
      #message for user
      message(paste0('Output: ', basename(output.dirs[i]), '\n',
                     'Completion: ',
                     round(100*j/length(foot.files), 2), '%',
                     '\n'))
      
      #obtain the date and time from the name string
      name.string <- gsub('.nc', '', basename(foot.files[j]))
      name.string <- str_split_fixed(name.string, pattern = '_', n = 3)
    
      #get the date.time, lon, and lat values
      date.time <- as.POSIXct(name.string[1], format = '%Y%m%d%H%M',
                              tz = 'UTC')
      lon <- as.numeric(name.string[2])
      lat <- as.numeric(name.string[3])
      
      #determine the observed value for this location
      which.obs <- which.min(sqrt((lon - obs$lon)^2 + (lat - obs$lat)^2))
      obs.value <- obs$xco2[which.obs]
      prior.value <- obs$modeled.xco2[which.obs]
      truth.value <- obs$true.xco2[which.obs]
      
      #read in the footprint
      footprint <- brick(foot.files[j])
      
      if(any(names(footprint) != names(prior)))
        stop('Mismatched footprints and prior')
      
      timestep.prior.xco2 <- NULL; timestep.prior.err.xco2 <- NULL
      timestep.truth.xco2 <- NULL; timestep.posterior.xco2 <- NULL
      timestep.posterior.err.xco2 <- NULL; footprint.layer <- NULL
      layer.name = NULL
      for(k in 1:nlayers(footprint)) {
        
        #get each layer name from the footprint brick
        layer.name[k] <- names(footprint)[k]
        
        #get the corresponding layers from the flux bricks
        eval(parse(text =
                     paste0('foot.layer <- footprint$',
                            layer.name[k])))
        eval(parse(text =
                     paste0('prior.layer <- prior$',
                            layer.name[k])))
        eval(parse(text = 
                     paste0('prior.uncert.layer <- prior.uncert$',
                            layer.name[k])))
        eval(parse(text =
                     paste0('truth.layer <- truth$',
                            layer.name[k])))
        eval(parse(text =
                     paste0('posterior.layer <- posterior$',
                            layer.name[k])))
        eval(parse(text =
                     paste0('posterior.uncert.layer <- posterior.uncert$',
                            layer.name[k])))
        
        #convolve each layer and save to the timestep vectors
        #calculate the area of the footprint first
        tmp.foot <- foot.layer; tmp.foot[tmp.foot > 0] <- 1
        foot.area <- raster::area(tmp.foot)
        
        footprint.layer[k] <-
          cellStats(foot.layer, sum)/cellStats(foot.area, sum)
        timestep.prior.err.xco2[k] <-
          suppressWarnings(cellStats(prior.uncert.layer*foot.layer,
                                     sum))
        timestep.posterior.xco2[k] <-
          suppressWarnings(cellStats(posterior.layer*foot.layer, sum))
        timestep.posterior.err.xco2[k] <-
          suppressWarnings(cellStats(posterior.uncert.layer*foot.layer,
                                     sum))
        
        #determine what part of the footprint is within the domain
        cropped.footprint <- crop(foot.layer, averaged.footprint)
        
        #add the value to the hourly.footprint.strength dataframe
        averaged.footprint[[k]] <-
          averaged.footprint[[k]] + cropped.footprint
        
      }; gc() #close the footprint layer loop
      
      if(exists('xco2_df.sectors')) {
        add.line <- 
          add.line.sectors(time = date.time, lon = lon, lat = lat,
                           obs = obs.value, prior = prior.value,
                           prior.err = sqrt(sum(timestep.prior.err.xco2)),
                           truth = truth.value,
                           posterior = sum(timestep.posterior.xco2),
                           posterior.err =
                             sqrt(sum(timestep.posterior.err.xco2)),
                           sum.foot = sum(footprint.layer),
                           foot.file = foot.files[j], 
                           sector.names =
                             names(xco2_df.sectors)[11:length(xco2_df.sectors)],
                           sector.dir = which.SAM)
        xco2_df.sectors <- rbind(xco2_df.sectors, add.line)
      } else if(exists('xco2_df')) {
        add.line <-
          data.frame(time = date.time,
                     lon = lon, lat = lat,
                     obs = obs.value,
                     prior = prior.value,
                     prior.err = sqrt(sum(timestep.prior.err.xco2)),
                     truth = truth.value,
                     posterior = sum(timestep.posterior.xco2),
                     posterior.err = sqrt(sum(timestep.posterior.err.xco2)),
                     sum.foot = sum(footprint.layer))
        xco2_df <- rbind(xco2_df, add.line)
      }
      

      
    } #close the footprint loop
    ##########################################
    
    if(exists('xco2_df')) {
      
      #remove values not considered (not interacting with the domain)
      xco2_df <- subset(xco2_df, prior > 0)
      
      #create the output directory
      if(!dir.exists(file.path(inversion.out.dir, analysis.output.dir.name)))
        dir.create(file.path(inversion.out.dir, analysis.output.dir.name))
      
      #write the csv file to the new directory
      write.csv(xco2_df,
                file = file.path(inversion.out.dir,
                                 analysis.output.dir.name,
                                 'all_xco2.csv'),
                row.names = FALSE)
      
    } else if(exists('xco2_df.sectors')) {
      
      #' Note: cases were prior == 0 are NOT filtered out since these
      #' output vectors need to match the dimensions of the R matrix.
      
      #create the output directory
      if(!dir.exists(file.path(inversion.out.dir, analysis.output.dir.name)))
        dir.create(file.path(inversion.out.dir, analysis.output.dir.name))
      
      #write the csv file to the new directory
      write.csv(xco2_df.sectors,
                file = file.path(inversion.out.dir,
                                 analysis.output.dir.name,
                                 'all_xco2.csv'),
                row.names = FALSE)
      
    }
    
    #complete averaging and name layers
    averaged.footprint <- averaged.footprint/length(foot.files)
    names(averaged.footprint) <- layer.name
    
    ### Save the averaged footprint raster ###
    #save the raster of averaged footprints
    writeRaster(averaged.footprint,
                filename = file.path(include.dir,
                                     'averaged_footprint.nc'),
                overwrite = TRUE)
    
    #open the raster as a NetCDF file and correct the information
    nc <- nc_open(file.path(include.dir, 'averaged_footprint.nc'),
                  write = TRUE)
    
    #' Define these dimensions as prescribed by L. Kunik and D. Mallia
    #' in previous code.
    #convert time
    time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                          as.numeric(gsub('X', '', layer.name)),
                          longname =
                            "seconds since R epoch: 1970-01-01 00:00:00")
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
    foot_var <- ncvar_def("foot", "ppm/(umol m-2 s-1)",
                          list(lon_dim, lat_dim, time_dim),
                          longname = "mean footprint")
    foot_vars <- list(foot_var)
    
    #get the variable values from the raster.
    #(defaults to "variable")
    if(ntime == 1) foot_values <- ncvar_get(nc, "layer")
    if(ntime > 1) foot_values <- ncvar_get(nc, "variable")
    nc_close(nc)
    
    #' Here, default name after `writeRaster()` is 'layer'.
    #' This must be changed to the proper footprint value 
    #' (Kunik & Mallia).
    nc_filename <- file.path(include.dir, 'averaged_footprint.nc')
    nc_foot <- nc_create(nc_filename, foot_vars)
    ncvar_put(nc_foot, foot_var, foot_values, count = c(nlon, nlat, ntime))
    nc_close(nc_foot)
    
  } #close the SAM loop
  
} #closes the function
