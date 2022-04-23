make_vulcan.sectors <- function(vulcan_emiss = NULL, sector_df = NULL, date.times = NULL,
                                ref.raster = NULL, output.path = NULL) {
  
  #create the output directory for sectors
  if(!dir.exists(file.path(output.path, 'sectors')))
    dir.create(file.path(output.path, 'sectors'))
  
  #list the category names
  cat.names <- sector_df$category
  
  #start building the category by timesteps
  for(i in 1:length(cat.names)) {
    
    for(j in 1:length(date.times)) {
      
      #get the category layer
      cat.layer <- get_vulcan_emiss_hourly(vulcan_emiss = vulcan_emiss,
                                           ref.raster = ref.raster,
                                           category = cat.names[i],
                                           sector_df = sector_df,
                                           date.time = date.times[j],
                                           nc.extent = outer.extent,
                                           include.maritime = FALSE,
                                           include.elec_prod = TRUE)
      
      if(j == 1)
        cat <- cat.layer
      if(j > 1)
        cat <- addLayer(cat, cat.layer)
      
    } #closes the date times
    
    ###########################
    ### Save the Categories ###
    ###########################
    #once all of the layers have been added, remove NAs
    #save the stack
    cat[is.na(cat)] <- 0
    
    #write the file
    file.name <- file.path(output.path, 'sectors',
                           paste0(gsub(' ', '', cat.names[i]),
                                  '.nc'))
    writeRaster(cat, filename = file.name,
                format = 'CDF', overwrite = TRUE)
    
    #open the raster as a NetCDF file and correct the information
    nc <- nc_open(file.name, write = TRUE)
    
    #' Define these dimensions as prescribed by L. Kunik and D. Mallia
    #' in previous code.
    #convert time
    time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                          as.numeric(date.times),
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
    sector_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                                  list(lon_dim, lat_dim, time_dim),
                                  longname = "prior emissions")
    sector_vars <- list(sector_emiss_var)
    
    #get the variable values from the raster.
    #(defaults to "variable")
    if(length(date.times) == 1) sector_values <- ncvar_get(nc, "layer")
    if(length(date.times) > 1) sector_values <- ncvar_get(nc, "variable")
    nc_close(nc)
    
    #' Here, default name after `writeRaster()` is 'layer'.
    #' This must be changed to the proper footprint value 
    #' (Kunik & Mallia).
    nc_filename <- file.name
    nc_emiss <- nc_create(nc_filename, sector_vars)
    ncvar_put(nc_emiss, sector_emiss_var, sector_values,
              count = c(nlon, nlat, ntime))
    nc_close(nc_emiss)
    
  } #closes the category loop
  
  return('Sectors complete!')
  
} #closes the function
