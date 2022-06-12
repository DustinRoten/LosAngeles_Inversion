get_vulcan_emiss_hourly <- function(vulcan_emiss = NULL, category = NA, sector_df = NA,
                                    date.time = NA, loc = 'US', nc.extent = NA,
                                    ref.raster = NULL) {
  
  #load required packages
  require(stringr)
  
  #list all available rasters
  all.rasters <- grep(list.files(vulcan_emiss, full.names = TRUE),
                      pattern = 'total',
                      value = TRUE, invert = TRUE)
  
  #subset by locations
  all.rasters.idx <-
    which(!is.na(str_match(all.rasters, pattern = loc)))
  all.rasters <- all.rasters[all.rasters.idx]
  
  #required year
  year.number <- year(date.time)
  
  #From the identified categories, get the available file paths
  #Find the closest year if current year is not available.
  idx <- grep(pattern = paste0(year.number),
              x = basename(all.rasters))
  available.rasters <- all.rasters[idx]
  if(length(available.rasters) == 0) {
    i <- 0
    while(length(available.rasters) == 0) {
      i <- i+1
      idx <- grep(pattern = paste0(year.number - i),
                  x = basename(all.rasters))
      available.rasters <- all.rasters[idx]
    }
    cat(paste0('Using Vulcan year ', year.number - i, ' instead.'),
        '\n')
    new.year <- year.number - i
  } else {new.year <- year.number}
  remove('i')
  
  #build the available raster dataframe
  all.rasters_df <- data.frame(raster.path = all.rasters)
  
  #combine data from file name
  all.rasters_df <-
    cbind(all.rasters_df,
          as.data.frame(str_split_fixed(basename(all.rasters),
                                        pattern = '\\.', n = 10)))
  names(all.rasters_df) <- c('Path', 'Inventory', 'Version',
                             'Location', 'Temp_Resolution', 'Spatial_Resolution',
                             'Sector', 'Calculation', 'Year', 'Day', 'Format')
  
  #convert 'Day' to numeric
  all.rasters_df$Day <- as.numeric(gsub('d', '', all.rasters_df$Day))
  
  #add additional information (Week of Year and Day)
  #first, convert to dates
  date.vector <- as.Date(all.rasters_df$Day,
                         origin = paste0(new.year, '-01-01'))

  #add the week number in. (Week starts at 0)  
  all.rasters_df$Week <- as.numeric(strftime(date.vector,
                                             format = '%W',
                                             tz = 'UTC'))
  all.rasters_df$Day_of_Week <- wday(date.vector)
  
  #determine the week of the year in question
  week.number <- as.numeric(strftime(date.time, format = '%W',
                                     tz = 'UTC'))
  day.value <- wday(date.time)
  
  #match the rasters
  sub.rasters <- subset(all.rasters_df,
                        Week == week.number & Day_of_Week == day.value)
  
  #in case there are no available inventory days
  if(nrow(sub.rasters) == 0)
    stop(paste0('Error mapping current day to an inventory day!'))
  

  
  ####################################
  ### Determine Category Emissions ###
  ####################################
  sel.cat <- category
  sector.list <- subset(sector_df, category == sel.cat)
  sector.list <- unlist(strsplit(sector.list$Vulcan_SectorList,
                                 split = ' '))
  
  #get the required sector paths
  which.idx <- sub.rasters$Sector %in% sector.list
  req.rasters <- sub.rasters$Path[which.idx]
  
  #retrieve each raster and layer names
  for(i in 1:length(req.rasters)) {
    
    vulcan.raster <- brick(req.rasters[i])
    
    pos.layer.names <- as.POSIXlt(gsub('X', '', names(vulcan.raster)),
                                  format = '%Y.%m.%d.%H.%M.%S',
                                  tz = 'UTC')
    
    #match the hours to detemrine which layer to use
    hour.idx <- which(hour(pos.layer.names) == hour(date.time))
    
    #grab the layer and add it to the stack
    if(i == 1) {
      tmp.layer <- vulcan.raster[[hour.idx]]
      tmp.layer[is.na(tmp.layer)] <- 0 #remove any NA's in the sectors
      vulcan.layer <- tmp.layer
    } else if(i > 1) {
      tmp.layer <- vulcan.raster[[hour.idx]]
      tmp.layer[is.na(tmp.layer)] <- 0 #remove any NA's in the sectors
      vulcan.layer <- vulcan.layer + tmp.layer
    }
    
  }; remove('vulcan.raster') #conserve RAM
  
  #project the raster onto the appropriate grid
  if(class(ref.raster) == "RasterLayer" | class(ref.raster) == "RasterBrick")
    Vulcan <- projectRaster(vulcan.layer, ref.raster)
  
  #cut down the size of the raster to match the domain
  Vulcan <- crop(Vulcan, nc.extent)
  
  ### convert units (code from Dien Wu)
  #convert the unit of CO2 emiss from tonne-C/km2/hr to umol/m2/s
  Vulcan <- Vulcan*(1E6/12*1E6) #First, convert tonne-C to umol-C (= umole-CO2)
  Vulcan <- Vulcan/3600 # convert per hour to per second
  Vulcan <- Vulcan/1E6	# convert per km2 to per m2
  # NOW sel.co2 has unit of umole-C02/m2/s, can be used directly with footprint
  
  #return the emission category
  return(Vulcan)
  
}
  

  

  

