obtain_XCO2 <- function(obs.dir, TCCON.bkg.path, OCO3.path, SMUrF.path,
                        sector.list) {
  
  #require libraries
  require(stringr); require(ncdf4)
  
  #read in the tccon background information
  TCCON.bkg <- read.csv(TCCON.bkg.path)
  
  #get the footprint list
  footprint.paths <- file.path(obs.dir, '../footprints')
  
  #nested directories present so list.files() is used twice
  footprint_list <- list.files(footprint.paths, full.names = TRUE)
  footprint_list <- list.files(footprint_list, full.names = TRUE)
  
  #read in SMUrF
  SMUrF <- brick(SMUrF.path)
  
  #create dataframes to store the data
  for(i in 1:length(footprint_list)) {
    
    ###################################
    ### Obtain Soundings from OCO-3 ###
    ###################################
    #get the date, time, lon, and lat from the file name
    foot.string <- str_split_fixed(basename(footprint_list[i]),
                                   pattern = '_', n = 3)
    
    #footprint date and time
    foot.date.time <- as.POSIXct(foot.string[1,1],
                                 format = '%Y%m%d%H%M',
                                 origin = '1970-01-01',
                                 tz = 'UTC')
    
    #footprint lon
    foot.lon <- as.numeric(foot.string[1,2])
    
    #footprint lat
    foot.lat <- as.numeric(gsub('.nc', '', foot.string[1,3]))
    
    #find the corresponding OCO3 file
    OCO3.file <- grep(list.files(OCO3.path, full.names = TRUE),
                      pattern = paste0('_', strftime(foot.date.time,
                                                     format = '%y%m%d',
                                                     tz = 'UTC'), '_'),
                      value = TRUE)
    
    #read in the oco3 file
    OCO3 <- oco3.nc_to_dataframe(OCO3.file)
    
    #find the corresponding OCO-3 sounding
    #quickly calculate the distances to soundings
    dists <- sqrt((OCO3$lon - foot.lon)^2 + (OCO3$lat - foot.lat)^2)
    OCO3.data <- OCO3[which(dists < 1e-4), c('xco2', 'xco2.uncert')]
    ###################################
    ###################################
    
    
    
    #####################################
    ### Begin Generating Modeled XCO2 ###
    #####################################
    #read in the footprint for convolving
    footprint <- brick(footprint_list[i])
    
    #' Convolving all at once can lead to tmp files.
    #' Breaking the process down into a for loop is
    #' slower but tmp files are less likely to be
    #' created.
    foot.names <- names(footprint)
    SMUrF.names <- names(SMUrF)
    for(j in 1:length(foot.names)) {
      
      #get the specific footprint layer
      footprint.layer <- footprint[[j]]
      
      #match the SMUrF layer in case the order is different
      smurf.idx <- which(SMUrF.names == foot.names[j])
      
      #get the SMUrF layer
      SMUrF.layer <- SMUrF[[smurf.idx]]
      
      #convolve and total
      SMUrF.xco2.layer <- cellStats(footprint.layer*SMUrF.layer, sum)
      
      #add the layer's total
      if(j == 1) total.SMUrF.xco2 <- SMUrF.xco2.layer
      if(j > 1) total.SMUrF.xco2 <- total.SMUrF.xco2 + SMUrF.xco2.layer
      
    }; remove('footprint.layer'); remove('SMUrF.layer')
    #####################################
    #####################################
    
    
    
    ########################################
    ### Obtain Modeled XCO2 from Sectors ###
    ########################################
    #get and clean up the sector names
    sector.names <- sector.list[,1]
    sector.names <- gsub(' ', '', sector.names)
    sector.names <- grep(sector.names, pattern = '_',
                         value = TRUE, invert = TRUE)
    
    #obtain the *.nc sector list
    sector.paths <- list.files(file.path(obs.dir, 'sectors'),
                               full.names = TRUE)
    
    #loop through the sector names and add up the XCO2
    sector.xco2 <- data.frame(matrix(NA, nrow = 1,
                                     ncol = length(sector.names)))
    names(sector.xco2) <- sector.names
    for(j in 1:length(sector.names)) {
      
      #get the accompanying sector path and read in the file
      sec.path <- grep(sector.paths, pattern = sector.names[[j]],
                       value = TRUE)
      sector <- brick(sec.path)
      sector.layer.names <- names(sector)
      
      #loop through each timestep of the footprint to convolve
      #(this helps reduces the number of tmp files created)
      for(k in 1:length(foot.names)) {
        
        #get the specific footprint layer
        footprint.layer <- footprint[[k]]
        
        #match the SMUrF layer in case the order is different
        sec.idx <- which(sector.layer.names == foot.names[k])
        
        #get the SMUrF layer
        sector.layer <- sector[[sec.idx]]
        
        #convolve and total
        sector.xco2.layer <-
          suppressWarnings(cellStats(footprint.layer*sector.layer, sum))
        
        #add the layer's total
        if(k == 1) total.sector.xco2 <- sector.xco2.layer
        if(k > 1) total.sector.xco2 <- total.sector.xco2 + sector.xco2.layer
        
      }; remove('footprint.layer'); remove('sector.layer')
      
      #add this to the temporary
      sector.xco2[1,j] <- total.sector.xco2
      
    } #closes the sector loop (j)
    ########################################
    ########################################
    
    
    #start constructing the output dataframe
    add.line <- data.frame(date.time = foot.date.time,
                           lon = foot.lon, lat = foot.lat,
                           OCO3.data, TCCON.bkg,
                           bio = total.SMUrF.xco2,
                           sector.xco2)
    
    #write each line to the output csv file
    if(i == 1) {
      write.table(add.line, file = file.path(obs.dir, 'xco2.csv'),
                  row.names = FALSE, col.names = TRUE, append = FALSE,
                  sep = ',')
    } else if(i > 1) {
      write.table(add.line, file = file.path(obs.dir, 'xco2.csv'),
                  row.names = FALSE, col.names = FALSE, append = TRUE,
                  sep = ',')
    }
    
  } #closes the footprint loop (i)
} #closes the function