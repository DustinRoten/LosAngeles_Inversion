convolve_categories <- function(work.dir, footprint.dir, category.dir,
                                output.path) {
  
  #set working directory and load dependencies
  setwd(work.dir); source('r/dependencies.r')
  
  #list all of the available categories from
  categories.nc <- list.files(category.dir, full.names = TRUE)
  
  ODIAC.nc <- grep(categories.nc, pattern = 'odiac')
  ODIAC.nc <- categories.nc[ODIAC.nc]
  odiac.categories <- str_split_fixed(basename(ODIAC.nc),
                                      pattern = '_', n = 2)[,1]
  
  Vulcan.nc <- grep(categories.nc, pattern = 'vulcan')
  Vulcan.nc <- categories.nc[Vulcan.nc]
  vulcan.categories <- str_split_fixed(basename(Vulcan.nc),
                                       pattern = '_', n = 2)[,1]
  
  #Do the categories match?
  if(all(odiac.categories != vulcan.categories))
    stop('Vulcan and ODIAC categories do not match!')
  num.cats <- length(odiac.categories)
  
  foot.nc <- list.files(list.files(footprint.dir,
                                   full.names = TRUE),
                        full.names = TRUE)
  
  category_xco2 <- data.frame(matrix(NA, nrow = 0, ncol = 11))
  names(category_xco2) <- c('SAM', 'lon', 'lat',
                            paste0('Vulcan.', vulcan.categories),
                            paste0('ODIAC.', odiac.categories))
  for(i in 1:length(foot.nc)) {
    
    date.time <- str_split_fixed(basename(foot.nc[i]),
                                 pattern = '_', n = 3)[1]
    date.time.pos <- as.POSIXct(date.time,
                                format = '%Y%m%d%H%M',
                                tz = 'UTC')
    
    lon <- str_split_fixed(basename(foot.nc[i]),
                           pattern = '_', n = 3)[2]
    lon <- as.numeric(lon)
    
    lat <- str_split_fixed(basename(foot.nc[i]),
                           pattern = '_', n = 3)[3]
    lat <- gsub('.nc', '', lat)
    lat <- as.numeric(lat)
    
    footprint <- brick(foot.nc[i])
    
    add.line <-
      data.frame(matrix(NA, nrow = 1, ncol = NCOL(category_xco2)))
    names(add.line) <- names(category_xco2)
    add.line$SAM <- date.time.pos
    add.line$lon <- lon; add.line$lat <- lat
    
    for(j in 1:num.cats) {
      
      #for matching dataframe column names
      which.vulcan.cat <- paste0('Vulcan.', vulcan.categories[j])
      which.odiac.cat <- paste0('ODIAC.', odiac.categories[j])
      
      ODIAC <- brick(ODIAC.nc[j])
      Vulcan <- brick(Vulcan.nc[j])
      
      hourly.values.ODIAC <- hourly.values.Vulcan <- NULL
      for(k in 1:nlayers(footprint)) {
        
        foot.layer <- footprint[[k]]
        layer.name <- names(foot.layer)
        
        eval(parse(text = paste0('odiac.layer <- ODIAC$',
                                 layer.name)))
        eval(parse(text = paste0('vulcan.layer <- Vulcan$',
                                 layer.name)))
        
        hourly.values.ODIAC[k] <-
          suppressWarnings(cellStats(foot.layer*odiac.layer, sum))
        hourly.values.Vulcan[k] <-
          suppressWarnings(cellStats(foot.layer*vulcan.layer, sum))
        
      } #close the hourly convolutions
      
      add.line[, which.odiac.cat] <- sum(hourly.values.ODIAC)
      add.line[, which.vulcan.cat] <- sum(hourly.values.Vulcan)
      
    } #close the category loop
    
    category_xco2 <- rbind(category_xco2, add.line)
    
  }
  
  write.csv(category_xco2,
            file = file.path(output.path, 'category_xco2.csv'),
            row.names = FALSE)
  
}