plot.categories <- function(site = NULL, tz = NULL, output.path = NULL,
                            analysis.ext = NULL, api.key = NULL,
                            plot.output = 'Out', plot.width = 9,
                            plot.height = 9, point.size = 2) {
  
  require(ggplot2); require(ggmap)
  register_google(api.key)
  map <- get_map(location = site, maptype = 'satellite', zoom = 10,
                 color = 'bw')
  
  #create the directory to save the plots
  if(!dir.exists('Out')) dir.create('Out')
  
  #list all of the output directories
  output.dirs <- list.files(output.path, full.names = TRUE)
  
  #list all of the "all_xco2.csv" file paths
  cat.xco2.path <- 
    list.files(file.path(output.dirs, analysis.ext),
               full.names = TRUE, pattern = 'category_xco2')
  
  #combine all of the csv files
  cat.xco2 <- do.call(rbind, lapply(cat.xco2.path, read.csv))
  ncats <- ncol(cat.xco2[,4:ncol(cat.xco2)])/2
  
  #separate ODIAC and Vulcan sectors
  col.names <- names(cat.xco2)
  which.odiac <- grep(col.names, pattern = 'ODIAC')
  which.vulcan <- grep(col.names, pattern = 'Vulcan')
  
  ODIAC.categories <- cat.xco2[,c(1:3, which.odiac)]
  Vulcan.categories <- cat.xco2[,c(1:3, which.vulcan)]
  
  melt.ODIAC.categories <-
    data.frame(rep(ODIAC.categories[,1], ncats),
               rep(ODIAC.categories[,2], ncats),
               rep(ODIAC.categories[,3], ncats),
               melt(ODIAC.categories[, 4:ncol(ODIAC.categories)]))
  names(melt.ODIAC.categories) <-
    c(names(cat.xco2)[1:3], 'Category', 'XCO2')
  melt.ODIAC.categories$Category <-
    gsub('ODIAC.', '', melt.ODIAC.categories$Category)
  
  melt.Vulcan.categories <-
    data.frame(rep(Vulcan.categories[,1], ncats),
               rep(Vulcan.categories[,2], ncats),
               rep(Vulcan.categories[,3], ncats),
               melt(Vulcan.categories[, 4:ncol(Vulcan.categories)]))
  names(melt.Vulcan.categories) <-
    c(names(cat.xco2)[1:3], 'Category', 'XCO2')
  melt.Vulcan.categories$Category <-
    gsub('Vulcan.', '', melt.Vulcan.categories$Category)
  
  if(all(unique(melt.ODIAC.categories[,1]) !=
         unique(melt.Vulcan.categories[,1])))
    stop('Mismatched timesteps!')
  
  cat.names <- unique(c(melt.ODIAC.categories$Category,
                        melt.Vulcan.categories$Category))
  
  #determine the limits of the plot legends
  cat.min.OR <- cat.max.OR <- F
  cat.min <- min(cat.xco2[,4:ncol(cat.xco2)])
  if(cat.min < 0) {cat.min <- 0; cat.min.OR <- T}
  cat.max <- max(cat.xco2[,4:ncol(cat.xco2)])
  if(cat.max > 2) {cat.max <- 2; cat.max.OR <- T}
  
  #setup breaks and labels for the plots
  breaks <- c(floor(cat.min):ceiling(cat.max))
  labels <- as.character(breaks)
  if(cat.min.OR)
    labels[1] <- paste0('\u2264', labels[1])
  if(cat.max.OR)
    labels[length(labels)] <- paste0('\u2265', labels[length(labels)])
  
  for(i in 1:length(cat.names)) {
    
    ### First, plot the ODIAC sector ###
    ODIAC.cat.vals <- subset(melt.ODIAC.categories,
                             Category == cat.names[i])

    #ensure that the times in the dataframe are POSIX
    ODIAC.cat.vals[,1] <- as.POSIXct(ODIAC.cat.vals[,1],
                                     tz = 'UTC')
    
    #convert to local time
    attr(ODIAC.cat.vals[,1], 'tzone') <- local.tz
    
    #generate and save the ODIAC plot
    odiac.cat.plot <- ggmap(map) +
      ggtitle(expression(paste('ODIAC Modeled Sectoral XCO'[2]))) +
      labs(subtitle = separate.camelback(cat.names[i])) +
      xlab('Longitude') +
      ylab('Latitude') +
      geom_point(data = ODIAC.cat.vals,
                 shape = 23, size = point.size,
                 aes(x = lon, y = lat, fill = XCO2)) +
      scale_fill_viridis(name =
                           expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                         limits = c(cat.min, cat.max),
                         breaks = breaks, labels = labels,
                         oob = scales::squish) +
      facet_wrap(. ~ strftime(SAM, tz = local.tz), ncol = 4) +
      theme_classic() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(0.75, 'in'),
            plot.title = element_text(hjust = 0.5))
    ggsave(odiac.cat.plot,
           filename = file.path(plot.output,
                                paste0('ODIAC_', cat.names[i], '.jpg')),
           device = 'jpg', width = plot.width, height = plot.height,
           units = 'in')
    
    ### Second, plot the Vulcan sector ###
    Vulcan.cat.vals <- subset(melt.Vulcan.categories,
                              Category == cat.names[i])
    
    #ensure that the times in the dataframe are POSIX
    Vulcan.cat.vals[,1] <- as.POSIXct(Vulcan.cat.vals[,1],
                                      tz = 'UTC')
    
    #convert to local time
    attr(Vulcan.cat.vals[,1], 'tzone') <- local.tz
    
    #generate and save the Vulcan plot
    vulcan.cat.plot <- ggmap(map) +
      ggtitle(expression(paste('Vulcan Modeled Sectoral XCO'[2]))) +
      labs(subtitle = separate.camelback(cat.names[i])) +
      xlab('Longitude') +
      ylab('Latitude') +
      geom_point(data = Vulcan.cat.vals,
                 shape = 23, size = point.size,
                 aes(x = lon, y = lat, fill = XCO2)) +
      scale_fill_viridis(name =
                           expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                         limits = c(cat.min, cat.max),
                         breaks = breaks, labels = labels,
                         oob = scales::squish) +
      facet_wrap(. ~ strftime(SAM, tz = local.tz), ncol = 4) +
      theme_classic() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(0.75, 'in'),
            plot.title = element_text(hjust = 0.5))
    ggsave(vulcan.cat.plot,
           filename = file.path(plot.output,
                                paste0('Vulcan_', cat.names[i], '.jpg')),
           device = 'jpg', width = plot.width, height = plot.height,
           units = 'in')
  }
            
  
}
