spatial.correlation <- function(spatial.error = NULL, vgm.binwidth = 1, vgm.cutoff = 10,
                                included.xco2.errors = NULL, plot.output.path = NULL) {
  
  library(geodist); library(gstat); library(sp)
  library(ggplot2)
  
  #convert the raster brick to a dataframe
  if(class(spatial.error) == 'RasterBrick') {
    raster_df <- raster::as.data.frame(spatial.error, xy = TRUE)
    num.sectors <- ncol(raster_df) - 2
  } else if(class(spatial.error) == 'data.frame') {
    raster_df <- data.frame(spatial.error$time,
                            spatial.error$lon,
                            spatial.error$lat,
                            spatial.error$error)
    names(raster_df) <- c('time', 'lon', 'lat', 'error')
    num.sectors <- 1
    
    #read in the xco2 error
    included.errors <- read.csv(included.xco2.errors)
    xco2_rmse <- sum(included.errors$RMSE^2)
  }
  
  #obtain the unique times associated with the measurements
  times <- unique(raster_df$time)
  
  #create the Q matrix
  Q <- matrix(0, nrow = nrow(raster_df),
              ncol = nrow(raster_df))
  Q.idx <- 1
  
  for(tt in 1:length(times)) {
    
    #obtain the subset SAM
    sub.time <- subset(raster_df, time == times[tt])
    
    #calculate the distance matrix and convert to km
    dists <- geodist::geodist(x = sub.time[,2:3],
                              measure = 'haversine')/1000
    
    lon.lat.values <- sub.time[,2:3]
    names(lon.lat.values) <- c('lon', 'lat')
    
    #calculate the spatial length scale for each sector
    for(i in 4:ncol(raster_df)) {
      
      lon.lat.values$values <- sub.time[,i]
      
      #create the sigma matrix
      sig <- as.matrix(sub.time[,i])
      sig[is.na(sig)] <- 0
      
      if(class(spatial.error) == 'RasterBrick') {
        #remove zero values and NAs
        sub.lon.lat.values <- subset(lon.lat.values,
                                     !is.na(values) &
                                       !is.nan(values) &
                                       values != 0 &
                                       values <= quantile(values, 0.95))
      } else if(class(spatial.error) == 'data.frame') {
        sub.lon.lat.values <- lon.lat.values
        R_ls <- ceiling(3*min(dists[dists != 0]))
      }
      
      #prepare the dataframe for the variogram calculataion
      coordinates(sub.lon.lat.values) <- ~lon+lat
      proj4string(sub.lon.lat.values) <- CRS('+proj=longlat +datum=WGS84')
      
      if(class(spatial.error) == 'RasterBrick') {
        
        vario.type <- 'Exponential'
        #inform the user of the type of fit being used.
        message('Fitting flux data. Exponential Semivariogram.')
        #calculate the variogram, fit the model
        vgm <- variogram(values~lon+lat, data = sub.lon.lat.values,
                         width = vgm.binwidth, 
                         cutoff = vgm.cutoff, cressie = F)
        model.fit <- fit.variogram(vgm, vgm('Exp'))
        
        #plot the variogram fit
        fitted.line <- variogramLine(model.fit, maxdist = vgm.cutoff)
        
        #prepare the plot title
        plot.title.2 <- paste0(plot.title, '\n', 'Averaged Flux - ',
                               names(raster_df)[i])
        
      } else if(class(spatial.error) == 'data.frame') {
        
        vario.type <- 'Gaussian'
        #inform the user of the type of fit being used.
        message('Fitting XCO2 data. Gaussian Semivariogram.')
        #calculate the variogram, fit the model
        vgm <- variogram(values~lon+lat, data = sub.lon.lat.values,
                         width = vgm.binwidth, 
                         cutoff = vgm.cutoff, cressie = F)
        model.fit <- fit.variogram(vgm, vgm('Gau'))
        
        #plot the variogram fit
        fitted.line <- variogramLine(model.fit, maxdist = vgm.cutoff)
        
        #prepare the plot title
        plot.title.2 <- paste0(strftime(times[tt]), '\n', 'XCO2')
        
      }
      
      #plot the variogram
      plot <- ggplot() +
        ggtitle(plot.title.2) +
        labs(subtitle = paste0(vario.type, ' Semivariogram')) +
        xlab('Distance (km)') +
        ylab(expression(paste(gamma, ' (ppm'^2, ')'))) +
        geom_point(data = vgm,
                   aes(x = dist, y = gamma)) +
        geom_line(data = fitted.line,
                  aes(x = dist, y = gamma),
                  color = 'blue') +
        geom_text(aes(x = Inf, y = -Inf,
                      label = paste0('Range = ', round(model.fit[2,3], 2),
                                     'km', '\n',
                                     if(model.fit[2,3]/3 <= vgm.cutoff) {
                                       paste0('L = ', round(model.fit[2,3]/3, 2))
                                     } else {
                                       paste0('Overridden L = ', round(R_ls, 2),
                                              'km', '\n',
                                              'L = ', round(model.fit[2,3]/3, 2))
                                     },
                                     'km', '\n',
                                     'Bin Width = ', vgm.binwidth, 'km', '\n')),
                  hjust = 1, vjust = 0) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
      
      if(class(spatial.error) == 'RasterBrick') {
        ggsave(plot,
               filename = paste0(plot.output, '_', names(raster_df)[i],
                                 '.jpg'),
               device = 'jpg', height = 4, width = 5,
               units = 'in')
      } else if(class(spatial.error) == 'data.frame') {
        ggsave(plot,
               filename = file.path(plot.output.path,
                                    paste0(gsub(' ', '_', times[tt]),
                                           '.jpg')),
               device = 'jpg', height = 4, width = 5, units = 'in')
      }
      
      #get the length scale
      ls <- model.fit[2,3]/3
      
      #perform the spatial decay calculation
      if(class(spatial.error) == 'RasterBrick')
        X <- exp(-dists/ls)
      if(class(spatial.error) == 'data.frame') {
        X <- xco2_rmse*exp(-dists/ls)
      }
      
      #add Q_part to the Q matrix
      Q[Q.idx:((Q.idx-1) + dim(X)[1]),
        Q.idx:((Q.idx-1) + dim(X)[2])] <- X
      Q.idx <- Q.idx + nrow(sub.time)
      
    } #close the sector loop
  } #close the time loop
  
  return(Q)
}