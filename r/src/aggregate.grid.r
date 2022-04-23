#' Aggregate values in an N x 4 dataframe where the format is:
#' `lon`, `lat`, `value`, `error`
aggregate.grid <- function(domain, spatial.data, resolution = 0.05) {
  
  if(dim(spatial.data)[2] != 4)
    stop('Check the format of spatial.data!')
  
  spatial.data <- subset(spatial.data,
                         !is.nan(spatial.data[,3]) &
                         !is.null(spatial.data[,3]) &
                         !is.na(spatial.data[,3]))
  
  #determine min/max lons
  min.lons <- domain$xmin; max.lons <- domain$xmax
  break.lons <- seq(min.lons, max.lons, resolution)
  
  #determine min/max lats
  min.lats <- domain$ymin; max.lats <- domain$ymax
  break.lats <- seq(min.lats, max.lats, resolution)
  
  agg_dataframe <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  names(agg_dataframe) <- c('lon', 'lat', 'value', 'error', 'count')
  for(j in 1:(length(break.lons)-1)) {
    for(i in 1:(length(break.lats)-1)) {
      
      tmp_df <-
        subset(spatial.data,
               lon > break.lons[j] & lon < break.lons[j+1] &
                 lat > break.lats[i] & lat < break.lats[i+1])
      
      if(nrow(tmp_df) != 0) {
        tmp.lon <- (break.lons[j] + break.lons[j+1])/2
        tmp.lat <- (break.lats[i] + break.lats[i+1])/2
        tmp.val <- mean(tmp_df[,3])
        if(nrow(tmp_df) >= 3) {
          tmp.error <- sd(tmp_df[,4])/sqrt(nrow(tmp_df))
        } else {
          tmp.error <- mean(tmp_df[,4])
        }
      } else if(nrow(tmp_df) == 0) {
        tmp.lon <- (break.lons[j] + break.lons[j+1])/2
        tmp.lat <- (break.lats[i] + break.lats[i+1])/2
        tmp.val <- NA; tmp.error <- NA
      }
      
      #add the data to the aggregated dataframe
      agg_dataframe <- rbind(agg_dataframe,
                             data.frame(lon = tmp.lon,
                                        lat = tmp.lat,
                                        values = tmp.val,
                                        error = tmp.error,
                                        count = nrow(tmp_df)))
      
    } #close the lats loop
  } #close the lons loop
 
  return(agg_dataframe)
   
}