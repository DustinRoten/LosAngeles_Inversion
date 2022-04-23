custom.grid <- function(lons, lats, values, layers, resolution = 0.05) {
  
  input_dataframe <- data.frame(lons, lats, values, layers)
  
  #determine min/max lons
  min.lons <- min(input_dataframe$lons)
  max.lons <- max(input_dataframe$lons)
  break.lons <- seq(min.lons, max.lons, resolution)
  
  #determine min/max lats
  min.lats <- min(input_dataframe$lats)
  max.lats <- max(input_dataframe$lats)
  break.lats <- seq(min.lats, max.lats, resolution)

  #determine the layers
  list.layers <- unique(input_dataframe$layers)
  
  agg_dataframe <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(agg_dataframe) <- names(input_dataframe)
  for(k in 1:length(list.layers)) {
    for(j in 1:(length(break.lons)-1)) {
      for(i in 1:(length(break.lats)-1)) {
        
        tmp_df <-
          subset(input_dataframe,
                 lons > break.lons[j] & lons < break.lons[j+1] &
                 lats > break.lats[i] & lats < break.lats[i+1] &
                 layers == list.layers[k])
        
        if(nrow(tmp_df) != 0) {
          tmp.lon <- (break.lons[j] + break.lons[j+1])/2
          tmp.lat <- (break.lats[i] + break.lats[i+1])/2
          tmp.val <- mean(tmp_df$values)
        } else if(nrow(tmp_df) == 0) {
          tmp.lon <- (break.lons[j] + break.lons[j+1])/2
          tmp.lat <- (break.lats[i] + break.lats[i+1])/2
          tmp.val <- 0
        }
        
        #add the data to the aggregated dataframe
        agg_dataframe <- rbind(agg_dataframe,
                               data.frame(lons = tmp.lon,
                                          lats = tmp.lat,
                                          values = tmp.val,
                                          layers = list.layers[k]))
        
      } #close the lats loop
    } #close the lons loop
  } #close the layer loop
  
  return(agg_dataframe)
  
}