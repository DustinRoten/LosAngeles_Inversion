time.series_flux.err.lt <- function(r1 = NULL, r2 = NULL, cutoff = 30,
                                     bins = 15, plot.label = NULL,
                                     output = NULL) {
  
  require(gstat); require(raster)
  require(stringr)
  
  r1 <- brick(r1); r2 <- brick(r2)
  
  if(nlayers(r1) != nlayers(r2))
    stop('Mismatched raster layers.')
  
  #get the layer names of raster 1
  layer.names <- names(r1)
  
  #' Check for the proper formatting of raster 1 (`r1`) names.
  #' These should be in *epoch time* and UTC.
  layer.names.pos <- as.POSIXct(as.numeric(gsub('X', '', layer.names)),
                                origin = '1970-01-01', tz = 'UTC')
  if(any(is.na(layer.names.pos)))
    stop('Mislabeled raster layers')
  
  #add to a dataframe and re-order the entries for the time series
  timesteps_df <- data.frame(layer.names = layer.names,
                             layer.names.pos = layer.names.pos)
  timesteps_df <- timesteps_df[order(layer.names.pos),]
  
  for.acf_df <- NULL
  for(i in 1:nrow(timesteps_df)) {
    
    #get the layer name for the iteration
    layer.name <- timesteps_df$layer.names[i]
    
    #read in the appropriate layer from rasters 1 and 2.
    eval(parse(text = paste0('r1.layer <- r1$', layer.name)))
    eval(parse(text = paste0('r2.layer <- r2$', layer.name)))
    
    #sum the layer
    sum.r1.layer <- cellStats(r1.layer, sum)
    sum.r2.layer <- cellStats(r2.layer, sum)
    
    #calculate the difference
    for.acf_df[i] <- sum.r1.layer - sum.r2.layer
        
  } #closes the layer
  
  #calculate the autocorrelation function
  out <- acf(for.acf_df, lag.max = nrow(timesteps_df),
             type = 'correlation', plot = FALSE)
  
  #save the plot here
  jpeg(paste0(output, '.jpg'), height = 300)
  plot(out, main = plot.label,
       xlab = 'Lag', ylab = 'ACF')
  dev.off()
  
  #check to see what values are above 2*S.E.
  t.step <- max(which(abs(out$acf) >= 2/sqrt(nrow(timesteps_df))))
  
  return(t.step)
  
} #closes the function