variogram_flux.err.ls <- function(r1 = NULL, r2 = NULL, cutoff = 30,
                                  bins = 15, plot.label = NULL,
                                  output = NULL) {
  
  require(ggplot2); require(gstat); require(raster)
  require(stringr)
  
  #read in both rasters
  if(class(r1) != 'RasterBrick' & class(r2) != 'RasterBrick')
    r1 <- brick(r1); r2 <- brick(r2)
  
  for(i in 1:nlayers(r1)) {
    if(i == 1) total.r1 <- r1[[i]]
    if(i > 1) total.r1 <- total.r1 + r1[[i]]
  }; total.r1 <- total.r1/nlayers(r1)
  
  for(i in 1:nlayers(r2)) {
    if(i == 1) total.r2 <- r2[[i]]
    if(i > 1) total.r2 <- total.r2 + r2[[i]]
  }; total.r2 <- total.r2/nlayers(r2)
  
  diff <- total.r1 - total.r2
  diff_df <- raster::as.data.frame(diff, xy = TRUE)
  names(diff_df) <- c('lon', 'lat', 'diff_flux')
  
  coordinates(diff_df) <- ~lon+lat
  proj4string(diff_df) <- CRS('+proj=longlat +datum=WGS84')
  
  vgm <- variogram(diff_flux~lon+lat, diff_df,
                   width = round(cutoff / bins), 
                   cutoff = cutoff, cressie = F)
  
  fit.vgm <- fit.variogram(vgm, vgm('Exp'))
  vgm.line <- variogramLine(fit.vgm, cutoff)
  range.val <- paste(round(fit.vgm[2,3], 2))
  
  #has a file extension already been given?
  dev <- str_match(basename(output), 'jpg|png')[1,1]
  dev <- is.na(dev)
  
  ggsave(
    ggplot() +
      ggtitle('Exponential Variogram') +
      xlab('Distance [km]') +
      ylab(expression(paste(gamma, ' [(', mu, 'mol/m'^2, '/s)'^2, ']'))) +
      labs(subtitle = plot.label) +
      geom_text(aes(x = Inf, y = -Inf), hjust = 1, vjust = 0,
                label = paste0('Range = ', round(fit.vgm[2,3], 2), 'km',
                               '\n',
                               'L = ', round(fit.vgm[2,3]/3, 2), 'km',
                           		'\n',
                           		'Bin width = ', round(cutoff/bins, 2), 'km',
                               '\n')) +
      geom_point(data = vgm, aes(x = dist, y = gamma)) +
      geom_line(data = vgm.line, aes(x = dist, y = gamma),
                color = 'blue') +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)),
    
    filename =
      if(is.na(dev)) {paste0(output, '.jpg')} else {output},
    device = 'jpg', width = 5, height = 4, units = 'in'
  )
  
  #range = 3L. L is needed.
  return(fit.vgm[2,3]/3)

}
