#' Require the sums of the layer totals from the standard xco2
#' values but calculate the sums from the sectors within the
#' function. This function was designed to adapt to the number
#' of supplied sectors.
add.line.sectors <- function(time, lon, lat, obs, prior, prior.err,
                             truth, posterior, posterior.err,
                             sum.foot, foot.file, sector.names,
                             sector.dir) {
  
  #get the list of raster files
  raster.list <- list.files(file.path(sector.dir, 'include', 'sectors'),
                            pattern = 'prior', full.names = TRUE)
  
  prior.sector.totals <-
    data.frame(matrix(NA, nrow = 1, ncol = length(sector.names)))
  names(prior.sector.totals) <- sector.names
  for(i in 1:length(sector.names)) {
    
    raster.file <- grep(gsub(' ', '', sector.names[i]),
                        raster.list, value = TRUE)
    sector.raster <- brick(raster.file)
    
    footprint <- brick(foot.file)
    
    if(all(names(footprint) == names(sector.raster))) {
      sector.total <- sum(cellStats(sector.raster*footprint, sum))
      eval(parse(text = paste0('prior.sector.totals$`',
                               sector.names[i], '` <- ',
                               sector.total)))
    } else {stop('Mismatched timesteps in sector!')}
    
  }
  
  add.line <- cbind(data.frame(time, lon, lat, obs, prior, prior.err,
                               truth, posterior, posterior.err,
                               sum.foot), prior.sector.totals)
  
  return(add.line)
                    
}