plot.averaged.contribution <- function(high.res.path = NULL, gg.map = NULL) {
  
  #list the available SAMs
  SAM.list <- list.files(high.res.path, full.names = TRUE)
  
  #average the footprints from each SAM
  for(i in 1:length(SAM.list)) {
    
    #read in the prior raster and averaged footprints
    prior <- brick(file.path(SAM.list[i], 'include',
                             'prior_emiss.nc'))
    avg.foot <- brick(file.path(SAM.list[i], 'include',
                                'averaged_footprint.nc'))
    
    ppm_contribution <- suppressWarnings(prior*avg.foot)
    
    if(i == 1)
      sum_contribution <- ppm_contribution
    if(i > 1)
      sum_contribution <- sum_contribution + ppm_contribution
    
  }; sum_contribution <- (1/length(SAM.list))*sum_contribution
  
  #convert to dataframe
  timestep <- 1
  for(i in nlayers(sum_contribution):1) {
    
    tmp_df <- raster::as.data.frame(sum_contribution[[i]], xy = TRUE)
    tmp_df$timestep <- timestep
    names(tmp_df) <- c('lon', 'lat', 'Contribution', 'Timestep')
    
    if(i == nlayers(sum_contribution))
      contributions <- tmp_df
    if(i < nlayers(sum_contribution))
      contributions <- rbind(contributions, tmp_df)
    
    #iterate forward one timestep
    timestep <- timestep + 1
    
  }
  
  #identify the timesteps for the output plot
  step.sequence <- seq(1, nlayers(sum_contribution), 3)
  which.rows <- contributions$Timestep %in% step.sequence
  sub.contributions <- contributions[which.rows,]
  sub.contributions <- subset(sub.contributions, Contribution > 0)
  
  ggmap(gg.map) +
    geom_raster(data = sub.contributions,
                aes(x = lon, y = lat, fill = Contribution)) +
    scale_fill_viridis(trans = 'log10') +
    coord_cartesian() +
    facet_wrap(. ~ Timestep)
  
}