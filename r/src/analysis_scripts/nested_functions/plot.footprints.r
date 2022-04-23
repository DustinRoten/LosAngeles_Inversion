plot.footprints <- function(site, tz = local.tz, output.path, prior.truth.ext,
                            posterior.ext, plot.output = 'Out/Native_Results',
                            plot.width = 9, plot.height = 9, p.caption = NULL,
                            gg.map = NULL) {
  
  library(ggplot2); library(ggmap); library(viridis)
  library(patchwork); library(ggallin)
  
  #create the directory to save the plots
  if(!dir.exists(plot.output))
    dir.create(plot.output, recursive = TRUE)
  
  #list all of the output directories
  output.dirs <- list.files(output.path, full.names = TRUE)
  for(i in 1:length(output.dirs)) {
    
    #get the path to the averaged footprint
    averaged.footprints <- list.files(file.path(output.dirs[i], prior.truth.ext),
                                      full.names = TRUE,
                                      pattern = 'averaged.footprint.nc')
    foot <- brick(averaged.footprints)
    
    #get the path to the prior raster
    prior.path <- list.files(file.path(output.dirs[i], prior.truth.ext),
                             full.names = TRUE,
                             pattern = 'prior_emiss.nc')
    prior <- brick(prior.path)
    
    #get the path to the truth raster
    truth.path <- list.files(file.path(output.dirs[i], prior.truth.ext),
                             full.names = TRUE,
                             pattern = 'truth_emiss.nc')
    truth <- brick(truth.path)
    
    differences <- prior - truth
    strength <- foot*differences
    
    if(i == 1) {
      summed.foot <- foot
      summed.strength <- strength
    } else if(i > 1) {
      summed.foot <- summed.foot + foot
      summed.strength <- summed.strength + strength
    }
    
  }; removeTmpFiles(h = 0)
  
  #convert the summed rasters to averages
  avg.foot <- (1/length(output.dirs))*summed.foot
  avg.strength <- (1/length(output.dirs))*summed.strength
  remove('summed.foot'); remove('summed.strength')
  
  #convert to a dataframe
  timestep.names <- c(nlayers(avg.foot):1)
  for(i in 1:nlayers(avg.foot)) {
    
    #get each foot layer
    tmp.foot <- crop(avg.foot[[i]], avg.strength[[i]])
    tmp.foot_df <- raster::as.data.frame(tmp.foot, xy = TRUE)
    tmp.foot_df <- subset(tmp.foot_df, tmp.foot_df[,3] != 0)
    names(tmp.foot_df) <- c('lon', 'lat', 'foot')
    
    #get each strength layer
    tmp.strength <- avg.strength[[i]]
    tmp.strength_df <- raster::as.data.frame(tmp.strength, xy = TRUE)
    tmp.strength_df <- subset(tmp.strength_df, tmp.strength_df[,3] != 0)
    names(tmp.strength_df) <- c('lon', 'lat', 'diff.flux')
    
    if(i == 1) {
      foot_df <- data.frame(backward.timestep = timestep.names[i],
                            tmp.foot_df)
      strength_df <- data.frame(backward.timestep = timestep.names[i],
                                tmp.strength_df)
    } else if(i > 1) {
      foot_df <- rbind(foot_df,
                       data.frame(backward.timestep = timestep.names[i],
                                  tmp.foot_df))
      strength_df <- rbind(strength_df,
                           data.frame(backward.timestep = timestep.names[i],
                                      tmp.strength_df))
    }
  } #closes the dataframe building loop

  #Determine discrete values
  strength_df$bin <- NA
  strength_df$bin[strength_df$diff.flux > 0] <- 'Positive'
  strength_df$bin[strength_df$diff.flux < 0] <- 'Negative'
  
  sub.strength_df <- subset(strength_df,
                         diff.flux > 1e-6 | diff.flux < -1e-6)
  
  timestep.contributions <- ggmap(gg.map) +
    ggtitle(expression('Flux Contribution to XCO'[2])) +
    labs('Hours Before SAM') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = sub.strength_df,
                aes(x = lon, y = lat, fill = diff.flux)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name = paste0('Contribution', '\n', '(ppm)'),
                         trans = pseudolog10_trans,
                         limits = c(-0.00001, 0.00001),
                         oob = scales::squish) +
    coord_cartesian() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ backward.timestep)
  ggsave(timestep.contributions,
         filename = file.path(plot.output,
                              'Timestep_Contributions.jpg'),
         device = 'jpg', height = 8.5, width = 8, units = 'in')
  
}