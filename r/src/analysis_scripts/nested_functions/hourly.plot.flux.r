hourly.plot.flux <- function(site, tz = local.tz, output.path = native.path,
                             prior.truth.ext, posterior.ext, api.key,
                             plot.output = output.path, p.caption = agg_line,
                             gg.map = gg.map) {
  
  library(ggplot2); library(ggallin); library(ggmap)
  library(viridis); library(patchwork)
  
  #create the directory to save the plots
  if(!dir.exists(plot.output))
    dir.create(plot.output, recursive = TRUE)
  
  #list all of the output directories
  output.dirs <- list.files(output.path, full.names = TRUE)
  
  #Contruct the averaged prior, truth, and difference rasters
  for(i in 1:length(output.dirs)) {
  
    #get each SAM output
    SAM.path <- output.dirs[i]
    
    ##########################
    ### Obtain the Rasters ###
    ##########################
    ### Get Prior ###
    prior.path <- list.files(file.path(SAM.path, prior.truth.ext),
                             full.names = TRUE,
                             pattern = 'prior_emiss.nc')
    prior <- brick(prior.path)
    
    ### Get Posterior ###
    posterior.path <- list.files(file.path(SAM.path, posterior.ext),
                             full.names = TRUE,
                             pattern = 'posterior.nc')
    posterior <- brick(posterior.path)
    ##########################
    ##########################
    
    
    
    ##############################
    ### Determine the SAM Time ###
    ##############################
    timesteps <- as.numeric(gsub('X', '', names(prior)))
    delta_t <- unique(abs(diff(timesteps)))
    SAM.date.time <- max(timesteps) + delta_t
    SAM.date.time <- as.POSIXct(SAM.date.time, tz = 'UTC',
                                origin = '1970-01-01')
    ##############################
    ##############################
    
    
    
    ##################################################
    ### Calculate the Percent Difference (per SAM) ###
    ##################################################
    Percent_diff <- 100*(posterior - prior)/prior
    Percent_diff.mean <- calc(Percent_diff, mean)
    Percent_diff_df <- raster::as.data.frame(Percent_diff.mean, xy = TRUE)
    names(Percent_diff_df) <- c('lon', 'lat', 'perc.diff.flux')
    if(i == 1)
      per.SAM.diff_df <- data.frame(SAM = SAM.date.time, Percent_diff_df)
    if(i > 1)
      per.SAM.diff_df <- rbind(per.SAM.diff_df,
                               data.frame(SAM = SAM.date.time,
                                          Percent_diff_df))
    ##################################################
    ##################################################
    
    
    
    #######################################
    ### Add the Timesteps from All SAMs ###
    #######################################
    if(i == 1)
      Percent_diff.hr <- Percent_diff
    if(i > 1)
      Percent_diff.hr <- Percent_diff.hr + Percent_diff
    #######################################
    #######################################
    
  } #close the SAM list
  
  
  
  ###############################################
  ### Convert the Timestep Percent Difference ###
  ###############################################
  timestep.names <- c(length(names(Percent_diff.hr)):1)
  for(i in timestep.names) {
    timestep_df <-
      raster::as.data.frame(Percent_diff.hr[[i]], xy = TRUE)
    names(timestep_df) <- c('lon', 'lat', 'perc.diff')
    if(i == max(timestep.names))
      Percent_diff.hr_df <-
        data.frame(backward.timestep = timestep.names[i],
                   timestep_df)
    if(i < max(timestep.names))
      Percent_diff.hr_df <- rbind(Percent_diff.hr_df,
                                  data.frame(backward.timestep =
                                               timestep.names[i],
                                             timestep_df))
  }
  Percent_diff.hr_df <- subset(Percent_diff.hr_df, abs(perc.diff) > 0.25)
  ###############################################
  ###############################################
  
  
  
  #######################
  ### Construct Plots ###
  #######################
  min.label <- paste0('\u2264', '-100')
  max.label <- paste0('\u2265', '100')
  labels <- c(min.label, '-10', '0', '10', max.label)
  
  attr(per.SAM.diff_df$SAM, 'tzone') <- tz
  SAM.flux <- ggmap(gg.map) +
    ggtitle('Percent Change from Prior') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = per.SAM.diff_df,
                aes(x = lon, y = lat, fill = perc.diff.flux)) +
    coord_cartesian() +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0, trans = pseudolog10_trans,
                         name = '(%)',
                         breaks = c(-100, -10, 0, 10, 100),
                         labels = labels,
                         limits = c(-100, 100),
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ strftime(SAM, tz = tz), ncol = 4)
  ggsave(SAM.flux,
         filename = file.path(plot.output, 'SAM_Flux.jpg'),
         device = 'jpg', height = 9, width = 9, units = 'in')
  
  timestep.plot <- ggmap(gg.map) +
    ggtitle('Percent Change by Timestep') +
    labs(subtitle = 'Hours Before SAM') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = Percent_diff.hr_df,
                aes(x = lon, y = lat, fill = perc.diff)) +
    coord_cartesian() +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0, trans = pseudolog10_trans,
                         name = '(%)',
                         breaks = c(-100, -10, 0, 10, 100),
                         labels = labels,
                         limits = c(-100, 100),
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(. ~ backward.timestep, ncol = 5)
  ggsave(timestep.plot,
         filename = file.path(plot.output, 'Timestep_Flux.jpg'),
         device = 'jpg', height = 8.5, width = 8, units = 'in')
    
}