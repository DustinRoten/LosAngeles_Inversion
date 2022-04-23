averaged.plot.flux <- function(site, tz = local.tz, output.path, prior.truth.ext,
                               posterior.ext, api.key, plot.output = 'Out',
                               plot.width = 9, plot.height = 9, p.caption = NULL,
                               gg.map = NULL) {
  
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
    
    ### Prior ###
    #average each prior raster over timesteps
    prior.path <- list.files(file.path(SAM.path, prior.truth.ext),
                             full.names = TRUE,
                             pattern = 'prior_emiss.nc')
    prior <- calc(brick(prior.path), mean)
    
    ### Truth ###
    #average each truth raster over timesteps
    truth.path <- list.files(file.path(SAM.path, prior.truth.ext),
                             full.names = TRUE,
                             pattern = 'truth_emiss.nc')
    truth <- calc(brick(truth.path), mean)
    
    ### Posterior ###
    #average each posterior raster over timesteps
    posterior.path <- list.files(file.path(SAM.path, posterior.ext),
                                 full.names = TRUE,
                                 pattern = 'posterior.nc')
    posterior <- calc(brick(posterior.path), mean) #co-variance matrix
    
    ### Average Rasters ###
    if(i == 1) {
      sum.prior <- prior
      sum.truth <- truth
      sum.posterior <- posterior
    } else if(i > 1) {
      sum.prior <- sum.prior + prior
      sum.truth <- sum.truth + truth
      sum.posterior <- sum.posterior + posterior
    }
  }
  
  #average the sums, convert to a dataframe
  #prior
  avg.prior <- (1/length(output.dirs))*sum.prior
  avg.prior_df <- raster::as.data.frame(avg.prior, xy = TRUE)
  names(avg.prior_df) <- c('lon', 'lat', 'flux')
  avg.prior_df <- subset(avg.prior_df, abs(flux) > 0.25)
  
  #truth
  avg.truth <- (1/length(output.dirs))*sum.truth
  avg.truth_df <- raster::as.data.frame(avg.truth, xy = TRUE)
  names(avg.truth_df) <- c('lon', 'lat', 'flux')
  avg.truth_df <- subset(avg.truth_df, abs(flux) > 0.25)
  
  #posterior
  avg.posterior <- (1/length(output.dirs))*sum.posterior
  avg.posterior_df <- raster::as.data.frame(avg.posterior, xy = TRUE)
  names(avg.posterior_df) <- c('lon', 'lat', 'flux')
  avg.posterior_df <- subset(avg.posterior_df, abs(flux) > 0.25)
  
  #combined prior and truth dataframe
  avg.prior.truth_df <- rbind(data.frame(source = 'Prior Flux',
                                         avg.prior_df),
                              data.frame(source = 'True Flux',
                                         avg.truth_df))
  
  #calculate the difference of the averages (prior, truth)
  avg.diff <- avg.prior - avg.truth
  avg.diff_df <- raster::as.data.frame(avg.diff, xy = TRUE)
  names(avg.diff_df) <- c('lon', 'lat', 'diff.flux')
  avg.diff_df <- subset(avg.diff_df, abs(diff.flux) > 0.25)
  
  #calculate the difference of the averages (posterior, prior)
  avg.post.diff <- 100*(avg.posterior - avg.prior)/avg.prior
  avg.post.diff_df <- raster::as.data.frame(avg.post.diff, xy = TRUE)
  names(avg.post.diff_df) <- c('lon', 'lat', 'diff.perc.flux')
  avg.post.diff_df <- subset(avg.post.diff_df, abs(diff.perc.flux) > 0.25)
  
  #plot the prior and truth rasters
  prior.truth.plot <- ggmap(gg.map) +
    ggtitle('Prior and True Flux') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = avg.prior.truth_df,
                aes(x = lon, y = lat, fill = flux)) +
    scale_fill_viridis(trans = 'log10',
                       name = expression(paste(mu, 'mol/m'^2, '/s'))) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    coord_cartesian() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ source, ncol = 2)
  
  #prepare the log10 scale labels
  #determine the limits of plot legends (differences)
  max.diff.OR <- min.diff.OR <- F
  max.diff <- max(avg.diff_df$diff.flux)
  if(max.diff > 1000) {max.diff <- 1000; max.diff.OR <- T}
  min.diff <- min(avg.diff_df$diff.flux)
  if(min.diff < -1000) {min.diff <- -1000; min.diff.OR <- T}
  
  #setup breaks and labels for the difference plots
  diff.breaks <- c(-1000, -100, -10, 0, 10, 100, 1000)
  diff.labels <- as.character(diff.breaks)
  if(max.diff.OR)
    diff.labels[1] <- paste0('\u2264', diff.labels[1])
  if(min.diff.OR)
    diff.labels[length(diff.labels)] <-
    paste0('\u2265', diff.labels[length(diff.labels)])
  
  #plot the difference rasters
  diff.flux.plot <- ggmap(gg.map) +
    ggtitle('Difference in Flux') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = avg.diff_df,
                aes(x = lon, y = lat, fill = diff.flux)) +
    coord_cartesian() +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0, trans = pseudolog10_trans,
                         name = expression(paste(mu, 'mol/m'^2, '/s')),
                         breaks = diff.breaks,
                         labels = diff.labels) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'))
  
  #combine the previous plots
  combo.plot <- (prior.truth.plot / diff.flux.plot)
  ggsave(combo.plot,
         filename = file.path(plot.output, 'PriorTruth_Flux.jpg'),
         device = 'jpg', height = 10, width = 5.5, units = 'in')
  
  #plot the percent difference between the posterior and prior
  perc.diff.plot <- ggmap(gg.map) +
    ggtitle('Percent Change from Prior') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = avg.post.diff_df,
                aes(x = lon, y = lat, fill = diff.perc.flux)) +
    coord_cartesian() +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0, trans = pseudolog10_trans,
                         name = '(%)',
                         breaks = c(-1000, -100, -10, 0, 10, 100, 1000),
                         labels = c('-1000', '-100', '-10', '0',
                                    '10', '100', '1000')) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'))
  
  #save the percent difference plot
  ggsave(perc.diff.plot,
         filename = file.path(plot.output, 'PercCorr_Flux.jpg'),
         device = 'jpg', height = 5, width = 5, units = 'in')
  
} #closes the function
