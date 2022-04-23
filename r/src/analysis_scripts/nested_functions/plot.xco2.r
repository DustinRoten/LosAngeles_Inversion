#generate multi-panel plots of all XCO2 data
plot.xco2 <- function(site = NULL, tz = NULL, output.path = NULL,
                      analysis.ext = NULL, obs.error = NULL, api.key = NULL,
                      plot.output = 'Out', plot.width = 9, plot.height = 9,
                      point.size = 2, reg.caption = NULL, gg.map = NULL) {
  
  require(Metrics)
  
  #create the directory to save the plots
  if(!dir.exists(plot.output))
    dir.create(plot.output, recursive = TRUE)
  
  #list all of the output directories
  output.dirs <- list.files(output.path, full.names = TRUE)
  
  #list all of the "all_xco2.csv" file paths
  all.xco2.path <- 
    list.files(file.path(output.dirs, analysis.ext),
               full.names = TRUE, pattern = 'all_xco2')
  
  #combine all of the csv files
  all.xco2 <- do.call(rbind, lapply(all.xco2.path, read.csv))
  
  ### Prior, Truth, and Posterior Plots ###
  #determine the limits of the plot legends
  max.OR <- min.OR <- F
  max.obs <- max(c(all.xco2$obs, all.xco2$prior, all.xco2$posterior))
  if(max.obs > 4) {max.obs <- 4; max.OR <- T}
  min.obs <- min(c(all.xco2$obs, all.xco2$prior, all.xco2$posterior))
  if(min.obs < 0) {min.obs <- 0; min.OR <- T}
  
  #setup breaks and labels for the plots
  breaks <- c(floor(min.obs):ceiling(max.obs))
  labels <- as.character(breaks)
  if(min.OR)
    labels[1] <- paste0('\u2264', labels[1])
  if(max.OR)
    labels[length(labels)] <- paste0('\u2265', labels[length(labels)])
  #########################################
  
  #determine the limits of plot legends (differences)
  max.diff.OR <- min.diff.OR <- F
  max.diff <- max(c(all.xco2$prior - all.xco2$obs,
                    all.xco2$posterior - all.xco2$obs,
                    all.xco2$posterior - all.xco2$prior))
  if(max.diff > 2) {max.diff <- 2; max.diff.OR <- T}
  min.diff <- min(c(all.xco2$prior - all.xco2$obs,
                    all.xco2$posterior - all.xco2$obs,
                    all.xco2$posterior - all.xco2$prior))
  if(min.diff < -2) {min.diff <- -2; min.diff.OR <- T}
  
  #setup breaks and labels for the difference plots
  diff.breaks <- c(floor(min.diff):ceiling(max.diff))
  diff.labels <- as.character(diff.breaks)
  if(max.diff.OR)
    diff.labels[1] <- paste0('\u2264', diff.labels[1])
  if(min.diff.OR)
    diff.labels[length(diff.labels)] <-
    paste0('\u2265', diff.labels[length(diff.labels)])
  
  ##### Synthetic Observations #####
  #isolate other synth obs
  obs <- all.xco2
  obs$time <- as.POSIXct(obs$time,
                         format = '%Y-%m-%d %H:%M:%S',
                         tz = 'UTC')
  attr(obs, 'tzone') <- local.tz
  
  #plot the synthetic observations
  synth.obs <- ggmap(gg.map) +
    ggtitle(expression(paste('Observed XCO'[2]))) +
    labs(subtitle = site) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = obs,
               shape = 23, size = point.size,
               aes(x = lon, y = lat, fill = obs)) +
    scale_fill_viridis(name =
                         expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                       limits = c(min.obs, max.obs),
                       breaks = breaks, labels = labels,
                       oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5))
  ggsave(synth.obs,
         filename = file.path(plot.output, 'Synthetic_Observations.jpg'),
         device = 'jpg', width = plot.width, height = plot.height,
         units = 'in')
  ##################################
  
  
  
  ##### Modeled Observations #####  
  mod.obs <- all.xco2
  mod.obs$time <- as.POSIXct(mod.obs$time,
                             format = '%Y-%m-%d %H:%M:%S',
                             tz = 'UTC')
  attr(mod.obs$time, 'tzone') <- local.tz
  
  #plot the modeled observations
  modeled.obs <- ggmap(gg.map) +
    ggtitle(expression(paste('Modeled XCO'[2]))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = mod.obs,
               shape = 23, size = point.size,
               aes(x = lon, y = lat, fill = prior)) +
    scale_fill_viridis(name =
                         expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                       limits = c(min.obs, max.obs),
                       breaks = breaks, labels = labels,
                       oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5))
  ggsave(modeled.obs,
         filename = file.path(plot.output, 'Modeled_Observations.jpg'),
         device = 'jpg', width = plot.width, height = plot.height,
         units = 'in')
  ################################
  
  
  
  ##### Plot the difference (modeled - synth) #####
  diff.obs <- all.xco2
  diff.obs$time <- as.POSIXct(diff.obs$time,
                              format = '%Y-%m-%d %H:%M:%S',
                              tz = 'UTC')
  attr(diff.obs$time, 'tzone') <- local.tz
  
  diff1.obs <- ggmap(gg.map) +
    ggtitle(expression(paste('Modeled - Observed XCO'[2]))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = diff.obs, shape = 23, size = point.size,
               aes(x = lon, y = lat, fill = prior - obs)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name =
                           expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                         limits = c(min.diff, max.diff),
                         breaks = diff.breaks, labels = diff.labels,
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5))
  ggsave(diff1.obs,
         filename = file.path(plot.output, 'Modeled.vs.Synthetic.jpg'),
         device = 'jpg', width = plot.width, height = plot.height,
         units = 'in')
  #################################################
  
  
  ##### Posterior Observations #####  
  #plot the modeled observations
  # posterior.obs <- ggmap(map) +
  #   ggtitle(expression(paste('Optimized XCO'[2], ' from Posterior'))) +
  #   labs(subtitle = site) +
  #   xlab('Longitude') +
  #   ylab('Latitude') +
  #   geom_point(data = mod.obs, shape = 23, size = point.size,
  #              aes(x = lon, y = lat, fill = posterior)) +
  #   scale_fill_viridis(name =
  #                        expression(paste(Delta, 'XCO'[2], ' (ppm)')),
  #                      limits = c(min.obs, max.obs),
  #                      breaks = breaks, labels = labels,
  #                      oob = scales::squish) +
  #   facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
  #   theme_classic() +
  #   theme(legend.position = 'bottom',
  #         legend.key.width = unit(0.75, 'in'),
  #         plot.title = element_text(hjust = 0.5))
  # ggsave(posterior.obs,
  #        filename = file.path(plot.output, 'Posterior_Observations.jpg'),
  #        device = 'jpg', width = plot.width, height = plot.height,
  #        units = 'in')
  ################################
  
  
  
  ##### Plot the difference (posterior - synth) #####
  #isolate large synth obs
  diff2.obs <- ggmap(gg.map) +
    ggtitle(expression(paste('Optimized - Observed XCO'[2]))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = diff.obs, shape = 23, size = point.size,
               aes(x = lon, y = lat, fill = posterior - obs)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name =
                           expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                         limits = c(min.diff, max.diff),
                         breaks = diff.breaks, labels = diff.labels,
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5))
  ggsave(diff2.obs,
         filename = file.path(plot.output, 'Posterior.vs.Synthetic.jpg'),
         device = 'jpg', width = plot.width, height = plot.height,
         units = 'in')
  #################################################
  
  
  
  ##### Plot the difference (posterior - prior) #####
  diff3.obs <- ggmap(gg.map) +
    ggtitle(expression(paste('Optimized - Modeled XCO'[2]))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = diff.obs, shape = 23, size = point.size,
               aes(x = lon, y = lat, fill = posterior - prior)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name =
                           expression(paste(Delta, 'XCO'[2], ' (ppm)')),
                         limits = c(min.diff, max.diff),
                         breaks = diff.breaks, labels = diff.labels,
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4) +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5))
  ggsave(diff3.obs,
         filename = file.path(plot.output, 'Posterior.vs.Modeled.jpg'),
         device = 'jpg', width = plot.width, height = plot.height,
         units = 'in')
  
  
  ##### Plot Regression Line #####
  #set the observation limit to remove any anomalously large values here.
  #(A note will be left on the plot)
  obs.limit <- 5 #ppm
  large.obs <- subset(all.xco2, obs > obs.limit)
  data.removed <- round(100*nrow(large.obs)/nrow(all.xco2), 3)
  
  all.xco2$time <- as.POSIXct(all.xco2$time, tz = 'UTC')
  attr(all.xco2$time, 'tzone') <- local.tz
  
  rmse_df <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  names(rmse_df) <- c('time', 'org.rmse', 'new.rmse')
  for(i in 1:length(unique(all.xco2$time))) {
    
    sub.all.xco2 <- subset(all.xco2, time == unique(all.xco2$time)[i])
    time <- unique(all.xco2$time)[i]
    org.rmse <- Metrics::rmse(sub.all.xco2$obs, sub.all.xco2$prior) 
    new.rmse <- Metrics::rmse(sub.all.xco2$obs, sub.all.xco2$posterior)
    
    add.line <- data.frame(time, org.rmse, new.rmse)
    rmse_df <- rbind(rmse_df, add.line)
    
  }
  
  xco2.regression.plot <- ggplot() +
    ggtitle(expression(paste('Prior/Posterior vs. Observed XCO'[2]))) +
    labs(caption = paste0('*Prior in gray, Posterior in black', '\n',
                          '(Removed ', data.removed, '% of data over 5ppm threshold)')) +
    xlab(expression(paste('Observed XCO'[2], ' (ppm)'))) +
    ylab(expression(paste('*Prior/Posterior XCO'[2], ' (ppm)'))) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    
    #prior vs. obs
    geom_errorbar(data = subset(all.xco2, obs <= obs.limit),
                  color = 'black',
                  aes(x = obs,
                      ymin = prior - prior.err,
                      ymax = prior + prior.err)) +
    geom_errorbarh(data = subset(all.xco2, obs <= obs.limit),
                   color = 'black',
                   aes(xmin = obs - obs.error,
                       xmax = obs + obs.error,
                       y = prior)) +
    geom_point(data = subset(all.xco2, obs <= obs.limit),
               fill = 'black', color = 'black', alpha = 0.5,
               aes(x = obs, y = prior)) +

    #posterior vs. obs
    geom_errorbar(data = subset(all.xco2, obs <= obs.limit),
                  color = 'red',
                  aes(x = obs,
                      ymin = posterior - posterior.err,
                      ymax = posterior + posterior.err)) +
    geom_errorbarh(data = subset(all.xco2, obs <= obs.limit),
                   color = 'red',
                   aes(xmin = obs - obs.error,
                       xmax = obs + obs.error,
                       y = posterior)) +
    geom_point(data = subset(all.xco2, obs <= obs.limit),
               shape = 21,
               fill = 'red', color = 'black',
               aes(x = obs, y = posterior)) +
    
    #text/labels
    geom_text(data = rmse_df,
              aes(x = -Inf, y = Inf, hjust = 0, vjust = 1,
                  label = paste('\n',
                                ' Orig. RMSE = ', round(org.rmse, 2),
                                'ppm', '\n',
                                ' New RMSE = ', round(new.rmse, 2),
                                'ppm'))) +
    
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 4)
  
  ggsave(xco2.regression.plot,
         filename = file.path(plot.output, 'XCO2_Regressions.jpg'),
         device = 'jpg', width = plot.width + 1, height = plot.height + 1,
         units = 'in')
  
}