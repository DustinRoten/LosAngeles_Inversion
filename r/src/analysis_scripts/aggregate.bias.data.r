aggregate.bias.data <- function(bias.paths = NULL, inversion.out.ext = 'inversion/out',
                                include.ext = 'include', gg.map = NULL, obs.error = 0,
                                plot.output = 'Out/Bias_Results', tz = local.tz) {
  
  if(!dir.exists(plot.output))
    dir.create(plot.output)
  
  require(stringr); require(raster)
  require(ggplot2); require(viridis)
  
  #calculate lambdas for sectors
  collect.bias.correction(bias.paths, gg.map, plot.output)
  
  #########################
  ### Find the LPS Data ###
  #########################
  LPS.idx <- grep(bias.paths, pattern = '_LPS')
  LPS.paths <- bias.paths[LPS.idx]
  for(i in 1:length(LPS.paths)) {
    
    ##################################
    ### Create the XCO2 data frame ###
    ##################################
    #Get the LPS path and its associated bias.
    LPS.path <- LPS.paths[i]
    bias <- as.numeric(str_split_fixed(basename(LPS.path),
                                       pattern = '_', n = 2)[1])
    
    #Obtain a list of relevant SAMs
    SAM.list <- list.files(LPS.path, full.names = TRUE)
    
    #list the csv files containing 
    csv.list <- list.files(file.path(SAM.list,
                                     inversion.out.ext,
                                     'Analysis'),
                           full.names = TRUE,
                           pattern = '.csv')
    
    #aggregate the data
    xco2 <- do.call('rbind', lapply(csv.list, read.csv))
    
    #add the data to the LPS dataframe
    add.lines <- data.frame(scenario = 'Large Point Sources',
                            bias = bias,
                            xco2)
    if(i == 1) LPS.data <- add.lines
    if(i > 1) LPS.data <- rbind(LPS.data, add.lines)
    ##################################
    ##################################
    
    
    
    #################################
    ### Aggregate the TRUTH files ###
    #################################
    truth.list <- list.files(file.path(SAM.list,
                                       include.ext),
                             pattern = 'truth_emiss.nc',
                             full.names = TRUE)
    for(j in 1:length(truth.list)) {
      
      #grab a raster
      truth.raster <- brick(truth.list[j])
      mean.truth <- calc(truth.raster, mean)
      
      if(j == 1) agg.truth <- mean.truth
      if(j > 1) agg.truth <- agg.truth + mean.truth
      
    }; agg.truth <- agg.truth/length(truth.list)
    #################################
    for(j in 1:length(truth.list)) {
      if(j == 1)
        hourly.truth <- brick(truth.list[j])
      if(j > 1)
        hourly.truth <- hourly.truth + brick(truth.list[j])
    }; hourly.truth <- hourly.truth/length(truth.list)
    #################################
    #################################
    
    
    
    #################################
    ### Aggregate the PRIOR files ###
    #################################
    prior.list <- list.files(file.path(SAM.list,
                                       include.ext),
                             pattern = 'prior_emiss.nc',
                             full.names = TRUE)
    for(j in 1:length(prior.list)) {
      
      #grab a raster
      prior.raster <- brick(prior.list[j])
      mean.prior <- calc(prior.raster, mean)
      
      if(j == 1) agg.prior <- mean.prior
      if(j > 1) agg.prior <- agg.prior + mean.prior
      
    }; agg.prior <- agg.prior/length(prior.list)
    #################################
    for(j in 1:length(prior.list)) {
      if(j == 1)
        hourly.prior <- brick(prior.list[j])
      if(j > 1)
        hourly.prior <- hourly.prior + brick(prior.list[j])
    }; hourly.prior <- hourly.prior/length(prior.list)
    #################################
    #################################
    
    
    
    #####################################
    ### Aggregate the POSTERIOR files ###
    #####################################
    posterior.list <- list.files(file.path(SAM.list,
                                           inversion.out.ext),
                                 pattern = 'posterior.nc',
                                 full.names = TRUE)
    for(j in 1:length(posterior.list)) {
      
      #grab a raster
      posterior.raster <- brick(posterior.list[j])
      mean.posterior <- calc(posterior.raster, mean)
      
      if(j == 1) agg.posterior <- mean.posterior
      if(j > 1) agg.posterior <- agg.posterior + mean.posterior
      
    }; agg.posterior <- agg.posterior/length(posterior.list)
    #####################################
    for(j in 1:length(posterior.list)) {
      if(j == 1)
        hourly.posterior <- brick(posterior.list[j])
      if(j > 1)
        hourly.posterior <-
          hourly.posterior + brick(posterior.list[j])
    }; hourly.posterior <- hourly.posterior/length(posterior.list)
    #####################################
    #####################################

    
        
    # Initial Error Calculations
    Abs.Err <- abs(agg.prior - agg.truth)
    Abs.Err_df <- raster::as.data.frame(Abs.Err, xy = TRUE)
    Abs.Err_df <- subset(Abs.Err_df, layer != 0)
    names(Abs.Err_df) <- c('lon', 'lat', 'Error')
    Abs.Err_df <- data.frame(scenario = 'Large Point Sources',
                             bias = bias,
                             Abs.Err_df)
    
    if(i == 1) LPS_Abs.Err <- Abs.Err_df
    if(i > 1) LPS_Abs.Err <- rbind(LPS_Abs.Err, Abs.Err_df)
    
    # % Diff Calculations
    Percent.Diff <- 100*(agg.posterior - agg.prior)/agg.prior
    Percent.Diff_df <- raster::as.data.frame(Percent.Diff, xy = TRUE)
    Percent.Diff_df <- subset(Percent.Diff_df, layer != 0)
    names(Percent.Diff_df) <- c('lon', 'lat', 'Diff')
    Percent.Diff_df <- data.frame(scenario = 'Large Point Sources',
                                  bias = bias,
                                  Percent.Diff_df)
    
    if(i == 1) LPS_Percent.Diff <- Percent.Diff_df
    if(i > 1) LPS_Percent.Diff <- rbind(LPS_Percent.Diff, Percent.Diff_df)
    
    # Hourly % Diff Calculations
    Percent.Diff.Hourly <-
      100*(hourly.posterior - hourly.prior)/hourly.prior
    for(j in 1:nlayers(Percent.Diff.Hourly)) {
      
      #calculate the backwards timestep
      backwards.timestep <-
        nlayers(Percent.Diff.Hourly)-(j-1)
      
      if(i == 1 & j == 1) {
        Percent.Diff.Hourly_df <-
          raster::as.data.frame(Percent.Diff.Hourly[[j]],
                                xy = TRUE)
        
        Percent.Diff.Hourly_df <-
          data.frame(scenario = 'Large Point Source Error',
                     bias = bias,
                     backwards.timestep,
                     Percent.Diff.Hourly_df)
        
        names(Percent.Diff.Hourly_df) <-
          c('scenario', 'bias', 'timestep',
            'longitude', 'latitude', 'Diff')
      } else if(i > 1 | j > 1) {
        add.lines <-
          data.frame(scenario = 'Large Point Source Error',
                     bias = bias,
                     backwards.timestep,
                     raster::as.data.frame(Percent.Diff.Hourly[[j]],
                                           xy = TRUE))
        
        names(add.lines) <-
          c('scenario', 'bias', 'timestep',
            'longitude', 'latitude', 'Diff')
        
        Percent.Diff.Hourly_df <-
          rbind(Percent.Diff.Hourly_df, add.lines)
      } # closes iterations where i > 1
    } #closes hourly timesteps
  } #closes the non-LPS loop
  LPS.Percent.Diff.Hourly_df <- Percent.Diff.Hourly_df
  remove('Percent.Diff.Hourly_df')
  
  #############################
  ### Find the Non-LPS Data ###
  #############################
  noLPS.idx <- grep(bias.paths, pattern = '_noLPS')
  noLPS.paths <- bias.paths[noLPS.idx]
  for(i in 1:length(noLPS.paths)) {
    
    ##################################
    ### Create the XCO2 data frame ###
    ##################################
    #Get the LPS path and its associated bias.
    noLPS.path <- noLPS.paths[i]
    bias <- as.numeric(str_split_fixed(basename(noLPS.path),
                                       pattern = '_', n = 2)[1])
    
    #Obtain a list of relevant SAMs
    SAM.list <- list.files(noLPS.path, full.names = TRUE)
    
    #list the csv files containing
    #aggregate the data
    csv.list <- list.files(file.path(SAM.list,
                                     inversion.out.ext,
                                     'Analysis'),
                           full.names = TRUE,
                           pattern = '.csv')
    xco2 <- do.call('rbind', lapply(csv.list, read.csv))
    
    #add the data to the LPS dataframe
    add.lines <- data.frame(scenario = 'No Large Point Source Error',
                            bias = bias,
                            xco2)
    if(i == 1) noLPS.data <- add.lines
    if(i > 1) noLPS.data <- rbind(noLPS.data, add.lines)
    ##################################
    ##################################
    
    
    
    #################################
    ### Aggregate the TRUTH files ###
    #################################
    truth.list <- list.files(file.path(SAM.list,
                                       include.ext),
                             pattern = 'truth_emiss.nc',
                             full.names = TRUE)
    for(j in 1:length(truth.list)) {
      
      #grab a raster
      truth.raster <- brick(truth.list[j])
      mean.truth <- calc(truth.raster, mean)
      
      if(j == 1) agg.truth <- mean.truth
      if(j > 1) agg.truth <- agg.truth + mean.truth
      
    }; agg.truth <- agg.truth/length(truth.list)
    #################################
    for(j in 1:length(truth.list)) {
      if(j == 1)
        hourly.truth <- brick(truth.list[j])
      if(j > 1)
        hourly.truth <- hourly.truth + brick(truth.list[j])
    }; hourly.truth <- hourly.truth/length(truth.list)
    #################################
    #################################
    
    
    
    #################################
    ### Aggregate the PRIOR files ###
    #################################
    prior.list <- list.files(file.path(SAM.list,
                                       include.ext),
                             pattern = 'prior_emiss.nc',
                             full.names = TRUE)
    for(j in 1:length(prior.list)) {
      
      #grab a raster
      prior.raster <- brick(prior.list[j])
      mean.prior <- calc(prior.raster, mean)
      
      if(j == 1) agg.prior <- mean.prior
      if(j > 1) agg.prior <- agg.prior + mean.prior
      
    }; agg.prior <- agg.prior/length(prior.list)
    #################################
    #################################
    for(j in 1:length(prior.list)) {
      if(j == 1)
        hourly.prior <- brick(prior.list[j])
      if(j > 1)
        hourly.prior <- hourly.prior + brick(prior.list[j])
    }; hourly.prior <- hourly.prior/length(prior.list)
    #################################
    
    
    
    #####################################
    ### Aggregate the POSTERIOR files ###
    #####################################
    posterior.list <- list.files(file.path(SAM.list,
                                           inversion.out.ext),
                                 pattern = 'posterior.nc',
                                 full.names = TRUE)
    for(j in 1:length(posterior.list)) {
      
      #grab a raster
      posterior.raster <- brick(posterior.list[j])
      mean.posterior <- calc(posterior.raster, mean)
      
      if(j == 1) agg.posterior <- mean.posterior
      if(j > 1) agg.posterior <- agg.posterior + mean.posterior
      
    }; agg.posterior <- agg.posterior/length(posterior.list)
    #####################################
    for(j in 1:length(posterior.list)) {
      if(j == 1)
        hourly.posterior <- brick(posterior.list[j])
      if(j > 1)
        hourly.posterior <-
          hourly.posterior + brick(posterior.list[j])
    }; hourly.posterior <- hourly.posterior/length(posterior.list)
    #####################################
    #####################################
    
    # % Error Calculations
    Abs.Err <- abs(agg.posterior - agg.truth)
    Abs.Err_df <- raster::as.data.frame(Abs.Err, xy = TRUE)
    Abs.Err_df <- subset(Abs.Err_df, layer != 0)
    names(Abs.Err_df) <- c('lon', 'lat', 'Error')
    Abs.Err_df <- data.frame(scenario = 'No Large Point Source Error',
                             bias = bias,
                             Abs.Err_df)
    
    if(i == 1) noLPS_Abs.Err <- Abs.Err_df
    if(i > 1) noLPS_Abs.Err <- rbind(noLPS_Abs.Err,
                                     Abs.Err_df)
    
    # % Diff Calculations
    Percent.Diff <- 100*(agg.posterior - agg.prior)/agg.prior
    Percent.Diff_df <- raster::as.data.frame(Percent.Diff, xy = TRUE)
    Percent.Diff_df <- subset(Percent.Diff_df, layer != 0)
    names(Percent.Diff_df) <- c('lon', 'lat', 'Diff')
    Percent.Diff_df <- data.frame(scenario = 'No Large Point Source Error',
                                  bias = bias,
                                  Percent.Diff_df)
    
    if(i == 1) noLPS_Percent.Diff <- Percent.Diff_df
    if(i > 1) noLPS_Percent.Diff <- rbind(noLPS_Percent.Diff,
                                          Percent.Diff_df)
    
    # Hourly % Diff Calculations
    Percent.Diff.Hourly <-
      100*(hourly.posterior - hourly.prior)/hourly.prior
    for(j in 1:nlayers(Percent.Diff.Hourly)) {
      
      #calculate the backwards timestep
      backwards.timestep <-
        nlayers(Percent.Diff.Hourly)-(j-1)
      
      if(i == 1 & j == 1) {
        Percent.Diff.Hourly_df <-
          raster::as.data.frame(Percent.Diff.Hourly[[j]],
                                xy = TRUE)
        
        Percent.Diff.Hourly_df <-
          data.frame(scenario = 'No Large Point Source Error',
                     bias = bias,
                     backwards.timestep,
                     Percent.Diff.Hourly_df)
        
        names(Percent.Diff.Hourly_df) <-
          c('scenario', 'bias', 'timestep',
            'longitude', 'latitude', 'Diff')
      } else if(i > 1 | j > 1) {
        add.lines <-
          data.frame(scenario = 'No Large Point Source Error',
                     bias = bias,
                     backwards.timestep,
                     raster::as.data.frame(Percent.Diff.Hourly[[j]],
                                           xy = TRUE))
        
        names(add.lines) <-
          c('scenario', 'bias', 'timestep',
            'longitude', 'latitude', 'Diff')
        
        Percent.Diff.Hourly_df <-
          rbind(Percent.Diff.Hourly_df, add.lines)
      } #closes iterations were i > 1
    } #closes hourly timesteps
  } #closes the non-LPS loop
  noLPS.Percent.Diff.Hourly_df <- Percent.Diff.Hourly_df
  remove('Percent.Diff.Hourly_df')
  
  #put separate dataframes together
  all.data <- rbind(LPS.data, noLPS.data)
  err.data <- rbind(LPS_Abs.Err, noLPS_Abs.Err)
  diff.data <- rbind(LPS_Percent.Diff, noLPS_Percent.Diff)
  hourly.diff.data <- rbind(LPS.Percent.Diff.Hourly_df,
                            noLPS.Percent.Diff.Hourly_df)
  
  
  
  ##############################################################
  ### Plot the comparsions between modeled and observed XCO2 ###
  ##############################################################
  Bias.Correlation <- ggplot() +
    ggtitle('Induced Biases') +
    labs(subtitle = 'Proportion of True Emissions') +
    xlab(expression(paste('Observed XCO'[2], ' (ppm)'))) +
    ylab(expression(paste('Modeled XCO'[2], ' (ppm)'))) +
    geom_point(data = all.data,
               aes(x = obs, y = truth, fill = 'True'),
               shape = 21, color = 'black', alpha = 0.25) +
    geom_point(data = all.data,
               aes(x = obs, y = posterior, fill = 'Posterior'),
               shape = 21, color = 'black', alpha = 0.5) +
    geom_point(data = all.data,
               aes(x = obs, y = prior, fill = 'Prior'),
               shape = 21, color = 'black', alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    scale_fill_manual(name = expression(paste('CO'[2], ' Emission Source')),
                      values = c('True' = 'gray',
                                 'Prior' = 'blue',
                                 'Posterior' = 'green')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom') +
    facet_grid(bias ~ scenario)
  
  #save the plot
  ggsave(Bias.Correlation,
         filename = file.path(plot.output,
                              paste0('XCO2_Plot.jpg')),
         height = 5.5, width = 5.5, units = 'in', device = 'jpg')
  ##############################################################
  ##############################################################
  
  
  
  ##########################
  ### Spatial XCO2 Plots ###
  ##########################
  #list the different scenarios that require plotting
  scenario_df <- distinct(all.data[,1:2])
  for(i in 1:nrow(scenario_df)) {
    
    #construct the subset
    sub.data <- subset(all.data,
                       scenario == scenario_df$scenario[i] &
                         bias == scenario_df$bias[i])
    
    #setup breaks and labels for the difference plots
    min.value <- min(sub.data$posterior - sub.data$prior)
    max.value <- max(sub.data$posterior - sub.data$prior)
    if(max.value - min.value > 1)
      breaks <- seq(ceiling(min.value/2), round(max.value/2, 1), 0.25)
    if(max.value - min.value < 1) {
      breaks <- seq(min.value, max.value, (max.value-min.value)/5)
      breaks <- round(breaks, 2)
    }
    labels <- as.character(breaks)
    
    if(min.value < breaks[1])
      labels[1] <- paste0('\u2264', labels[1])
    if(max.value > breaks[length(breaks)])
      labels[length(labels)] <-
      paste0('\u2265', labels[length(labels)])
    
    scenario.iteration <- ggmap(gg.map) +
      ggtitle(paste0('Correction of Induced Biases', '\n',
                     scenario_df[i,1], ' Included')) +
      labs(subtitle = paste0('Proportion of True Emissions: ',
                             scenario_df[i,2])) +
      xlab('Longitude') +
      ylab('Latitude') +
      geom_point(data = sub.data, shape = 23, color = 'black',
                 aes(x = lon, y = lat, fill = posterior - prior)) +
      scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                           midpoint = 0,
                           limits = c(breaks[1],
                                      breaks[length(breaks)]),
                           breaks = breaks, labels = labels,
                           oob = scales::squish,
                           name = expression(paste('Difference in XCO'[2],
                                                   ' (ppm)'))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'bottom',
            legend.key.width = unit(0.75, 'in')) +
      facet_wrap(. ~ time, ncol = 4)
    
    #save each plot
    if(length(which(breaks == 0)) <= 1)
      ggsave(scenario.iteration,
             filename = file.path(plot.output,
                                  paste0('XCO2_',
                                         gsub(' ', '', scenario_df[i,1]),
                                         '_', scenario_df[i,2], '.jpg')),
             height = 8, width = 8, units = 'in', device = 'jpg')
  }
  
  #################################
  ### Plot the Flux Differences ###
  #################################
  Bias.AbsErr <- ggmap(gg.map) +
    ggtitle(expression(paste('CO'[2], ' Emission Absolute Differences'))) +
    labs(subtitle = 'Plotted with Maximum Bias',
         caption = expression(paste('>1', mu, 'mol/m'^2, '/s Only'))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = subset(err.data, Error > 1 &
                                bias == min(err.data$bias)),
                aes(x = lon, y = lat, fill = Error)) +
    scale_fill_viridis(name = expression(paste(mu, 'mol/m'^2, '/s')),
                       option = 'inferno', trans = 'log10') +
    coord_cartesian() +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(. ~ scenario)
  
  #save the plot
  ggsave(Bias.AbsErr,
         filename = file.path(plot.output,
                              paste0('Flux_AbsErr.jpg')),
         height = 5, width = 7, units = 'in', device = 'jpg')
  #################################
  #################################
  
  
  
  ###############################
  ### Plot Hourly Corrections ###
  ###############################
  hourly.diff.data$Bins <- NA
  hourly.diff.data$Bins[hourly.diff.data$Diff < 0 &
                          hourly.diff.data$Diff >= -5] <-
    '[-5,0)'
  hourly.diff.data$Bins[hourly.diff.data$Diff > 0 &
                          hourly.diff.data$Diff <= 5] <-
    '(0,5]'
  hourly.diff.data$Bins[hourly.diff.data$Diff > 5 &
                          hourly.diff.data$Diff <= 15] <-
    '(5,15]'
  hourly.diff.data$Bins[hourly.diff.data$Diff > 15 &
                          hourly.diff.data$Diff <= 30] <-
    '(15,30]'
  hourly.diff.data$Bins[hourly.diff.data$Diff > 30] <-
    '(30,Inf)'
  
  hourly.diff.data$Bins <-
    factor(hourly.diff.data$Bins,
           levels = c('[-5,0)', '(0,5]', '(5,15]', '(15,30]', '(30,Inf)'))
  Disc.Data <- subset(hourly.diff.data, !is.na(Bins))
  
  scenarios <- unique(Disc.Data$scenario)
  select.timesteps <- c(1, 5, 9, 13, 17)
  for(i in 1:length(scenarios)) {
    hourly.diff <- ggmap(gg.map) +
      ggtitle(paste0('Correction of Induced Biases', '\n',
                     scenarios[i],' Included')) +
      labs(subtitle = 'Backwards Timestep and Proportion of True Emissions') +
      xlab('Longitude') +
      ylab('Latitude') +
      geom_raster(data = Disc.Data[Disc.Data$timestep %in% select.timesteps &
                                     !is.na(Disc.Data$Bins) &
                                     abs(Disc.Data$Diff) > 0 &
                                     Disc.Data$scenario == scenarios[i],],
                  aes(x = longitude, y = latitude, fill = Bins)) +
      scale_fill_manual(values = c('lightblue', 'yellow', 'orange', 'red', 'purple'),
                        name = '(%)') +
      coord_cartesian() +
      theme_classic() +
      facet_grid(bias ~ timestep) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'bottom',
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    #save each plot
    ggsave(hourly.diff,
           filename = file.path(plot.output,
                                paste0('FluxCorrections_',
                                       gsub(' ', '', scenarios[i]),
                                       '.jpg')),
           height = 7.5, width = 8, units = 'in', device = 'jpg')
    
  }
  ###############################
  ###############################
  bias_percent_uncert <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(bias_percent_uncert) <- c('Bias', 'Point.Sources',
                                  'Date', 'Perc.Uncert')
  for(i in 1:length(bias.paths)) {
    
    bias.value <- str_split_fixed(basename(bias.paths[i]),
                                  pattern = '_', n = 2)[1]
    bias.value <- as.numeric(bias.value)
    
    LPS.value <- str_split_fixed(basename(bias.paths[i]),
                                 pattern = '_', n = 2)[2]
    if(LPS.value == 'LPS') LPS.value <- 'Large Point Sources'
    if(LPS.value == 'noLPS') LPS.value <- 'No Large Point Sources'
    
    SAM.paths <- list.files(bias.paths[i], full.names = TRUE)
    
    for(j in 1:length(SAM.paths)) {
      
      date <- str_split_fixed(basename(SAM.paths[j]),
                              pattern = '_', n = 5)[3]
      date.pos <- as.POSIXct(date, format = '%Y%m%d%H',
                             tz = 'UTC')
      attr(date.pos, 'tzone') <- tz
      Date <- strftime(date.pos, tz = tz,
                       format = '%Y-%m-%d')
      
      perc_unc_red.path <- list.files(file.path(SAM.paths[j], inversion.out.ext),
                                      pattern = 'perc_unc_red.rds',
                                      full.names = TRUE)
      perc_unc_red <- readRDS(perc_unc_red.path)
      
      add.line <- data.frame(Bias = bias.value,
                             Point.Sources = LPS.value,
                             date = Date,
                             Perc.Uncert = perc_unc_red[1,1])
      bias_percent_uncert <- rbind(bias_percent_uncert, add.line)
    }
  }
  
  bias_percent_uncert$Bias <- as.character(bias_percent_uncert$Bias)
  
  perc_unc.plot <- ggplot() +
    ggtitle('Uncertainty Reduction vs. Bias') +
    xlab('SAM Date') +
    ylab('Percent Uncertainty Reduction (%)') +
    geom_line(data = bias_percent_uncert,
              aes(x = date, y = Perc.Uncert,
                  group = Bias,
                  color = Bias)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(. ~ Point.Sources)
  ggsave(perc_unc.plot, device = 'jpg',
         filename = file.path(plot.output, 'PercentUncert.jpg'),
         height = 4.5, width = 6, units = 'in')
  
  
  
} #closes function
