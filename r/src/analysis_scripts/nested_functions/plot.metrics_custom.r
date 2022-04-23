plot.metrics_custom <- function(custom.path = NULL, inversion.out.ext = 'inversion/out',
                                include.ext = 'include', output.path = 'Out/Custom_Results',
                                local.tz = NULL) {
  
  library(ggplot2); library(ggmap); library(viridis)
  library(patchwork)
  
  #for all aggregation factors
  hourly.flux <- data.frame(matrix(NA, nrow = 0, ncol = 9))
  names(hourly.flux) <- c('SAM', 'scenario', 'timestep',
                          'Prior', 'Prior.Err', 'True', 'Posterior',
                          'Posterior.Err', 'Footprint')
  
  #for only the highest resolution
  summed.flux <- data.frame(matrix(NA, nrow = 0, ncol = 13))
  names(summed.flux) <- c('SAM', 'scenario', 'nobs', 'uq.nobs', 'mean.xco2',
                          'uq.mean.xco2', 'Prior', 'Prior.Err', 'True',
                          'Posterior', 'Posterior.Err', 'Footprint',
                          'Uncert.Red')
  for(i in 1:length(custom.path)) {
    
    #obtain the aggregation factor
    scenario <- basename(custom.path[i])

    SAM.list <- list.files(custom.path[i], full.names = TRUE)
    
    #aggregate the observation files
    observation.list <- list.files(file.path(SAM.list,
                                             inversion.out.ext,
                                             'Analysis'),
                                   pattern = 'all_xco2',
                                   full.names = TRUE)
    
    for(j in 1:length(SAM.list)) {
      
      #add observations to the aggregated dataframe      
      observations <- read.csv(observation.list[j])
      
      ### Gather some stats ###
      #number of soundings
      nobs <- nrow(observations)
      
      #number of upper quartile soundings
      uq.obs <- nrow(subset(observations,
                            obs > summary(observations$obs)[5]))
      
      #mean xco2
      mean.xco2 <- mean(observations$obs)
      
      #mean of the upper quartile xco2
      uq.mean.xco2 <- mean(subset(observations,
                                  obs > summary(observations$obs)[5])$obs)
      #########################
      
      if(i == 1 & j == 1) {
        all.observations <- data.frame(scenario, observations)
      } else {
        all.observations <- rbind(all.observations,
                                  data.frame(scenario, observations))
      }
      
      #list the prior files
      prior.files <- list.files(file.path(SAM.list[j],
                                          include.ext),
                                pattern = 'prior_emiss.nc',
                                full.names = TRUE)
      
      #list the "truth" files
      truth.files <- list.files(file.path(SAM.list[j],
                                          include.ext),
                                pattern = 'truth_emiss.nc',
                                full.names = TRUE)
      
      #list the posterior files
      posterior.files <- list.files(file.path(SAM.list[j],
                                              inversion.out.ext),
                                    pattern = 'posterior.nc',
                                    full.names = TRUE)
      
      #list the prior uncertainty files
      prior.uncert.files <- list.files(file.path(SAM.list[j],
                                                 inversion.out.ext),
                                       pattern = 'uncertainty_prior.nc',
                                       full.names = TRUE)
      
      #list the posterior uncertainty files
      posterior.uncert.files <- list.files(file.path(SAM.list[j],
                                                     inversion.out.ext),
                                           pattern = 'uncertainty.nc',
                                           full.names = TRUE)
      
      #list the error reduction files
      uncert.reduction.files <- list.files(file.path(SAM.list[j],
                                                     inversion.out.ext),
                                           pattern = 'perc_unc_red.rds',
                                           full.names = TRUE)
      
      #list the averaged footprint files
      foot.files <- list.files(file.path(SAM.list[j],
                                         include.ext),
                               pattern = 'averaged_footprint.nc',
                               full.names = TRUE)
      
      #read in the rasters
      prior <- brick(prior.files); truth <- brick(truth.files)
      posterior <- brick(posterior.files)
      prior.uncert <- sqrt(brick(prior.uncert.files)) #co-variance matrix
      posterior.uncert <- sqrt(brick(posterior.uncert.files)) #co-variance matrix
      avg.foot <- brick(foot.files)
      avg.foot <- crop(avg.foot, prior)
      
      #percent uncertainty reduction
      perc_unc_red <- readRDS(uncert.reduction.files)
      
      #get a list of timesteps from the prior file
      timesteps <- gsub('X', '', names(prior))
      
      #convert the timesteps to POSIX time
      timesteps.pos <- as.POSIXct(as.numeric(timesteps),
                                  origin = '1970-01-01',
                                  tz = 'UTC')
      
      #calculate the timestep size
      diff.time <- unique(diff(as.numeric(timesteps.pos)))
      
      #add the timestep size to the max time for the SAM time.
      SAM.time <- max(timesteps.pos) + diff.time
      
      #' Perform some quick checks to make sure the approrpriate
      #' raster files are present.
      #Test #1
      if(nlayers(prior) != nlayers(truth) |
         nlayers(prior) != nlayers(posterior))
        stop('Mismatched raster files!')
      
      #Test #2
      if(names(prior) != names(truth) ||
         names(prior) != names(posterior))
        stop('Mismatched raster files!')
      
      #Test #3
      if(timesteps.pos[1] >= timesteps.pos[length(timesteps.pos)])
        stop('Timesteps must be in chronological order!')
      
      #####################################
      ### Begin the Metric Calculations ###
      #####################################
      #prior layer converted from umol/m2/s to ktCO2
      area.raster <- raster::area(prior)
      CO2.prior <- prior*area.raster*(44.01/1e9)*diff.time
      
      #truth layer converted from umol/m2/s to ktCO2      
      area.raster <- raster::area(truth)
      CO2.truth <- truth*area.raster*(44.01/1e9)*diff.time
      
      #prior error layer converted from umol/m2/s to ktCO2
      area.raster <- raster::area(prior.uncert)
      CO2.prior.err <- prior.uncert*area.raster*(44.01/1e9)*diff.time
      
      #posterior layer converted from umol/m2/s to ktCO2
      area.raster <- raster::area(posterior)
      CO2.posterior <- posterior*area.raster*(44.01/1e9)*diff.time
      
      #posterior uncertainty
      area.raster <- raster::area(posterior.uncert)
      CO2.posterior.uncert <- posterior.uncert*area.raster*(44.01/1e9)*diff.time
      
      #determine the total area of the mean footprint still in the domain
      avg.foot[values(avg.foot) > 0] <- 1
      area.raster <- raster::area(avg.foot) #leave in km
      area.foot <- area.raster*avg.foot
      
      #no longer needed
      remove('area.raster')
      
      #add the line to the accumulation dataframe
      add.line <- data.frame(SAM = SAM.time,
                             scenario = scenario,
                             timestep = nlayers(CO2.prior):1,
                             Prior = cellStats(CO2.prior, sum),
                             Prior.Err = cellStats(CO2.prior.err, sum),
                             True = cellStats(CO2.truth, sum),
                             Posterior = cellStats(CO2.posterior, sum),
                             Posterior.Err = cellStats(CO2.posterior.uncert, sum),
                             Footprint = cellStats(area.foot, sum))
      hourly.flux <- rbind(hourly.flux, add.line)
      
      add.line2 <- data.frame(SAM = SAM.time,
                              scenario = scenario,
                              nobs = nobs,
                              uq.obs = uq.obs,
                              mean.xco2 = mean.xco2,
                              uq.mean.xco2 = uq.mean.xco2,
                              Prior = sum(add.line$Prior),
                              Prior.Err = sqrt(sum(add.line$Prior.Err^2)),
                              True = sum(add.line$True),
                              Posterior = sum(add.line$Posterior),
                              Posterior.Err = sqrt(sum(add.line$Posterior.Err^2)),
                              Footprint = sum(add.line$Footprint),
                              Uncert.Red = perc_unc_red)
      summed.flux <- rbind(summed.flux, add.line2)

    } #closes j
  } #closes the output directory list loop
  
  #change to local time
  attr(hourly.flux$SAM, 'tzone') <- local.tz
  attr(summed.flux$SAM, 'tzone') <- local.tz
  
  #add a local tz hour column
  summed.flux$hour <- hour(summed.flux$SAM) + minute(summed.flux$SAM)/60

  #change the labels
  summed.flux$scenario[summed.flux$scenario != 'Custom_F1'] <-
    'Unmodified'
  summed.flux$scenario[summed.flux$scenario == 'Custom_F1'] <-
    'Modified'
  
  ##########################################
  ### Plot Percent Uncertainty Reduction ###
  ##########################################
  Perc_Unc_Red_plot <- ggplot() +
    ggtitle('Percent Uncertainty Reduction') +
    xlab('Time of SAM [hr]') +
    ylab('Uncertainty Reduction [%]') +
    geom_line(data = summed.flux, linetype = 'dashed',
              aes(x = hour,
                  y = Uncert.Red,
                  color = scenario)) +
    geom_point(data = summed.flux, shape = 21,
               aes(x = hour,
                   y = Uncert.Red,
                   fill = scenario)) +
    scale_fill_viridis(discrete = TRUE,
                       name = 'Prior Flux: ') +
    scale_color_viridis(discrete = TRUE) +
    guides(color = FALSE, 
           fill = guide_legend(label.position = "bottom")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  ggsave(Perc_Unc_Red_plot, device = 'jpg',
         height = 4, width = 4, units = 'in',
         filename = file.path(output.path, 'Custom_Perc_Unc_Red.jpg'))

############################
############################
############################



####################
### Hourly Plots ###
####################
# First, determine blues and reds
  hourly.flux$bias <- NA
  tmp.SAM.list <- unique(hourly.flux$SAM)
  for(sam in tmp.SAM.list) {
  
    tmp.high <- subset(hourly.flux,
                       SAM == sam & timestep == 1 &
                         (abs(Posterior - True) -
                            abs(Prior - True))/Footprint > 0)
    
    tmp.low <- subset(hourly.flux,
                      SAM == sam & timestep == 1 &
                        (abs(Posterior - True) -
                           abs(Prior - True))/Footprint < 0)
    
    #apply these results to the whole dataframe
    if(nrow(tmp.high) > 0) {
      for(tmp.idx in 1:nrow(tmp.high)) {
        idx <- which(hourly.flux$SAM == tmp.high$SAM[tmp.idx] &
                       hourly.flux$scenario == tmp.high$scenario[tmp.idx])
        
        hourly.flux$bias[idx] <- 'High'
      }
    }
    
    if(nrow(tmp.low) > 0) {
      for(tmp.idx in 1:nrow(tmp.low)) {
        idx <- which(hourly.flux$SAM == tmp.low$SAM[tmp.idx] &
                       hourly.flux$scenario == tmp.low$scenario[tmp.idx])
        
        hourly.flux$bias[idx] <- 'Low'
      }
    }
  }
  ###################
  
  #change the labels
  hourly.flux$scenario[hourly.flux$scenario != 'Custom_F1'] <-
    'Unmodified'
  hourly.flux$scenario[hourly.flux$scenario == 'Custom_F1'] <-
    'Modified'
  
  ### Hourly Reduction Plot ###
  # (Bias)
  hourly.reduction.bias.plot <- ggplot() +
    ggtitle('Change in Error vs. Backward Timestep') +
    xlab('Backward Timestep [h]') +
    ylab(expression(paste(Delta, ' ', frac('Abs Error', 'Footprint Area'),
                          ' [kg km'^-2, ']'))) +
    geom_line(data = hourly.flux, linetype = 'dashed',
              aes(x = timestep,
                  y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
                  group = paste0(substr(as.character(SAM), 1, 10), '\n',
                                 substr(as.character(SAM), 12, 16)),
                  color = bias)) +
    geom_point(data = hourly.flux, shape = 21, size = 2, alpha = 0.5,
               aes(x = timestep,
                   y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
                   fill = bias),
               color = 'black') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_fill_discrete(name = 'SAM') +
    guides(color = FALSE, fill = FALSE) +
    scale_y_continuous() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(0.5, 'in')) +
    facet_wrap(. ~ scenario, ncol = 2, scales = 'free')
  
  ggsave(hourly.reduction.bias.plot, device = 'jpg',
         height = 3, width = 6, units = 'in',
         filename = file.path(output.path, 'Hourly_Corrections_Bias.jpg'))
  ####################
  
  
  
  # ### Plot the Power Fit ###
  # power.fit_df <-
  #   data.frame(matrix(NA, nrow = length(agg.factors), ncol = 6))
  # names(power.fit_df) <- c('agg.factor', 'slope', 'y.int',
  #                          'RMSE', 'R^2', 'p.value')
  # for(agg in 1:length(agg.factors)) {
  #   power.fit <- lm(data = subset(hourly.flux,
  #                                 agg.factor == agg.factors[agg]),
  #                   log10(1e6*abs(Posterior - Prior)/Footprint) ~ log10(timestep + 1)) #convert to kg
  #   power.fit_summary <- summary(power.fit)
  #   
  #   power.fit_df$agg.factor[agg] <- agg.factors[agg]
  #   power.fit_df$slope[agg] <- power.fit_summary$coefficients[2,1]
  #   power.fit_df$y.int[agg] <- power.fit_summary$coefficients[1,1]
  #   power.fit_df$RMSE[agg] <- sqrt(mean(power.fit_summary$residuals^2))
  #   power.fit_df$`R^2`[agg] <- power.fit_summary$r.squared
  #   power.fit_df$p.value[agg] <- power.fit_summary$coefficients[2,4]
  #   
  # }
  # 
  # #plot hourly corrections
  # power.fit <- ggplot() +
  #   ggtitle('Decay of Corrective Power') +
  #   labs(subtitle = expression(paste('Grid Resolutions (Log'[10], '-', 'Log'[10], ' Plot)')),
  #        caption = expression('Linearized Power Function: '*y == frac(k[0],'('*Delta*'t'+1*')'^'m'))) +
  #   xlab(expression(paste(Delta, 't + 1 [h]'))) +
  #   ylab(expression(paste(frac('|Posterior - Prior|', 'Footprint Area'),
  #                         ' [kg ', 'km'^-2, ']'))) +
  #   geom_line(data = hourly.flux, linetype = 'dashed',
  #             aes(x = log10(timestep + 1),
  #                 y = log10(abs(Posterior - Prior)/Footprint),
  #                 group = paste0(substr(as.character(SAM), 1, 10), '\n',
  #                                substr(as.character(SAM), 12, 16))),
  #                 color = 'gray') +
  #   geom_point(data = hourly.flux, shape = 21, alpha = 0.65,
  #              aes(x = log10(timestep + 1),
  #                  y = log10(abs(Posterior - Prior)/Footprint)),
  #                  fill = 'gray') +
  #   geom_text(data = power.fit_df, hjust = 0, vjust = 0,
  #             aes(x = -Inf, y = -Inf,
  #                 label = paste0(' k', '\u2080', '=', round(10^y.int, 2),
  #                                '\u338F', '/', '\u33A2', '\n',
  #                                ' m=', -1*round(slope, 2), '\n'))) +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.key.height = unit(0.5, 'in')) +
  #   facet_wrap(. ~ agg.factor, ncol = 2) +
  #   geom_smooth(data = hourly.flux, method = 'lm', se = FALSE,
  #               linetype = 'dashed', color = 'black',
  #             aes(x = log10(timestep + 1),
  #                 y = log10(abs(Posterior - Prior)/Footprint))) +
  #   scale_fill_discrete(name = 'SAM') +
  #   scale_color_discrete() +
  #   guides(color = FALSE) +
  #   theme_classic() +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   facet_wrap(. ~ agg.factor)
  # 
  # ggsave(power.fit, device = 'jpg',
  #        height = 6, width = 7, unit = 'in',
  #        filename = file.path(output.path, 'Power_Fit.jpg'))
  ####################
  ####################
  ####################
  
} #closes the function
