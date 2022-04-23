plot.metrics <- function(coarsen.paths = NULL, inversion.out.ext = 'inversion/out',
                         include.ext = 'include', output.path = 'Out/Coarse_Results',
                         local.tz = NULL) {
  
  library(ggplot2); library(ggmap); library(viridis)
  library(patchwork)
  
  #for all aggregation factors
  hourly.flux <- data.frame(matrix(NA, nrow = 0, ncol = 9))
  names(hourly.flux) <- c('SAM', 'agg.factor', 'timestep',
                          'Prior', 'Prior.Err', 'True', 'Posterior',
                          'Posterior.Err', 'Weighted.Avg')
  
  #for only the highest resolution
  summed.flux <- data.frame(matrix(NA, nrow = 0, ncol = 13))
  names(summed.flux) <- c('SAM', 'agg.factor', 'nobs', 'uq.nobs', 'mean.xco2',
                          'uq.mean.xco2', 'Prior', 'Prior.Err', 'True',
                          'Posterior', 'Posterior.Err', 'Weighted.Avg',
                          'Uncert.Red')
  
  for(i in 1:length(coarsen.paths)) {
    
    #obtain the aggregation factor
    agg.factor <- as.numeric(gsub('F', '', basename(coarsen.paths[i])))
    agg.factor <- paste0(agg.factor, 'km x ', agg.factor, 'km')
    
    SAM.list <- list.files(coarsen.paths[i], full.names = TRUE)
    
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
        all.observations <- data.frame(agg.factor, observations)
      } else {
        all.observations <- rbind(all.observations,
                                  data.frame(agg.factor, observations))
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
      Weighted.Avg <- 
        (cellStats(area.raster, sum)*(44.01/1e9)*diff.time)*
        cellStats(avg.foot*(posterior - prior), sum)/
        (cellStats(avg.foot, sum))
      
      #no longer needed
      remove('area.raster')
      
      #add the line to the accumulation dataframe
      add.line <- data.frame(SAM = SAM.time,
                             agg.factor = agg.factor,
                             timestep = nlayers(CO2.prior):1,
                             Prior = cellStats(CO2.prior, sum),
                             Prior.Err = cellStats(CO2.prior.err, sum),
                             True = cellStats(CO2.truth, sum),
                             Posterior = cellStats(CO2.posterior, sum),
                             Posterior.Err = cellStats(CO2.posterior.uncert, sum),
                             Weighted.Avg = Weighted.Avg)
      hourly.flux <- rbind(hourly.flux, add.line)
      
      add.line2 <- data.frame(SAM = SAM.time,
                              agg.factor = agg.factor,
                              nobs = nobs,
                              uq.obs = uq.obs,
                              mean.xco2 = mean.xco2,
                              uq.mean.xco2 = uq.mean.xco2,
                              Prior = sum(add.line$Prior),
                              Prior.Err = sqrt(sum(add.line$Prior.Err^2)),
                              True = sum(add.line$True),
                              Posterior = sum(add.line$Posterior),
                              Posterior.Err = sqrt(sum(add.line$Posterior.Err^2)),
                              Weighted.Avg = sum(add.line$Weighted.Avg),
                              Uncert.Red = perc_unc_red)
      summed.flux <- rbind(summed.flux, add.line2)

    } #closes j
  } #closes the output directory list loop
  
  #change to local time
  attr(hourly.flux$SAM, 'tzone') <- local.tz
  attr(summed.flux$SAM, 'tzone') <- local.tz
  
  #add a local tz hour column
  summed.flux$hour <- hour(summed.flux$SAM) + minute(summed.flux$SAM)/60
  
  #add the agg labels
  summed.flux$label <- str_split_fixed(summed.flux$agg.factor,
                                       pattern = ' ', n = 3)[,1]
  
  #Mean XCO2 and Time of Day
  fit0 <- lm(data = summed.flux, uq.mean.xco2 ~ hour)
  summary.fit0 <- summary(fit0)
  Enhancement <- ggplot() +
    ggtitle(expression(paste('(a) Upper Quartile XCO'[2], ' vs. Time of SAM'))) +
    xlab('Time of SAM [h]') +
    ylab(expression(paste('Upper Quartile XCO'[2], ' [ppm]'))) +
    geom_smooth(data = summed.flux, method = 'lm',
                se = FALSE,
                aes(x = hour,
                    y = uq.mean.xco2),
                color = 'black', linetype = 'dashed') +
    geom_point(data = summed.flux,
               aes(x = hour,
                   y = uq.mean.xco2,
                   shape = label)) +
    geom_text(aes(x = Inf, y = Inf,
                  label = paste0('\n', 'R^2=',
                                 round(summary.fit0$r.squared, 3),
                                 '\n', 'p-value=',
                                 round(summary.fit0$coefficients[2,4], 3))),
              hjust = 1, vjust = 1) +
    scale_fill_discrete(name = 'Grid Resolution: ') +
    guides(color = FALSE, 
           fill = guide_legend(label.position = "bottom")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  #Correction vs. Observations
  fit1 <- lm(data = summed.flux, Posterior - Prior ~ nobs)
  summary.fit1 <- summary(fit1)
  fit.Observations <- ggplot() +
    ggtitle('(b) Correction vs. Soundings') +
    xlab('Number of Soundings') +
    ylab(expression(paste('Posterior Flux', ' - ', 'Prior Flux ',
                          ' [ktCO'[2], ']'))) +
    geom_line(data = summed.flux,
              linetype = 'solid', alpha = 0.75,
              aes(x = nobs,
                  y = Posterior - Prior,
                  group = label,
                  color = label)) +
    geom_point(data = summed.flux,
               aes(x = nobs,
                   y = Posterior - Prior,
                   shape = label)) +
    geom_smooth(data = summed.flux,
                method = 'lm', se = FALSE,
                aes(x = nobs,
                    y = Posterior - Prior),
                color = 'black',
                linetype = 'dashed') +
    geom_text(aes(x = Inf, y = Inf,
                  label = paste0('\n', 'R^2=',
                                 round(summary.fit1$r.squared, 3),
                                 '\n', 'p-value=',
                                 round(summary.fit1$coefficients[2,4], 3))),
              hjust = 1, vjust = 1) +
    scale_fill_discrete(name = 'Grid Resolution: ') +
    guides(color = FALSE, 
           fill = guide_legend(label.position = "bottom")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  #Correction vs. Upper Quartile Mean
  fit2 <- lm(data = summed.flux, Posterior - Prior ~ uq.mean.xco2)
  summary.fit2 <- summary(fit2)
  fit.mean.XCO2 <- ggplot() +
    ggtitle(expression(paste('(c) Correction vs. Upper Quartile XCO'[2]))) +
    xlab(expression(paste('Upper Quartile XCO'[2], ' [ppm]'))) +
    ylab(expression(paste('Posterior Flux', ' - ', 'Prior Flux ',
                          '[ktCO'[2], ']'))) +
    geom_line(data = summed.flux,
              linetype = 'solid', alpha = 0.75,
              aes(x = uq.mean.xco2,
                  y = Posterior - Prior,
                  group = label,
                  color = label)) +
    geom_point(data = summed.flux,
               aes(x = uq.mean.xco2,
                   y = Posterior - Prior,
                   shape = label)) +
    geom_smooth(data = summed.flux,
                method = 'lm', se = FALSE,
                aes(x = uq.mean.xco2,
                    y = Posterior - Prior),
                color = 'black', linetype = 'dashed') +
    geom_text(aes(x = -Inf, y = Inf,
                  label = paste0('\n', ' R^2=',
                                 round(summary.fit2$r.squared, 3),
                                 '\n', ' p-value=',
                                 round(summary.fit2$coefficients[2,4], 3))),
              hjust = 0, vjust = 1) +
    scale_fill_discrete(name = 'Grid Resolution: ') +
    guides(color = FALSE, 
           fill = guide_legend(label.position = "bottom")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  #Correction vs. Time of SAM
  fit3 <- lm(data = summed.flux, Posterior - Prior ~ hour)
  summary.fit3 <- summary(fit3)
  time.of.SAM <- ggplot() +
    ggtitle('(d) Correction vs. Time of SAM') +
    xlab('Time of SAM [h]') +
    ylab(expression(paste('Posterior Flux', ' - ', 'Prior Flux ',
                          '[ktCO'[2], ']'))) +
    geom_line(data = summed.flux,
              linetype = 'solid', alpha = 0.75,
              aes(x = hour,
                  y = Posterior - Prior,
                  group = label,
                  color = label)) +
    geom_point(data = summed.flux,
               aes(x = hour,
                   y = Posterior - Prior,
                   shape = label)) +
    geom_smooth(data = summed.flux,
                method = 'lm', se = FALSE,
                aes(x = hour,
                    y = Posterior - Prior),
                color = 'black', linetype = 'dashed') +
    geom_text(aes(x = Inf, y = Inf,
                  label = paste0('\n', 'R^2=',
                                 round(summary.fit3$r.squared, 3),
                                 '\n', 'p-value=',
                                 round(summary.fit3$coefficients[2,4], 3))),
              hjust = 1, vjust = 1) +
    scale_fill_discrete(name = 'Grid Resolution: ') +
    guides(color = FALSE, 
           fill = guide_legend(label.position = "bottom")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  #Build the multi-panel plot
  combo.plot <-
    (Enhancement | fit.Observations) /
    (fit.mean.XCO2 | time.of.SAM)
  ggsave(combo.plot, filename = file.path(output.path, 'Regressions.jpg'),
         device = 'jpg', width = 8.2, height = 7.5, units = 'in')
  
  
  
  ############################
  ### Table of Summed Fits ###
  ############################
  agg.factors <- unique(summed.flux$agg.factor)
  
  #upper quartile and hour
  fit0_summed.table <-
    data.frame(matrix(NA, nrow = length(agg.factors), ncol = 5))
  names(fit0_summed.table) <-
    c('agg.factor', 'slope', 'y.int', 'R2', 'p.value')
  
  #correction vs. obs
  fit1_summed.table <-
    data.frame(matrix(NA, nrow = length(agg.factors), ncol = 5))
  names(fit1_summed.table) <-
    c('agg.factor', 'slope', 'y.int', 'R2', 'p.value')
  
  #correction vs. upper quartile
  fit2_summed.table <-
    data.frame(matrix(NA, nrow = length(agg.factors), ncol = 5))
  names(fit2_summed.table) <-
    c('agg.factor', 'slope', 'y.int', 'R2', 'p.value')
  
  #correction vs. hour
  fit3_summed.table <-
    data.frame(matrix(NA, nrow = length(agg.factors), ncol = 5))
  names(fit3_summed.table) <-
    c('agg.factor', 'slope', 'y.int', 'R2', 'p.value')
  
  for(agg in 1:length(agg.factors)) {
    
    #upper quartile vs. hour
    write.table(x = 'Upper Quartile vs. Hour',
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = FALSE, sep = '\t')
    
    fit0 <- lm(data = subset(summed.flux, agg.factor == agg.factors[agg]),
               uq.mean.xco2 ~ hour)
    summary.fit0 <- summary(fit0)
    fit0_summed.table$agg.factor[agg] <- agg.factors[agg]
    fit0_summed.table$slope[agg] <- summary.fit0$coefficients[2,1]
    fit0_summed.table$y.int[agg] <- summary.fit0$coefficients[1,1]
    fit0_summed.table$R2[agg] <- summary.fit0$r.squared
    fit0_summed.table$p.value[agg] <- summary.fit0$coefficients[2,4]
    
    write.table(x = fit0_summed.table,
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    ##########
    
    #correction vs. obs
    write.table(x = 'Correction vs. Obs',
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    
    fit1 <- lm(data = subset(summed.flux, agg.factor == agg.factors[agg]),
               Posterior - Prior ~ nobs)
    summary.fit1 <- summary(fit1)
    fit1_summed.table$agg.factor[agg] <- agg.factors[agg]
    fit1_summed.table$slope[agg] <- summary.fit1$coefficients[2,1]
    fit1_summed.table$y.int[agg] <- summary.fit1$coefficients[1,1]
    fit1_summed.table$R2[agg] <- summary.fit1$r.squared
    fit1_summed.table$p.value[agg] <- summary.fit1$coefficients[2,4]
    
    write.table(x = fit1_summed.table,
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    ##########
    
    #correction vs. upper quartile
    write.table(x = 'Correction vs. Upper Quartile',
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    
    fit2 <- lm(data = subset(summed.flux, agg.factor == agg.factors[agg]),
               Posterior - Prior ~ uq.mean.xco2)
    summary.fit2 <- summary(fit2)
    fit2_summed.table$agg.factor[agg] <- agg.factors[agg]
    fit2_summed.table$slope[agg] <- summary.fit2$coefficients[2,1]
    fit2_summed.table$y.int[agg] <- summary.fit2$coefficients[1,1]
    fit2_summed.table$R2[agg] <- summary.fit2$r.squared
    fit2_summed.table$p.value[agg] <- summary.fit2$coefficients[2,4]
    
    write.table(x = fit2_summed.table,
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    ##########
    
    #correction vs. hour
    write.table(x = 'Correction vs. Hour',
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    
    fit3 <- lm(data = subset(summed.flux, agg.factor == agg.factors[agg]),
               Posterior - Prior ~ hour)
    summary.fit3 <- summary(fit3)
    fit3_summed.table$agg.factor[agg] <- agg.factors[agg]
    fit3_summed.table$slope[agg] <- summary.fit3$coefficients[2,1]
    fit3_summed.table$y.int[agg] <- summary.fit3$coefficients[1,1]
    fit3_summed.table$R2[agg] <- summary.fit3$r.squared
    fit3_summed.table$p.value[agg] <- summary.fit3$coefficients[2,4]
    
    write.table(x = fit3_summed.table,
                file = file.path(output.path, 'Metrics_Table.tsv'),
                append = TRUE, sep = '\t')
    ##########
    
    
    
    ##########################################
    ### Plot Percent Uncertainty Reduction ###
    ##########################################
    Perc_Unc_Red_plot <- ggplot() +
      ggtitle('Percent Uncertainty Reduction') +
      xlab('Time of SAM [h]') +
      ylab('Uncertainty Reduction [%]') +
      geom_line(data = summed.flux, linetype = 'dashed',
                aes(x = hour,
                    y = Uncert.Red,
                    color = label)) +
      geom_point(data = summed.flux, shape = 21,
                 aes(x = hour,
                     y = Uncert.Red,
                     fill = label)) +
      scale_fill_viridis(discrete = TRUE,
                         name = 'Grid Resolution: ') +
      scale_color_viridis(discrete = TRUE) +
      guides(color = FALSE, 
             fill = guide_legend(label.position = "bottom")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'bottom')
    ggsave(Perc_Unc_Red_plot, device = 'jpg',
           height = 4, width = 4, units = 'in',
           filename = file.path(output.path, 'Perc_Unc_Red.jpg'))
  }
  ############################
  ############################
  ############################
  
  
  ##!! These section for Hourly Plots is no longer working !!##
  ####################
  ### Hourly Plots ###
  ####################
  # # First, determine blues and reds
  # hourly.flux$bias <- NA
  # tmp.SAM.list <- unique(hourly.flux$SAM)
  # for(sam in tmp.SAM.list) {
  # 
  #   tmp.high <- subset(hourly.flux,
  #                      SAM == sam & timestep == 1 &
  #                        (abs(Posterior - True) -
  #                           abs(Prior - True))/Footprint > 0)
  #   
  #   tmp.low <- subset(hourly.flux,
  #                     SAM == sam & timestep == 1 &
  #                       (abs(Posterior - True) -
  #                          abs(Prior - True))/Footprint < 0)
  #   
  #   #apply these results to the whole dataframe
  #   if(nrow(tmp.high) > 0) {
  #     for(tmp.idx in 1:nrow(tmp.high)) {
  #       idx <- which(hourly.flux$SAM == tmp.high$SAM[tmp.idx] &
  #                      hourly.flux$agg.factor == tmp.high$agg.factor[tmp.idx])
  #       
  #       hourly.flux$bias[idx] <- 'High'
  #     }
  #   }
  #   
  #   if(nrow(tmp.low) > 0) {
  #     for(tmp.idx in 1:nrow(tmp.low)) {
  #       idx <- which(hourly.flux$SAM == tmp.low$SAM[tmp.idx] &
  #                      hourly.flux$agg.factor == tmp.low$agg.factor[tmp.idx])
  #       
  #       hourly.flux$bias[idx] <- 'Low'
  #     }
  #   }
  # }
  ###################
  
  
  # ### Hourly Reduction Plot ###
  # # (All colors)
  # hourly.reduction.plot <- ggplot() +
  #   ggtitle('Change in Error vs. Backward Timestep') +
  #   xlab('Backward Timestep [h]') +
  #   ylab(expression(paste(Delta, ' ', frac('Abs Error', 'Footprint Area'),
  #                         ' [kg', ' km'^-1, ']'))) +
  #   geom_line(data = hourly.flux, linetype = 'dashed',
  #             aes(x = timestep,
  #                 y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
  #                 group = paste0(substr(as.character(SAM), 1, 10), '\n',
  #                                substr(as.character(SAM), 12, 16)),
  #                 color = paste0(substr(as.character(SAM), 1, 10), '\n',
  #                                substr(as.character(SAM), 12, 16)))) +
  #   geom_point(data = hourly.flux, shape = 21, size = 2,
  #              aes(x = timestep,
  #                  y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
  #                  fill = paste0(substr(as.character(SAM), 1, 10), '\n',
  #                                substr(as.character(SAM), 12, 16))),
  #              color = 'black') +
  #   geom_hline(yintercept = 0, linetype = 'dashed') +
  #   scale_fill_discrete(name = 'SAM') +
  #   guides(color = FALSE) +
  #   theme_classic() +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.key.height = unit(0.5, 'in')) +
  #   facet_wrap(. ~ agg.factor, ncol = 2)
  # 
  # ggsave(hourly.reduction.plot, device = 'jpg',
  #        height = 7, width = 8, units = 'in',
  #        filename = file.path(output.path, 'Hourly_Corrections.jpg'))
  # ####################
  
  
  
  ### Hourly Reduction Plot ###
  # (Bias)
  # hourly.reduction.bias.plot <- ggplot() +
  #   ggtitle('Change in Error vs. Backward Timestep') +
  #   xlab('Backward Timestep [h]') +
  #   ylab(expression(paste(Delta, ' ', frac('Abs Error', 'Footprint Area'),
  #                         ' [kg', ' km'^-2, ']'))) +
  #   geom_line(data = hourly.flux, linetype = 'dashed',
  #             aes(x = timestep,
  #                 y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
  #                 group = paste0(substr(as.character(SAM), 1, 10), '\n',
  #                                substr(as.character(SAM), 12, 16)),
  #                 color = bias)) +
  #   geom_point(data = hourly.flux, shape = 21, size = 2, alpha = 0.5,
  #              aes(x = timestep,
  #                  y = 1e6*(abs(Posterior - True) - abs(Prior - True))/Footprint, #convert to kg
  #                  fill = bias),
  #              color = 'black') +
  #   geom_hline(yintercept = 0, linetype = 'dashed') +
  #   scale_fill_discrete(name = 'SAM') +
  #   guides(color = FALSE, fill = FALSE) +
  #   scale_y_continuous() +
  #   theme_classic() +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.key.height = unit(0.5, 'in')) +
  #   facet_wrap(. ~ agg.factor, ncol = 2)
  # 
  # ggsave(hourly.reduction.bias.plot, device = 'jpg',
  #        height = 6, width = 7, units = 'in',
  #        filename = file.path(output.path, 'Hourly_Corrections_Bias.jpg'))
  ####################
  
  
  
  ### Plot the Power Fit ###
  power.fit_df <-
    data.frame(matrix(NA, nrow = length(agg.factors), ncol = 6))
  names(power.fit_df) <- c('agg.factor', 'slope', 'y.int',
                           'RMSE', 'R^2', 'p.value')
  for(agg in 1:length(agg.factors)) {
    power.fit <- lm(data = subset(hourly.flux,
                                  agg.factor == agg.factors[agg]),
                    log10(abs(Weighted.Avg)) ~ log10(timestep + 1)) #convert to kg
    power.fit_summary <- summary(power.fit)
    
    power.fit_df$agg.factor[agg] <- agg.factors[agg]
    power.fit_df$slope[agg] <- power.fit_summary$coefficients[2,1]
    power.fit_df$y.int[agg] <- power.fit_summary$coefficients[1,1]
    power.fit_df$RMSE[agg] <- sqrt(mean(power.fit_summary$residuals^2))
    power.fit_df$`R^2`[agg] <- power.fit_summary$r.squared
    power.fit_df$p.value[agg] <- power.fit_summary$coefficients[2,4]
    
  }
  
  #plot hourly corrections
  power.fit <- ggplot() +
    ggtitle('Decay of Corrective Power') +
    xlab(expression(paste('log'[10], '(', Delta, 't + 1) [h]'))) +
    ylab(expression(paste('Correction [kg km'^-2, ']'))) +
    geom_line(data = hourly.flux, linetype = 'dashed',
              aes(x = log10(timestep + 1),
                  y = log10(abs(Weighted.Avg)),
                  group = paste0(substr(as.character(SAM), 1, 10), '\n',
                                 substr(as.character(SAM), 12, 16))),
              color = 'gray') +
    geom_point(data = hourly.flux, shape = 21, alpha = 0.65,
               aes(x = log10(timestep + 1),
                   y = log10(abs(Weighted.Avg)))) +
    geom_text(data = power.fit_df, hjust = 0, vjust = 0,
              aes(x = -Inf, y = -Inf,
                  label = paste0(' k', '\u2080', '=', round(10^(y.int), 2),
                                 '\u338F', '/', '\u33A2', '\n',
                                 ' m=', -1*round(slope, 3), '\n'))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(0.5, 'in')) +
    facet_wrap(. ~ agg.factor, ncol = 2) +
    geom_smooth(data = hourly.flux, method = 'lm', se = FALSE,
                linetype = 'dashed', color = 'black',
                aes(x = log10(timestep + 1),
                    y = log10(abs(Weighted.Avg)))) +
    scale_fill_discrete(name = 'SAM') +
    scale_color_discrete() +
    guides(color = FALSE) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(. ~ agg.factor)
  
  ggsave(power.fit, device = 'jpg',
         height = 6, width = 7, unit = 'in',
         filename = file.path(output.path, 'Power_Fit.jpg'))
  
  
  
  ### Decay with Footprint Area ###
  new.power.fit <- ggplot() +
    ggtitle('Decay of Corrective Power') +
    labs(subtitle = expression(paste('Grid Resolutions [Log'[10], '-', 'Log'[10], ' Plot]')),
         caption = expression('Linearized Power Function: '*y == frac(k[0],'('*Delta*'t'+1*')'^'m'))) +
    xlab(expression(paste('log'[10], '(', Delta, 't + 1) [h]'))) +
    ylab(expression(paste('log'[10], '(Footprint Weighted Correction) [mtCO'[2], ']'))) +
    geom_line(data = hourly.flux, linetype = 'dashed',
              aes(x = log10(timestep + 1),
                  y = log10(abs(Weighted.Avg)),
                  group = paste0(substr(as.character(SAM), 1, 10), '\n',
                                 substr(as.character(SAM), 12, 16))),
              color = 'gray') +
    geom_point(data = hourly.flux, shape = 21, alpha = 0.65,
               aes(x = log10(timestep + 1),
                   y = log10(abs(Weighted.Avg))),
               fill = 'gray') +
    geom_text(data = power.fit_df, hjust = 0, vjust = 0,
              aes(x = -Inf, y = -Inf,
                  label = paste0(' k', '\u2080', '=', round(10^y.int, 2),
                                 ' mtCO', '\u2082', '\n',
                                 ' m=', -1*round(slope, 2), '\n'))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(0.5, 'in')) +
    facet_wrap(. ~ agg.factor, ncol = 2) +
    geom_smooth(data = hourly.flux, method = 'lm', se = FALSE,
                linetype = 'dashed', color = 'black',
                aes(x = log10(timestep + 1),
                    y = log10(abs(Weighted.Avg)))) +
    scale_fill_discrete(name = 'SAM') +
    scale_color_discrete() +
    guides(color = FALSE) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(. ~ agg.factor)
  
  ggsave(new.power.fit, device = 'jpg',
         height = 6, width = 7, unit = 'in',
         filename = file.path(output.path, 'Updated_Power_Fit.jpg'))
  
  ####################
  ####################
  ####################
  
} #closes the function
