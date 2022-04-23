reduced.error.grids_flux.totals <- function(reduced_err.path = NULL,
                                            inversion.out.ext = 'inversion/out',
                                            include.ext = 'include',
                                            output.path = 'Out/RedErr_Results') {
  
  library(scales)
  
  #First, let's take a look at the XCO2 error
  total_df <- data.frame(matrix(NA, nrow = 0, ncol = 7))
  names(total_df) <- c('scenario', 'SAM', 'CO2', 'CO2.diff',
                       'CO2.err', 'source', 'over.est')
  for(i in 1:length(reduced_err.path)) {
    
    #obtain the aggregation factor
    scenario <- basename(reduced_err.path[i])
    
    #list the available SAMs
    SAM.list <- list.files(reduced_err.path[i], full.names = TRUE)
    
    #prior paths
    prior.path <- list.files(file.path(SAM.list, include.ext),
                             pattern = 'prior_emiss.nc',
                             full.names = TRUE)
    
    #prior error paths
    prior.err.path <- list.files(file.path(SAM.list, inversion.out.ext),
                                 pattern = 'uncertainty_prior.nc',
                                 full.names = TRUE)
    
    #truth paths
    truth.path <- list.files(file.path(SAM.list, include.ext),
                             pattern = 'truth_emiss.nc',
                             full.names = TRUE)
    
    #posterior paths
    posterior.path <- list.files(file.path(SAM.list, inversion.out.ext),
                                 pattern = 'posterior.nc',
                                 full.names = TRUE)
    
    #posterior error paths
    posterior.err.path <- list.files(file.path(SAM.list,
                                               inversion.out.ext),
                                     pattern = 'uncertainty.nc',
                                     full.names = TRUE)
    
    prior.tot <- prior.err.tot <- truth.tot <- NULL
    posterior.tot <- posterior.err.tot <- SAM <- NULL
    for(j in 1:length(prior.path)) {
      
      SAM[j] <- j
      
      ###################
      ### Prior Total ###
      ###################
      prior.brick <- brick(prior.path[j])
      prior.area <- raster::area(prior.brick)
      
      #lots of nested functions to get the timestep
      #(this is a pretty simple calculation, just dense)
      timestep <-
        unique(diff(as.numeric(gsub('X', '', names(prior.brick)))))
      
      prior.total <- prior.brick*prior.area*timestep
      
      #|umol to mol|mol to g|g to kg|kg to mT|
      prior.total <- prior.total*(1/1e9)*(44.01/1)
      prior.tot[j] <- sum(cellStats(prior.total, sum))
      ###################
      ###################
      
      
      
      #########################
      ### Prior Uncertainty ###
      #########################
      prior.err.brick <- sqrt(brick(prior.err.path[j]))
      prior.area <- raster::area(prior.err.brick)
      
      #lots of nested functions to get the timestep
      #(this is a pretty simple calculation, just dense)
      timestep <-
        unique(diff(as.numeric(gsub('X', '', names(prior.err.brick)))))
      
      prior.err.total <- prior.err.brick*prior.area*timestep
      
      #|umol to mol|mol to g|g to kg|kg to mT|
      prior.err.total <-
        prior.err.total*(1/1e9)*(44.01/1)
      prior.err.tot[j] <- sum(cellStats(prior.err.total, sum))
      #########################
      #########################
      
      
      
      ###################
      ### Truth Total ###
      ###################
      truth.brick <- brick(truth.path[j])
      truth.area <- raster::area(truth.brick)
      
      #lots of nested functions to get the timestep
      #(this is a pretty simple calculation, just dense)
      timestep <-
        unique(diff(as.numeric(gsub('X', '', names(truth.brick)))))
      
      truth.total <- truth.brick*truth.area*timestep
      
      #|umol to mol|mol to g|g to kg|kg to mT|
      truth.total <- truth.total*(1/1e9)*(44.01/1)
      truth.tot[j] <- sum(cellStats(truth.total, sum))
      ###################
      ###################
      
      
      
      #######################
      ### Posterior Total ###
      #######################
      posterior.brick <- brick(posterior.path[j])
      posterior.area <- raster::area(posterior.brick)
      
      #lots of nested functions to get the timestep
      #(this is a pretty simple calculation, just dense)
      timestep <-
        unique(diff(as.numeric(gsub('X', '', names(posterior.brick)))))
      
      posterior.total <- posterior.brick*posterior.area*timestep
      
      #|umol to mol|mol to g|g to kg|kg to mT|
      posterior.total <-
        posterior.total*(1/1e9)*(44.01/1)
      posterior.tot[j] <- sum(cellStats(posterior.total, sum))
      #######################
      #######################
      
      
      
      #############################
      ### Posterior Uncertainty ###
      #############################
      posterior.err.brick <- sqrt(brick(posterior.err.path[j]))
      posterior.err.area <- raster::area(posterior.err.brick)
      
      #lots of nested functions to get the timestep
      #(this is a pretty simple calculation, just dense)
      timestep <-
        unique(diff(as.numeric(gsub('X', '',
                                    names(posterior.err.brick)))))
      
      posterior.err.total <-
        posterior.err.brick*posterior.err.area*timestep
      
      #|umol to mol|mol to g|g to kg|kg to mT|
      posterior.err.total <-
        posterior.err.total*(1/1e9)*(44.01/1)
      posterior.err.tot[j] <- sum(cellStats(posterior.err.total, sum))
      #############################
      #############################
      
    } #close j loop (SAM)
    
    #add the prior data
    add.prior <- data.frame(scenario = scenario,
                            SAM = SAM,
                            CO2 = prior.tot,
                            CO2.diff = prior.tot - truth.tot,
                            CO2.err = prior.err.tot,
                            source = 'Prior',
                            over.est = FALSE)
    total_df <- rbind(total_df, add.prior)
    
    #add the truth data
    add.truth <- data.frame(scenario = scenario,
                            SAM = SAM,
                            CO2 = truth.tot,
                            CO2.diff = NA,
                            CO2.err = NA,
                            source = 'Truth',
                            over.est = FALSE)
    total_df <- rbind(total_df, add.truth)
    
    #determine if the posterior is over estimated or not
    tmp.over.est <-
      abs(posterior.tot - truth.tot) - abs(prior.tot - truth.tot)
    tmp.over.est <- tmp.over.est > 0
    
    #add the posterior data
    add.posterior <- data.frame(scenario = scenario,
                                SAM = SAM,
                                CO2 = posterior.tot,
                                CO2.diff = posterior.tot - truth.tot,
                                CO2.err = posterior.err.tot,
                                source = paste0('Posterior_',
                                                scenario),
                                over.est = tmp.over.est)
    total_df <- rbind(total_df, add.posterior)
    
  }; gc() #close i loop (aggregation factor)
  
  total_df$source[total_df$source == 'Posterior_F1'] <-
    '100% Error'
  total_df$source[total_df$source == 'Posterior_Reduced_Err'] <-
    '50% Error'
  
  total_df$source <- factor(total_df$source,
                            levels = c('Truth',
                                       'Prior',
                                       '100% Error',
                                       '50% Error'))
  
  diff_df <- subset(total_df, source != 'Truth')
  
  total.diffs_RedErr <- ggplot(data = diff_df,
                               aes(x = as.factor(SAM),
                                   y = CO2.diff,
                                   fill = source)) +
    ggtitle(paste0('Prior/Posterior Differences from True Emissions', '\n',
                   'from Reductions in Atmospheric Transport Error')) +
    geom_bar(stat = 'identity',
             position = 'dodge') +
    xlab('SAM') +
    ylab(expression(paste('Total Emissions [mtCO'[2], ']'))) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_bar(data = diff_df,
             stat = 'identity',
             aes(color = over.est),
             position = position_dodge()) +
    geom_errorbar(aes(ymin = CO2.diff - CO2.err,
                      ymax = CO2.diff + CO2.err),
                  position = position_dodge(0.9),
                  width = 0.5) +
    scale_fill_manual('Source:',
                      values = c('Prior' = 'yellow',
                                 '100% Error' = 'lightblue',
                                 '50% Error' = 'lightgreen')) +
    scale_color_manual('Over Estimate:',
                       values = c('TRUE' = 'red',
                                  'FALSE' = 'black')) +
    guides(fill = guide_legend(title = 'Source:')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  ggsave(total.diffs_RedErr, device = 'jpg',
         height = 4, width = 7, units = 'in',
         filename = file.path(output.path, 'TotalCO2.RedErr_diffs.jpg'))
  
} #close the function