corr.vs.diff_err <- function(reduced_err.path = NULL, inversion.out.ext = 'inversion/out',
                             obs.error = 0, output.path = 'Out/RedErr_Results') {
  
  #create the output directory if it doesn't exist
  if(!dir.exists(output.path)) dir.create(output.path)
  
  library(ggallin); library(Metrics); library(patchwork)
  
  #First, let's take a look at the XCO2 error
  for(i in 1:length(reduced_err.path)) {
    
    #obtain the aggregation factor
    scenario <- basename(reduced_err.path[i])
    
    #list the available SAMs
    SAM.list <- list.files(reduced_err.path[i], full.names = TRUE)
    xco2.path <- list.files(file.path(SAM.list, inversion.out.ext,
                                      'Analysis'), full.names = TRUE)
    
    #list the xco2.files
    all.xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
    all.xco2$scenario <- scenario
    
    if(i == 1)
      all.scenarios_xco2 <- all.xco2
    if(i > 1)
      all.scenarios_xco2 <- rbind(all.scenarios_xco2, all.xco2)
  }
  
  disp.values <- subset(all.scenarios_xco2, obs > obs.error & obs <= 5)
  disp.values <- subset(all.scenarios_xco2, obs <= 5)
  omitted.data <-
    100*(1 - nrow(subset(all.scenarios_xco2, obs <= 5))/nrow(all.scenarios_xco2))
  
  lm.data <- data.frame(matrix(NA, nrow = 0, ncol = 6))
  names(lm.data) <- c('scenario', 'R2', 'slope',
                      'y.int', 'y.sign', 'RMSE')
  for(i in 1:length(unique(all.scenarios_xco2$scenario))) {
    
    tmp <- subset(all.scenarios_xco2,
                  scenario == unique(all.scenarios_xco2$scenario)[i])
    
    lm.sum <- summary(lm(data = tmp,
                         formula = abs(posterior - prior) ~
                           abs(obs - prior)))
    RMSE <- rmse(abs(tmp$obs - tmp$prior),
                 abs(tmp$posterior - tmp$prior))
    
    add.line <- data.frame(scenario =
                             unique(all.scenarios_xco2$scenario)[i],
                           R2 = lm.sum$r.squared,
                           slope = lm.sum$coefficients[2,1],
                           y.int = lm.sum$coefficients[1,1],
                           y.sign = NA,
                           RMSE = RMSE)
    lm.data <- rbind(lm.data, add.line)
  }
  
  lm.data$y.sign[lm.data$y.int < 0] <- '-'
  lm.data$y.sign[lm.data$y.int >= 0] <- '+'
  
  disp.values$scenario[disp.values$scenario != 'Reduced_Err'] <-
    '100% Error'
  disp.values$scenario[disp.values$scenario == 'Reduced_Err'] <-
    '50% Error'
  
  lm.data$scenario[lm.data$scenario != 'Reduced_Err'] <-
    '100% Error'
  lm.data$scenario[lm.data$scenario == 'Reduced_Err'] <-
    '50% Error'
  
  Corr.Vs.Diff.plot <- ggplot() +
    ggtitle(paste0('Changes in Corrective Power from Reduction', '\n',
                   'in Atmospheric Transport Error')) +
    xlab(expression('|'*italic(z)-bold(H)*italic(s[p])*'| [ppm]')) +
    ylab(expression('|'*bold(H)*(hat(italic(s)) - italic(s[p]))*'| [ppm]')) +
    geom_point(data = disp.values, shape = 1, alpha = 0.5,
               color = 'black',
               aes(x = (abs(obs - prior)),
                   y = (abs(posterior - prior)))) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    geom_smooth(data = disp.values, se = FALSE,
                color = 'gold', linetype = 'dashed',
                aes(x = (abs(obs - prior)),
                    y = (abs(posterior - prior)))) +
    geom_abline(data = lm.data, linetype = 'solid',
                color = 'blue',
                aes(slope = slope, intercept = y.int)) +
    geom_text(data = lm.data, hjust = 0, vjust = 1,
              size = 3,
              aes(x = -Inf, y = Inf,
                  label = paste0('\n',
                                 ' y=', round(slope, 2), 'x',
                                 y.sign, round(abs(y.int), 2), '\n',
                                 ' R^2=', round(R2, 2)))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(. ~ scenario, ncol = 2)
  
  
  
  #############################
  ### Averaged Spatial XCO2 ###
  #############################
  gridded.SAMs <- custom.grid(lons = all.scenarios_xco2$lon,
                              lats = all.scenarios_xco2$lat,
                              values =
                                abs(all.scenarios_xco2$posterior -
                                      all.scenarios_xco2$obs) -
                                abs(all.scenarios_xco2$prior -
                                      all.scenarios_xco2$obs),
                              layers = all.scenarios_xco2$scenario)
  gridded.SAMs <- subset(gridded.SAMs, values != 0)
  
  gridded.SAMs$layers[gridded.SAMs$layers != 'Reduced_Err'] <-
    '100% Error'
  gridded.SAMs$layers[gridded.SAMs$layers == 'Reduced_Err'] <-
    '50% Error'
  
  xco2.error.reduction.plot <- ggmap(gg.map.large) +
    labs(subtitle = expression('|Posterior - Obs| - |Prior - Obs|')) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = gridded.SAMs,
                aes(x = lons, y = lats, fill = values)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name = expression(paste(Delta, 'XCO'[2], ' [ppm]'))) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    coord_cartesian() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(. ~ layers, ncol = 2)
  
  combo.plot <- Corr.Vs.Diff.plot/xco2.error.reduction.plot
  ggsave(combo.plot, device = 'jpg',
         height = 7, width = 6, units = 'in',
         filename = file.path(output.path, 'XCO2_ErrRed.jpg'))
}
