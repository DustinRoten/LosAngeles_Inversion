coarsened.grids_xco2 <- function(coarsen.paths, inversion.out.ext = 'inversion/out',
                                 gg.map = NULL, output.path = 'Out/Coarse_Results') {
  
  for(i in 1:length(coarsen.paths)) {
    
    #obtain the aggregation factor
    agg.factor <- as.numeric(gsub('F', '', basename(coarsen.paths[i])))
    agg.factor <- paste0(agg.factor, 'km x ', agg.factor, 'km')
    
    SAM.list <- list.files(coarsen.paths[i], full.names = TRUE)
    
    xco2.path <- file.path(SAM.list, inversion.out.ext,
                           'Analysis', 'all_xco2.csv')
    xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
    
    # prior - obs
    grid.out.prior <- custom.grid(lons = xco2$lon, lats = xco2$lat,
                                  values = xco2$prior - xco2$obs,
                                  layers = agg.factor)
    
    # posterior - prior
    grid.out.posterior <- custom.grid(lons = xco2$lon, lats = xco2$lat,
                                      values =
                                        xco2$posterior - xco2$prior,
                                      layers = agg.factor)
    
    grid.out.abs.change <-
      custom.grid(lons = xco2$lon, lats = xco2$lat,
                  values = (100*abs(xco2$posterior - xco2$truth) -
                              100*abs(xco2$prior - xco2$truth))/xco2$truth,
                  layers = agg.factor)
    
    #add to the master dataframe
    if(i == 1) {
      all_grids.prior <- grid.out.prior
      all_grids.posterior <- grid.out.posterior
      all_grids.abs.change <- grid.out.abs.change
    } else if(i > 1) {
      all_grids.prior <-
        rbind(all_grids.prior, grid.out.prior)
      all_grids.posterior <-
        rbind(all_grids.posterior, grid.out.posterior)
      all_grids.abs.change <-
        rbind(all_grids.abs.change, grid.out.abs.change)
    }
    
  }
  
  #get rid of cells containing zero
  sub.all_grids.prior <- subset(all_grids.prior, values != 0)
  sub.all_grids.posterior <- subset(all_grids.posterior, values != 0)
  sub.all_grids.abs.change <- subset(all_grids.abs.change, values != 0)
  
  #plot the original prior - obs first
  prior.obs_plot <- ggmap(gg.map) +
    ggtitle(expression(paste('Prior - Pseudo-Enhancement ',
                             '[Spatially Averaged XCO'[2], ']'))) +
    labs(subtitle = 'Surface Flux Resolution' ) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = sub.all_grids.prior,
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
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ layers, ncol = 2)
  
  ggsave(prior.obs_plot, device = 'jpg', units = 'in',
         height = 6, width = 6,
         filename = file.path(output.path, 'PriorVsObs_XCO2.jpg'))
  
  
  
  #plot the correction: posterior - prior
  posterior.prior_plot <- ggmap(gg.map) +
    ggtitle(expression(paste('Posterior - Prior ',
                             '[Spatially Averaged XCO'[2], ']'))) +
    labs(subtitle = 'Surface Flux Resolution') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = sub.all_grids.posterior,
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
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ layers, ncol = 2)
  
  ggsave(posterior.prior_plot, device = 'jpg', units = 'in',
         height = 6, width = 6,
         filename = file.path(output.path, 'PosteriorVsPrior_XCO2.jpg'))
  
  
  
  #plot the correction: posterior - prior
  abs.change_plot <- ggmap(gg.map) +
    ggtitle(expression(paste('Change in Error ',
                             '[Spatially Averaged XCO'[2], ']'))) +
    labs(subtitle = 'Surface Flux Resolution') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = sub.all_grids.abs.change,
                aes(x = lons, y = lats, fill = values)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name = expression(paste(Delta,'%'['rel']))) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    coord_cartesian() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in')) +
    facet_wrap(. ~ layers, ncol = 2)
  
  ggsave(abs.change_plot, device = 'jpg', units = 'in',
         height = 6, width = 6,
         filename = file.path(output.path, 'AbsChange_XCO2.jpg'))
  
}
