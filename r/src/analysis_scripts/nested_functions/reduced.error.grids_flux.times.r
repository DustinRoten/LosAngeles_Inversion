reduced.error.grids_flux.times <- function(reduced_err.path = NULL,
                                           inversion.out.ext = 'inversion/out',
                                           include.ext = 'include',
                                           output.path = 'Out/RedErr_Results',
                                           gg.map = NULL) {
  
  #collect data into comprehensive dataframe
  flux_df <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  names(flux_df) <- c('scenario', 'backwards.timestep',
                      'lon', 'lat', 'percent.change')
  
  #collect average spatial differneces dataframe
  flux.diff_df <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(flux.diff_df) <- c('scenario', 'lon', 'lat', 'difference')
  
  #Moran's I information
  moran.i <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  names(moran.i) <- c('scenario', 'I', 'p.value')
  
  #First, let's take a look at the XCO2 per-sounding error
  for(i in 1:length(reduced_err.path)) {
    
    scenario <- basename(reduced_err.path[i])
    
    #list the available SAMs
    SAM.list <- list.files(reduced_err.path[i], full.names = TRUE)
    
    avg.prior <-
      brick(file.path(SAM.list[1], include.ext, 'prior_emiss.nc'))
    values(avg.prior) <- 0
    avg.truth <- avg.posterior <- avg.prior
    for(j in 1:length(SAM.list)) {
      
      #read in the prior raster
      prior.path <- list.files(file.path(SAM.list[j], include.ext),
                               pattern = 'prior_emiss.nc',
                               full.names = TRUE)
      prior <- brick(prior.path)
      
      #read in the truth raster
      truth.path <- list.files(file.path(SAM.list[j], include.ext),
                               pattern = 'truth_emiss.nc',
                               full.names = TRUE)
      truth <- brick(truth.path)
      
      #read in the posterior raster
      posterior.path <- list.files(file.path(SAM.list[j], inversion.out.ext),
                                   pattern = 'posterior.nc',
                                   full.names = TRUE)
      posterior <- brick(posterior.path)

      #add to the list      
      avg.prior <- avg.prior + prior
      avg.truth <- avg.truth + truth
      avg.posterior <- avg.posterior + posterior
      
    } #close the SAM list loop
    
    avg.prior <- (1/length(SAM.list))*avg.prior
    avg.truth <- (1/length(SAM.list))*avg.truth
    avg.posterior <- (1/length(SAM.list))*avg.posterior
    
    avg.calculation <-
      100*abs(avg.posterior - avg.truth)/avg.truth -
      100*abs(avg.prior - avg.truth)/avg.truth
    
    #grab each layer of the averaged raster and add them to dataframe
    for(j in nlayers(avg.calculation):1) {
      tmp_df <- raster::as.data.frame(avg.calculation[[j]], xy = TRUE)
      tmp_df.2<- data.frame(scenario = scenario,
                            backwards.timestep =
                              (nlayers(avg.calculation) + 1) - j,
                            lon = tmp_df$x,
                            lat = tmp_df$y,
                            percent.change = tmp_df[,3])
      flux_df <- rbind(flux_df, tmp_df.2)
    } #close the individual layer loop (j)
    
    #average the difference between the prior and truth rasters
    avg.diff <- calc(avg.prior - avg.truth, mean)
    avg.diff_df <- raster::as.data.frame(avg.diff, xy = TRUE)
    
    avg.diff_df.2 <- data.frame(scenario = scenario,
                                lon = avg.diff_df$x,
                                lat = avg.diff_df$y,
                                difference = avg.diff_df$layer)
    flux.diff_df <- rbind(flux.diff_df, avg.diff_df.2)
    
    #calculating Moran's I
    dist.diff <- as.matrix(dist(cbind(avg.diff_df.2$lon,
                                      avg.diff_df.2$lat)))
    dist.diff.inv <- 1/dist.diff
    diag(dist.diff.inv) <- 0
    
    #get the I information
    m.i.info <- ape::Moran.I(avg.diff_df.2$difference, dist.diff.inv,
                             scaled = TRUE, alternative = 'greater')
    
    #add to the Moran's I dataframe
    moran.i <- rbind(moran.i, data.frame(scenario = scenario,
                                         I = m.i.info$observed,
                                         p.value = m.i.info$p.value))
    
  } #close the aggregation factor loop (i)
  
  
  
  ###############################
  ### Plot Averaged Timesteps ###
  ###############################
  #filter out small changes and NA cells
  sub.flux_df <- subset(flux_df,
                        abs(percent.change) > 0.01 &
                          !is.na(percent.change))
  
  #skip some timesteps so a larger timeframe is included
  selected.times <- sub.flux_df$backwards.timestep %in%
    seq(1, max(sub.flux_df$backwards.timestep), 4)
  sub.flux_df <- sub.flux_df[selected.times,]
  
  #rename the scenarios
  sub.flux_df$scenario[sub.flux_df$scenario != 'Reduced_Err'] <-
    '100% Error'
  sub.flux_df$scenario[sub.flux_df$scenario == 'Reduced_Err'] <-
    '50% Error'
  
  #generate the plot
  abs.change_plot <- ggmap(gg.map) +
    ggtitle('Spatial Changes in Error') +
    labs(subtitle = 'Backward Timesteps and Atmospheric Transport') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = sub.flux_df,
                aes(x = lon, y = lat, fill = percent.change)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         trans = ggallin::pseudolog10_trans,
                         name = expression(paste(Delta,'%'[rel])),
                         breaks = c(-100, -10, 0, 10, 100),
                         labels = c(paste0('\u2264', '-100'),
                                    '-10', '0', '10',
                                    paste0('\u2265', '100')),
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
    coord_cartesian() +
    facet_grid(scenario ~ backwards.timestep)
  
  #save the plot
  ggsave(abs.change_plot, device = 'jpg', units = 'in',
         height = 5, width = 8,
         filename = file.path(output.path, 'AbsDiff_Change_RedErr.jpg'))
  ###############################
  ###############################
  ###############################
   
}
