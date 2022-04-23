diff.flux.data <- function(file_path = NULL, gg.map = NULL) {
  
  library(ggallin); library(scales); library(patchwork)
  
  #list all SAMs
  SAM.list <- list.files(file_path, full.names = TRUE)
  
  for(i in 1:length(SAM.list)) {
    
    #grab the prior emissions
    prior <- list.files(file.path(SAM.list[i], 'include'),
                        pattern = 'prior_emiss.nc',
                        full.names = TRUE)
    prior <- brick(prior)
    
    #grab the true emissions
    truth <- list.files(file.path(SAM.list[i], 'include'),
                        pattern = 'truth_emiss.nc',
                        full.names = TRUE)
    truth <- brick(truth)
    
    if(i == 1) {
      sum.prior <- prior
      sum.truth <- truth
    } else if(i > 1) {
      sum.prior <- sum.prior + prior
      sum.truth <- sum.truth + truth
    }
  } #closes the SAM.list
  
  #average the rasters
  avg.prior <- calc(sum.prior/length(SAM.list), mean)
  avg.truth <- calc(sum.truth/length(SAM.list), mean)
  remove('sum.prior'); remove('sum.truth')
  
  #calculate the difference
  diff <- avg.prior - avg.truth
  
  #change to dataframe
  p_df <- raster::as.data.frame(avg.prior, xy = TRUE)
  t_df <- raster::as.data.frame(avg.truth, xy = TRUE)
  d_df <- raster::as.data.frame(diff, xy = TRUE)
  
  p_df$Source <- 'ODIAC-VIIRS (Prior)'
  t_df$Source <- 'Vulcan 3.0 ("True")'
  
  inventories <- rbind(p_df, t_df)
  inventories <- subset(inventories, abs(layer) > 1)
  
  d_df <- subset(d_df, abs(layer) > 1)
  
  #plot the inventories
  avg.p.t <- ggmap(gg.map) +
    ggtitle('Averaged Prior and "True" Flux') +
    labs(caption = '(Abs values < 1 excluded.)') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = inventories,
                aes(x = x, y = y, fill = layer)) +
    scale_fill_viridis(name =
                         expression(paste('Flux [', mu, 'mol m'^-2, ' s'^-1, ']   ')),
                       trans = pseudolog10_trans,
                       limits = c(0, 1000),
                       breaks = c(0, 10, 100, 1000),
                       labels = c('0', '10', '100',
                                  paste0('\u2265', '1000')),
                       oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    coord_cartesian() +
    facet_wrap(. ~ Source, ncol = 2) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'))
    
  #plot the differences
  avg.diff <- ggmap(gg.map) +
    ggtitle('Difference in Flux') +
    labs(caption = '(Abs values < 1 excluded.)') +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_raster(data = d_df,
                aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                         midpoint = 0,
                         name =
                           expression(paste('Flux [', mu, 'mol m'^-2, ' s'^-1, ']   ')),
                         trans = pseudolog10_trans,
                         limits = c(-1000, 1000),
                         breaks = c(-1000, -100, -10,
                                    0, 10, 100, 1000),
                         labels = c(paste0('\u2264', '-1000'),
                                    '-100', '-10',
                                    '0', '10', '100',
                                    paste0('\u2265', '1000')),
                         oob = scales::squish) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1)) +
    coord_cartesian() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.75, 'in'))
  
  combo.plot <- avg.p.t / avg.diff
  
  ggsave(combo.plot, filename = file.path('Out', 'Flux_Plots.jpg'),
         height = 9, width = 6.5, units = 'in', device = 'jpg')
  
}