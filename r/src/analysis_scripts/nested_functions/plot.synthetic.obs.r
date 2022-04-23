plot.synthetic.obs <- function(high.res.path = NULL, tz = NULL, gg.map = NULL,
                               plot.output = 'Out') {
  
  library(tidyverse)
  
  #list the available SAMs
  SAM.list <- list.files(high.res.path, full.names = TRUE)
  
  #list the available all_xco2 files
  xco2.path <- list.files(file.path(SAM.list, 'inversion/out/Analysis'),
                          full.names = TRUE)
  
  #read in the files
  xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
  xco2$time <- as.POSIXct(xco2$time, tz = 'UTC')
  
  #change to local time
  attr(xco2$time, 'tzone') <- tz

  synth.obs_plot <- ggmap(gg.map) +
    ggtitle(expression(paste('Pseudo-Enhancements (',
                             Delta, 'XCO'[2], ')'))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_point(data = xco2, shape = 23, size = 0.75,
               aes(x = lon, y = lat, fill = obs)) +
    scale_fill_gradient2(low = 'yellow', mid = 'orange', high = 'red',
                         midpoint = 2,
                       name = expression(paste('XCO'[2], ' [ppm]')),
                       breaks = c(0, 1, 2, 3, 4),
                       labels = c(paste0('\u2264', 0),
                                  '1', '2', '3',
                                  paste0('\u2265', 4)),
                       limits = c(0, 4),
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
    facet_wrap(. ~ strftime(time, tz = tz), ncol = 4)
  
  ggsave(synth.obs_plot, device = 'jpg',
         height = 8, width = 8, units = 'in',
         filename = file.path(plot.output, 'Synthetic_Obs.jpg'))
  
}
