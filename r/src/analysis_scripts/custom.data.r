custom.data <- function(custom.path = NULL, inversion.out.ext = 'inversion/out',
                        include.ext = 'include', gg.map.zoom = NULL, gg.map.large = NULL,
                        obs.error = 0, output.path = 'Out/Custom_Results', local.tz = NULL) {
  
  #create the appropriate directory
  if(!dir.exists(output.path)) dir.create(output.path)
  
  #deal with the xco2 first
  #output correction vs. input difference
  corr.vs.diff_custom(custom.path, inversion.out.ext, obs.error, output.path)
  
  #investigate the total flux
  custom.grids_flux.totals(custom.path)
  
  #investigate the gridded flux
  custom.grids_flux.times(custom.path, gg.map = gg.map.large,
                          gg.map.zoom = gg.map.zoom)
  
  #prepare metrics for the custom gridding scheme
  plot.metrics_custom(custom.path = custom.path,
                      local.tz = local.tz)
  
}