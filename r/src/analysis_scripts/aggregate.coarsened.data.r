aggregate.coarsened.data <- function(coarsen.paths = NULL, inversion.out.ext = 'inversion/out',
                                     include.ext = 'include', gg.map.zoom = NULL, gg.map.large = NULL,
                                     obs.error = 0, output.path = 'Out/Coarse_Results', tz = 'UTC') {
  
  #output correction vs. input difference
  corr.vs.diff_coarsen(coarsen.paths, inversion.out.ext, obs.error, output.path)
  
  #output the comparisons of total emissions
  coarsened.grids_flux.totals(coarsen.paths, inversion.out.ext,
                              include.ext, obs.error, output.path)
  
  #look at averaged timesteps in flux space
  #also plot averaged spatial data
  coarsened.grids_flux.times(coarsen.paths, inversion.out.ext,
                             include.ext, output.path, gg.map = gg.map.zoom)
  
  #investigate the aggregated SAMs
  coarsened.grids_xco2(coarsen.paths, gg.map = gg.map.large)
  
  #plot metrics and regressions for this data
  plot.metrics(coarsen.paths, local.tz = tz)
  
}