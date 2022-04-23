reduced.error.data <- function(reduced_err.path = NULL, inversion.out.ext = 'inversion/out',
                               include.ext = 'include', gg.map.zoom = NULL, gg.map.large = NULL,
                               obs.error = 0, output.path = 'Out/RedErr_Results') {
  
  #deal with the xco2 first
  #output correction vs. input difference
  corr.vs.diff_err(reduced_err.path, inversion.out.ext, obs.error, output.path)
  
  #investigate the total flux
  reduced.error.grids_flux.totals(reduced_err.path)
  
  #investigate the gridded flux
  reduced.error.grids_flux.times(reduced_err.path, gg.map = gg.map.large)
  
}