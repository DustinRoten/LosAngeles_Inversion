aggregate.native.data <- function(native.path = NULL, gg.map = NULL, site, tz = local.tz,
                                  inversion.out.ext = 'inversion/out', include.ext = 'include',
                                  obs.error = 0, plot.output = 'Out/Native_Results',
                                  analysis.ext, api.key = NULL, agg_line = NULL,
                                  reg.caption = NULL) {
  #plot xco2 values
  plot.xco2(site, tz = local.tz, output.path = native.path,
            analysis.ext, obs.error, api.key,
            reg.caption = agg_line, plot.output = plot.output,
            gg.map = gg.map)
  
  #plot all of the fluxes available
  averaged.plot.flux(site, tz = local.tz, output.path = native.path,
                     prior.truth.ext, posterior.ext, api.key,
                     plot.output = plot.output, p.caption = agg_line,
                     gg.map = gg.map)
  
  #plot an hourly breakdown of the fluxes
  hourly.plot.flux(site, tz = local.tz, output.path = native.path,
                   prior.truth.ext, posterior.ext, api.key,
                   plot.output = plot.output, p.caption = agg_line,
                   gg.map = gg.map)
  
  #plot various metrics
  plot.metrics(site, tz = local.tz, output.path = native.path,
               prior.truth.ext, posterior.ext,
               plot.output = plot.output, p.caption = agg_line)
  
  #plot the footprint strength and footprint*difference
  plot.footprints(site, tz = local.tz, output.path = native.path,
                  prior.truth.ext, posterior.ext,
                  plot.output = plot.output, gg.map = gg.map)
  
}