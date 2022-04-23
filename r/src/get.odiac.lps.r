get.odiac.lps <- function(CarMA = NULL, ODIAC = NULL) {
  
  # require the raster package
  require(raster)
  
  # determine the spatial domain necessary
  lims <- extent(ODIAC)
  
  # Grab the large points sources within the domain
  sub.CarMA <- subset(CarMA,
                      (longitude >= lims[1] & longitude <= lims[2]) &
                        (latitude >= lims[3] & latitude <= lims[4]))
  
  # Locations are all that are needed.
  lps.locations <- data.frame(sub.CarMA$longitude, sub.CarMA$latitude)
  names(lps.locations) <- c('lon', 'lat')
  
  # Use the rasterize functions to move these points to a grid
  # The relevant ODIAC subdomain will serve as the resolution guide
  lps.grid <- rasterize(x = lps.locations, y = ODIAC,
                        field = rep(1, nrow(lps.locations)))
  lps.grid[is.na(lps.grid)] <- 0
  
  # Multiply to obtain only the large point sources.
  ODIAC.lps <- ODIAC*lps.grid
  
  return(ODIAC.lps)
  
}