#' A function that will drop the resolution of a raster brick using the
#' `aggregate` function from the raster package. This function is
#' slower than direct application of the `aggregate` function because
#' it iteratively drops the resolution of each layer in the raster brick.
#' This keeps the process from writing temporary files that eat up memory.

drop.resolution_inRAM <- function(raster.path = NULL, factor = NULL,
                                  is.footprint = FALSE) {
  
  library(raster); library(ncdf4)
  
  #read the file in again as a NetCDF file
  nc <- nc_open(raster.path, write = TRUE)
  
  if(length(names(nc$var)) > 1)
    stop('Cannot operate on files with more than one variable!')
  
  #' The raster's current characteristics must first be obtained.
  #' This scheme has been set up to be as general as possible. This
  #' will allow several different types of rasters to be processed.
  
  #Dimension 1
  dim.1.name <- names(nc$dim)[1]
  dim.1.units <- ncatt_get(nc, names(nc$dim)[1])$units
  dim.1.longname <- ncatt_get(nc, names(nc$dim)[1])$long_name
  
  #Dimension 2
  dim.2.name <- names(nc$dim)[2]
  dim.2.units <- ncatt_get(nc, names(nc$dim)[2])$units
  dim.2.longname <- ncatt_get(nc, names(nc$dim)[2])$long_name
  
  #Dimension 3
  dim.3.name <- names(nc$dim)[3]
  dim.3.vals <- ncvar_get(nc, names(nc$dim)[3])
  dim.3.units <- ncatt_get(nc, names(nc$dim)[3])$units
  dim.3.longname <- ncatt_get(nc, names(nc$dim)[3])$long_name
  
  #Variable 1
  var.1.name <- names(nc$var)[1]
  var.1.units <- ncatt_get(nc, names(nc$var))$units
  var.1.longname <- ncatt_get(nc, names(nc$var))$long_name
  
  #close the nc file
  nc_close(nc)
  
  #read in the raster
  raster.brick <- brick(raster.path)
  
  #save the original raster's origin
  raster.origin <- raster::origin(raster.brick)
  raster.crs <- projection(raster.brick)
  
  layer.names <- NULL
  for(i in 1:nlayers(raster.brick)) {

    #reducee the resolution of each layer and stack them.
    tmp.layer <- raster.brick[[i]]; layer.names[i] <- names(tmp.layer)
    
    #remove the 'per area' component
    area.tmp <- raster::area(tmp.layer)
    tmp.layer <- tmp.layer*area.tmp

    #aggregate the contributions and reapply the 'per area'.
    agg.tmp.layer <- aggregate(tmp.layer, fact = factor,
                               fun = sum)
    agg.tmp.layer <- agg.tmp.layer/raster::area(agg.tmp.layer)
    agg.tmp.layer[is.na(agg.tmp.layer)] <- 0
    
    #form the new raster brick here
    if(i == 1) new.brick <- agg.tmp.layer
    if(i > 1) new.brick <- addLayer(new.brick, agg.tmp.layer)

  }; names(new.brick) <- layer.names
  raster::origin(new.brick) <- raster.origin
  
  if(is.footprint) new.brick <- (factor^2)*new.brick
  
  #overwrite the old raster with the new one
  writeRaster(new.brick, filename = raster.path, overwrite = TRUE)
  remove('new.brick')
  
  #open the newly saved raster brick and re-label its attributes
  nc <- nc_open(raster.path, write = TRUE)
  
  #New Dimension 1 information
  dim.1 <- ncdim_def(name = dim.1.name,
                     vals = ncvar_get(nc, names(nc$dim)[1]),
                     units = dim.1.units,
                     longname = dim.1.longname)
  ndim.1 <- length(dim.1$vals)
  
  #New Dimension 2 information
  dim.2 <- ncdim_def(name = dim.2.name,
                     vals = ncvar_get(nc, names(nc$dim)[2]),
                     units = dim.2.units,
                     longname = dim.2.longname)
  ndim.2 <- length(dim.2$vals)
  
  #New Dimension 3 information
  dim.3 <- ncdim_def(name = dim.3.name,
                     vals = dim.3.vals,
                     units = dim.3.units,
                     longname = dim.3.longname)
  ndim.3 <- length(dim.3$vals)
  
  #Variable 1
  var.1 <- ncvar_def(name = var.1.name,
                     list(dim.1, dim.2, dim.3),
                     units = var.1.units,
                     longname = var.1.longname)
  nc.vars <- list(var.1)
  
  #get the raster values
  #first name is likely "crs" by default
  nc.values <- ncvar_get(nc, names(nc$var)[2])
  
  #create the new nc file and save the variables/attributes
  nc_filename <- raster.path
  nc_new <- nc_create(nc_filename, nc.vars)
  ncvar_put(nc_new, var.1, nc.values,
            count = c(ndim.1, ndim.2, ndim.3))
  nc_close(nc_new)
  
}
