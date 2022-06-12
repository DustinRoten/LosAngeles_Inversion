oco3.nc_to_dataframe <- function(OCO3.file, QF = T) {
  
  #load required library
  require(ncdf4)
  
  #open the nc file
  nc <- nc_open(OCO3.file)
  
  #get the OCO3 dates
  oco3.times <- ncvar_get(nc, 'time')
  oco3.times <- as.POSIXct(oco3.times, origin = '1970-01-01',
                           tz = 'UTC')
  
  #get the OCO3 longitudes
  oco3.lons <- ncvar_get(nc, 'longitude')
  
  #get the OCO3 latitudes
  oco3.lats <- ncvar_get(nc, 'latitude')
  
  #get the OCO3 quality flags
  oco3.qf <- ncvar_get(nc, 'xco2_quality_flag')
  
  #get the OCO3 xco2
  oco3.xco2 <- ncvar_get(nc, 'xco2')
  
  #get the OCO3 xco2 uncertainty
  oco3.xco2.uncert <- ncvar_get(nc, 'xco2_uncertainty')
  
  #close the *.nc file
  nc_close(nc)
  
  #construct the OCO-3 dataframe
  oco3_df <- data.frame(date.time = oco3.times,
                        lon = oco3.lons,
                        lat = oco3.lats,
                        qf = oco3.qf,
                        xco2 = oco3.xco2,
                        xco2.uncert = oco3.xco2.uncert)
  
  #filter qf flags as specified
  if(QF) oco3_df <- subset(oco3_df, qf == 0)
  
  return(oco3_df)
  
}