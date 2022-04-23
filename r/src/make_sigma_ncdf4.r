# compile all the netcdf sector subfiles into 1 netcdf file for gridded uncertainty 
# calculation. This script computes uncertainty by taking differences between
# ODIAC and HESTIA.
#
# Written by LTK, modified by DVM 6/04/2019
# Modifed (again) by DDR 5/25/2021
make_sigma_ncdf4 <- function(prior_emiss = NULL, truth_emiss = NULL,
                             sigma.output = NULL) {

  # load the necessary packages
  library(ncdf4)
  library(lubridate)
  
  # ~~~~~~~~~~~~~~~~~~~~~~ Prior ~~~~~~~~~~~~~~~~~~~~~~~~#
  #' The `prior_emiss` variable is transferred to the
  #' `prior_nc_file` variable to be consistent with the
  #' original inversion code. (D. Roten)
  prior_nc_file <- prior_emiss
  nc_prior <- nc_open(prior_nc_file)
  lon_prior <- as.numeric(ncvar_get(nc_prior, "lon"))
  lat_prior <- as.numeric(ncvar_get(nc_prior, "lat"))
  time_prior <- ncvar_get(nc_prior, "time")
  emiss_prior <- ncvar_get(nc_prior, "emiss")
  nc_close(nc_prior)
  
  # set up the time zone
  class(time_prior) <- c("POSIXt", "POSIXct")
  attributes(time_prior)$tzone <- "UTC"
  
  # construct the grid
  nlon_prior <- length(lon_prior)
  nlat_prior <- length(lat_prior)
  lonlat_prior <- expand.grid(lon_prior, lat_prior)
  
  #averaging timescale for emissions
  #this value is guided by the timesteps in the prior
  #(Lots of nested functions but its self-explanatory)
  t_avg <- as.numeric(round(min(abs(diff(time_prior))), 0))
  
  #the hour at which to start averaging (UTC hour)
  #this value is guided by the timesteps in the prior
  utc_start_hr <- hour(min(time_prior))
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~ Truth ~~~~~~~~~~~~~~~~~~~~~~~~#
  #' In order to keep consistency with previous code,
  #' "Outer" emissions have been replaced with truth
  #' emissions yet the code retains the "Outer" designation
  #' throughout.
  outer_nc_file <- truth_emiss
  nc_outer <- nc_open(outer_nc_file)
  lon_outer <- as.numeric(ncvar_get(nc_outer, "lon"))
  lat_outer <- as.numeric(ncvar_get(nc_outer, "lat"))
  time_outer <- ncvar_get(nc_outer, "time")
  emiss_outer <- ncvar_get(nc_outer, "emiss")
  nc_close(nc_outer)
  
  class(time_outer) <- c("POSIXt", "POSIXct")
  attributes(time_outer)$tzone <- "UTC"
  
  nlon_outer <- length(lon_outer)
  nlat_outer <- length(lat_outer)
  lonlat_outer_emiss <- expand.grid(lon_outer, lat_outer)
  
  
  # it can be assumed that the prior grid falls within the outer grid
  ibig <- array(1:nrow(lonlat_outer_emiss), dim = c(nlon_outer, nlat_outer))
  sub_lon <- unique(lonlat_prior[, 1])
  sub_lat <- unique(lonlat_prior[, 2])
  imin_lat <- which(lat_outer == min(sub_lat))
  imax_lat <- which(lat_outer == max(sub_lat))
  imin_lon <- which(lon_outer == min(sub_lon))
  imax_lon <- which(lon_outer == max(sub_lon))
  
  iouter_mat <- ibig[imin_lon:imax_lon, imin_lat:imax_lat]  #indexed as [lon, lat]
  iInner <- as.vector(iouter_mat)
  
  #each timestamp will be a bin here
  time_bins <- time_outer
  
  #this is the number of timesteps covered by the prior emissions
  nbins <- length(time_bins)
  
  #Average the emissions to obtain mean vals for each time bin
  outer_binned <- array(0, dim = c(nlon_outer, nlat_outer, nbins))
  prior_binned <- array(0, dim = c(nlon_prior, nlat_prior, nbins))
  sigma_final <- array(0, dim = c(nlon_prior, nlat_prior, nbins))
  
  #find the absolute difference between the two emission inventories
  for (ii in 1:nbins) {
    
    #absolute difference between the two inventories
    diff_vec <- abs(as.vector(emiss_outer[,,ii]) - as.vector(emiss_prior[,,ii]))
    sigma_final[,,ii] <- t(matrix(diff_vec, nrow = nlat_prior, byrow = T))
    
  }
  
  sigma_times <- time_bins
  ntimes_sig <- length(sigma_times)
  
  time_dim <- ncdim_def("time", "seconds_since_1970_01_01", as.numeric(sigma_times),
                        longname = "seconds since R epoch: 1970-01-01 00:00:00")
  lat_dim <- ncdim_def("lat", "degrees_north", lat_prior,
                       longname = "latitude (center of cell)")
  lon_dim <- ncdim_def("lon", "degrees_east", lon_prior,
                       longname = "longitude (center of cell)")
  uncert_var <- ncvar_def("uncertainty", "umol m-2 s-1",
                          list(lon_dim, lat_dim, time_dim),longname = "prior uncertainty")
  
  sigma_vars <- list(uncert_var)
  
  nc_filename <- sigma.output
  nc_sigma <- nc_create(nc_filename, sigma_vars)
  ncvar_put(nc_sigma, uncert_var, sigma_final, count = c(nlon_prior, nlat_prior, ntimes_sig))
  
  #save the sigma file that was created
  nc_close(nc_sigma)

}