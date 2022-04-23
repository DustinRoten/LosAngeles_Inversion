# convert posterior optimized emissions and uncertainty to netcdf format
# author: Lewis Kunik

## prerequisite scripts:
##	inversion.r
##	make_Vshat.r
##
## output files:
## 	posterior.nc - netcdf file formatted to give posterior emissions
##	as well as posterior uncertainty in lon x lat x time format
##

#load package dependencies
library(ncdf4)

#load the appropriate configuration file
if(file.exists("X_config.RData")) {
  load('X_config.RData')
} else {source('config.r')}

source("src/prior_uncert_func.r") #function to extract prior gridded uncertainty
source("src/post_uncert_func.r") #function to extract posterior gridded uncertainty

# specify start/end from config.r variables
y1 <- flux_year_start; y2 <- flux_year_end
m1 <- flux_month_start; m2 <- flux_month_end
d1 <- flux_day_start; d2 <- flux_day_end
h1 <- flux_hour_start; h2 <- flux_hour_end
mn1 <- flux_min_start; mn2 <- flux_min_end

# establish the centers of time steps
# line up times exactly for satellite data
if(file.exists('X_config.RData')) {
  flux_times <-
    seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC"),
        to = ISOdatetime(y2, m2, d2, h2, mn2, 0, tz = "UTC"),
        by = flux_t_res) #flux_t_res defined in config.r
} else {
  flux_times <-
    seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC") + flux_t_res/2,
        to = ISOdatetime(y2, m2, d2, h2, mn2, 0, tz = "UTC"),
        by = flux_t_res) #flux_t_res defined in config.r
}

# load spatial covariance matrix E
sp_cov_file <- paste0(out_path, "sp_cov.rds")
E <- readRDS(sp_cov_file)

# load temporal covariance matrix D
tmp_cov_file <- paste0(out_path, "tmp_cov.rds")
D <- readRDS(tmp_cov_file)

# load prior uncertainty and convert to ntime x ncell matrix
sigma_file <- paste0(out_path, "sigma.rds")
sigma <- readRDS(sigma_file)
sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T) # ntimes defined in config.r

# load model data mismatch matrix
R_file <- paste0(out_path, "R.rds")
R <- readRDS(R_file)
nobs <- nrow(R)

# load HQHtranspose file
HQHt_file <- paste0(out_path, "HQHt.rds")
HQHt <- readRDS(HQHt_file)

# load lonlat_domain grid info
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# load posterior emissions
shat_file <- paste0(out_path, "s_hat.rds")
shat_vec <- readRDS(shat_file)
shat_mat <- matrix(shat_vec, nrow = ntimes, byrow = T)

# ~~~~~~~~~~~~~~~~~~~~~~ set up grid formatting ~~~~~~~~~~~~~~~~~~~~~~~~#
# define smallest possible rectangle to hold the posterior grid
lon <- round(seq(min(lonlat_domain[,1]), max(lonlat_domain[,1]), by = lon_res), round_digs)
lat <- round(seq(min(lonlat_domain[,2]), max(lonlat_domain[,2]), by = lat_res), round_digs)
lonlat_post <- as.matrix(expand.grid(lon, lat)) #get all lon-lat pairs in this new grid

nlon <- length(lon)
nlat <- length(lat)

# get index mapping of non-zeo domain cells in the new rectangular grid
iDomain <- apply(lonlat_domain,
                 FUN = function(x) which.min(sqrt((lonlat_post[, 1] - x[1])^2 +
                                                    (lonlat_post[, 2] - x[2])^2)),
                 MARGIN = 1)

# Define aggregation operator to weight each cell by its fractional area
cellareas <- abs(2 * pi * (6371009^2) * (sin((lonlat_domain[, 2] - lat_res/2) *
                                               pi/180) - sin((lonlat_domain[, 2] + lat_res/2) * pi/180))/(360/lon_res))
total_area <- sum(cellareas)
W <- matrix(cellareas/total_area, nrow = 1)  #W is the 'aggregation operator'

# set up prior/posterior gridded uncertainty array for each time step
prior_uncert_mat <- array(0, dim = c(ntimes, ncells))
post_uncert_mat <- array(0, dim = c(ntimes, ncells))

# ~~~~~~~~~~~~~~~~~~ obtain posterior gridded uncert ~~~~~~~~~~~~~~~~~~~#
# loop through times to get gridded uncertainties of each time step
for(ii in 1:ntimes) {
  print(paste("extracting gridded uncertainty of timestep", ii))
  # use the provided functions to obtain prior/posterior uncertainty mat for this time
  Qsum_time <- get_prior_uncert_time(ii, ii, E, D, sigma)
  prior_uncert_mat[ii, ] <- W %*% Qsum_time
  Vshat_time <- get_post_uncert_time(ii, ii, E, D, sigma, R, HQHt, lonlat_domain)
  post_uncert_mat[ii, ] <- W %*% Vshat_time #get the gridded uncertainty for this time
}

# set up 2d arrays with NA values (same dims for all)
post_grid_2d <- array(NA, dim = c(ntimes, nrow(lonlat_post)))
uncert_grid_2d <- array(NA, dim = c(ntimes, nrow(lonlat_post)))
uncert.prior_grid_2d <- array(NA, dim = c(ntimes, nrow(lonlat_post)))


# assign non-NA values to the proper spots in arrays
post_grid_2d[, iDomain] <- shat_mat
uncert_grid_2d[, iDomain] <- post_uncert_mat
uncert.prior_grid_2d[, iDomain] <- prior_uncert_mat


# set up 3-d arrays with NA values
post_grid_3d <- array(NA, dim = c(nlon, nlat, ntimes))
uncert_grid_3d <- array(NA, dim = c(nlon, nlat, ntimes))
uncert.prior_grid_3d <- array(NA, dim = c(nlon, nlat, ntimes))

# loop through times to assign lon x lat matrices to time dimension
for (ii in 1:ntimes) {
  
  post_grid_3d[, , ii] <- t(matrix(post_grid_2d[ii,], nrow = nlat, byrow = T))
  uncert_grid_3d[, , ii] <- t(matrix(uncert_grid_2d[ii,], nrow = nlat, byrow = T))
  uncert.prior_grid_3d[, , ii] <- t(matrix(uncert.prior_grid_2d[ii,], nrow = nlat, byrow = T))
  
}


# ~~~~~~~~~~~~~~~~~ save to netcdf file ~~~~~~~~~~~~~~~~~~~ #
# define dimension variables
time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                      as.numeric(flux_times),
                      longname = "center of time-step, seconds since R epoch: 1970-01-01 00:00:00")
lat_dim <- ncdim_def("lat", "degrees_north", lat, longname="latitude (center of cell)")
lon_dim <- ncdim_def("lon", "degrees_east", lon, longname="longitude (center of cell)")

# define emissions and uncertainty variables
post_emiss_var <- ncvar_def("emiss", flux_units, list(lon_dim, lat_dim, time_dim),
                            longname="gridded posterior emissions, optimized")
uncert_var <- ncvar_def("uncertainty", flux_units, list(lon_dim, lat_dim, time_dim),
                        longname="gridded covariance-aggregated posterior uncertainty")
uncert.prior_var <- ncvar_def("uncertainty", flux_units, list(lon_dim, lat_dim, time_dim),
                              longname = "gridded covariance-aggregated prior uncertainty")

# define list of variables to save
post_vars_emiss <- list(post_emiss_var)
post_vars_uncert <- list(uncert_var)
prior_vars_uncert <- list(uncert.prior_var)

### open files for writing
#file names
nc_filename <- paste0(out_path, "posterior.nc")
nc_filename.uncert <- paste0(out_path, "uncertainty.nc")
nc_filename.uncert.prior <- paste0(out_path, 'uncertainty_prior.nc')

#create the files
nc_post_emiss <- nc_create(nc_filename, post_vars_emiss)
nc_post_uncert <- nc_create(nc_filename.uncert, post_vars_uncert)
nc_post_uncert.prior <- nc_create(nc_filename.uncert.prior, prior_vars_uncert)

#add to files
ncvar_put(nc_post_emiss, post_emiss_var, post_grid_3d, count = c(nlon, nlat, ntimes)) #add emissions
ncvar_put(nc_post_uncert, uncert_var, uncert_grid_3d, count = c(nlon, nlat, ntimes)) #add uncertainty
ncvar_put(nc_post_uncert.prior, uncert.prior_var, uncert.prior_grid_3d, count = c(nlon, nlat, ntimes)) #add prior uncertainty

#close netcdf files
nc_close(nc_post_emiss); nc_close(nc_post_uncert); nc_close(nc_post_uncert.prior)