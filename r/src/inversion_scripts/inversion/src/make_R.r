# Assemble the model data mismatch matrix
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptor.s
##  make_z.rds
##
## output files:
##  R.rds - file containing a (#obs x #obs) matrix with values of model data mismatch
## (error covariance of observational constraints)

#load required libraries
library(lubridate)

# run dependent scripts
source("config_R_uncert.r")

## add_off_diags() Function declaration:
## add off-diagonals to the model data mismatch matrix for a particular site and
## a particular R-component
## Inputs:
# R - the existing model-data mismatch matrix you wish to add off-diagonals to
# site_vec - vector of site strings (or filepaths) corresponding to the receptor list
# datetimes - vector of POSIX times corresponding to the receptor list
# site - string corresponding to the site desired to add off-diagonals for
# rmse - the RMSE value of the particular R component
# correlation_scale - temporal correlation scale (in days) for this component
# non_decaying - boolean. Do you want to treat the correlation as non-decaying?
#                (i.e. rmse is constant down the off-diagonals)
# correlate_btwn_days - boolean. Do you want to extend the correlation of errors
#                       across all days (T) or limit to obs within same day (F)?

add_off_diags <- function(R, site_vec, datetimes, site, rmse, correlation_scale,
                          non_decaying = F, correlate_btwn_days = F) {
  
  # get the indices on the diagonal that correspond to this site
  # ORIGINAL CODE: siteTF <- grepl(site, site_vec)
  siteTF <- grepl(paste0(site, '/'), site_vec)
  
  # unique list of the dates that appear in the receptor list
  posix_dates <- unique(ISOdate(year(datetimes), month(datetimes),
                                day(datetimes), tz = "UTC"))
  
  # list of each date-time's matching days for this site
  isitedays <-
    lapply(posix_dates,
           FUN = function(x) which(year(datetimes) %in% year(x) &
                                     month(datetimes) %in% month(x) &
                                     day(datetimes) %in% day(x) & siteTF))
  
  # form the index-pairs matrix (will append to this next)
  match_pairs <- array(NA, dim = c(0,2))
  if(correlate_btwn_days){
    # if correlating between all days, get the site-matching indices
    match_pairs <- expand.grid(which(siteTF), which(siteTF))
  } else{
    # otherwise get those pairs which are only within the given day
    for(ii in 1:length(isitedays)){
      isiteday <- isitedays[[ii]]
      match_pairs <- rbind(match_pairs, expand.grid(isiteday, isiteday))
    }
  }
  
  # define the diagonal indices, which should NOT be altered
  idiag <- which(apply(match_pairs, FUN = function(x) x[1] == x[2], MARGIN = 1))
  irow <- match_pairs[-idiag, 1] #filter out diagonals
  icol <- match_pairs[-idiag, 2]
  ioff_diag <- cbind(irow, icol)
  
  if(!non_decaying){
    # get a matrix of nobs x nobs with time differences of receptors, in secs
    t_diff <- sapply(datetimes, FUN = function(x) abs(datetimes - x))
    R_Xt <- t_diff * (1/(24 * 3600)) #convert from seconds to days
    R_Dt <- exp(-R_Xt/correlation_scale) # establish the decay weights
    R_decay <- R_Dt * rmse^2 #mutliply decay weights by the rmse error
    R[ioff_diag] <- R[ioff_diag] + R_decay[ioff_diag] #add to off-diags of R
  } else{
    R[ioff_diag] <- R[ioff_diag] + rmse^2 #apply the unaltered RMSE value to the off-diagonals desired
  }
  #return R
  R
}

## aggregate_R Function declaration:
## aggregate the R matrix to daily observations (errors are generally reduced because
## 1 daily representation of mixing ratio is more accurate than individual obs)
## Inputs:
# R - the existing model-data mismatch matrix you wish to aggregate
# site_vec - vector of site strings corresponding to the receptor list
# datetimes - vector of POSIX times corresponding to the receptor list
aggregate_R <- function(R, site_vec, datetimes, site_vec_aggr, dates_aggr){
  
  nobs <- length(site_vec) #get number of total obs, pre-aggregation
  nobs_aggr <- length(site_vec_aggr) #number of total obs after aggregation
  R_aggr_diag <- rep(0, nobs_aggr) #this will hold the new aggregated R
  
  # vector of unique dates that appear in the receptor list
  posix_dates <- unique(ISOdate(year(datetimes), month(datetimes),
                                day(datetimes), tz = "UTC"))
  
  ndays <- length(posix_dates)
  
  # loop through the sites
  for(site in unique(site_vec_aggr)){
    
    # get the indices on the diagonal that correspond to this site
    siteTF <- grepl(paste0(site, '/'), site_vec)
    siteTF_aggr <- grepl(paste0(site, '/'), site_vec_aggr)
    
    # list of each date-time's matching days for this site
    isitedays <-
      lapply(posix_dates,
             FUN = function(x) which(year(datetimes) %in% year(x) &
                                       month(datetimes) %in% month(x) &
                                       day(datetimes) %in% day(x) & siteTF))
    
    # same as above, but for the aggregated matrix
    isitedays_aggr <-
      lapply(posix_dates,
             FUN = function(x) which(year(dates_aggr) %in% year(x) &
                                       month(dates_aggr) %in% month(x) &
                                       day(dates_aggr) %in% day(x) & siteTF_aggr))
    
    
    # loop through all days and aggregate
    for (ii in 1:ndays) {
      
      idaysite <- isitedays[[ii]] #indices of this site's obs on this day
      iaggr <- isitedays_aggr[[ii]]
      
      # if num obs for this site on this day is below minimum defined in config.r, skip
      if (length(idaysite) < min_agg_obs)
        next
      
      # establish an aggregation operator for this site's obs on this day
      # in general (if no missing data/footprints) it's 1/length(subset_hours_utc)
      W <- rep(0, nobs)
      W[idaysite] <- 1/length(idaysite)
      new_diag <- W %*% R #flatten R into a vector, where only vals from this site/day are non-zero
      
      R_chunk <- new_diag[idaysite]
      
      # NOTE: take mean
      R_aggr_diag[iaggr] <- mean(R_chunk)
    }
  } #end sites for-loop
  
  R_aggregated <- diag(R_aggr_diag)
  
  #return R_aggr
  R_aggregated
}

## Initialize a plotting function for spatial R plots when spatially
## correlated off diagonals are present
plot_spatial_off_diags <- function(R, R_ls, variogram.model, variogram.fit,
                                   vgm.binwidth, vgm.cutoff, output.path = 'out') {
  
  require(ggplot2); require(patchwork)
  
  #convert the incoming R matrix into a dataframe
  R_df <- reshape2::melt(R)
  
  #plot the R matrix
  R.plot <- ggplot() +
    ggtitle(expression(paste('Spatial Correlation Matrix of XCO'[2]))) +
    xlab('Columns') +
    ylab('Rows') +
    geom_tile(data = R_df, aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradient(low = 'white', high = 'red',
                        name = 'RMSE (ppm)') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom')
  
  #plot the variogram fit
  fitted.line <- variogramLine(variogram.fit, maxdist = vgm.cutoff)
  
  #plot the variogram
  vgm.plot <- ggplot() +
    ggtitle(expression(paste('Spatial Correlation of XCO'[2]))) +
    labs(subtitle = 'Exponential Semivariogram') +
    xlab('Distance (km)') +
    ylab(expression(paste(gamma, ' (ppm'^2, ')'))) +
    geom_point(data = variogram.model,
               aes(x = dist, y = gamma)) +
    geom_line(data = fitted.line,
              aes(x = dist, y = gamma),
              color = 'blue') +
    geom_text(aes(x = Inf, y = -Inf,
                  label = paste0('Range = ', round(variogram.fit[2,3], 2),
                                 'km', '\n',
                                 if(variogram.fit[2,3]/3 <= vgm.cutoff) {
                                   paste0('L = ', round(variogram.fit[2,3]/3, 2))
                                 } else {
                                   paste0('Overridden L = ', round(R_ls, 2),
                                          'km', '\n',
                                          'L = ', round(variogram.fit[2,3]/3, 2))
                                 },
                                 'km', '\n',
                                 'Bin Width = ', vgm.binwidth, 'km', '\n')),
              hjust = 1, vjust = 0) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  #combine and save the two plots
  combined.plots <- vgm.plot | R.plot
  ggsave(combined.plots,
         filename = file.path(output.path,
                              'SpatialCorrelation_XCO2.jpg'),
         device = 'jpg', width = 7, height = 4, units = 'in')
}

## Spatial R function declaration:
## apply an exponentially decaying correlation coefficient across dense observations
## Inputs:
# R, xco2 observations, z, and receptors.
add_spatial_off_diags <- function(R, obs.path, z.path, receptors.path,
                                  output.path) {
  
  #check that the files exist
  obs.exists <- file.exists(obs.path)
  z.exists <- file.exists(z.path)
  receptors.exists <- file.exists(receptors.path)
  all.present <- all(c(obs.exists, z.exists, receptors.exists))
  
  #stop if any files are missing
  if(!all.present) stop('Required files are missing!')
  
  #read in all of the appropriate files
  obs <- read.csv(obs.path)  
  z <- readRDS(z.path)
  receptors <- readRDS(receptors.path)
  
  #create a dataframe to store matched observations
  reordered.obs <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  names(reordered.obs) <- c('lon', 'lat', 'diff')
  for(i in 1:nrow(obs)) {
    
    #' The column of the `receptors`, `z`, and `obs` dataframes
    #' should be identical. Just in case, this for loop will
    #' step through and ensure that they are matched up. The
    #' output dataframe will be used in the spatial correlation
    #' calculations.
    obs.idx <- grep(pattern = receptors[i,2],
                    x = obs[,7])
    
    #calculated with the differences in XCO2
    add.lines <- data.frame(lon = obs$lon[obs.idx],
                            lat = obs$lat[obs.idx],
                            diff = obs$modeled.xco2[obs.idx] - 
                              obs$xco2[obs.idx])
    reordered.obs <- rbind(reordered.obs, add.lines)
    
  }
  
  #List the distances between all soundings
  aggregated.dist.diffs <- NULL
  for(j in 1:nrow(reordered.obs)) {
    dist.diffs <- geodist::geodist(x = reordered.obs[j,1:2],
                                   y = reordered.obs[,1:2],
                                   measure = 'haversine')
    aggregated.dist.diffs <- c(aggregated.dist.diffs,
                               dist.diffs)
  }
  
  #determine the median and min distances across the SAM
  median.dist <-
    summary(aggregated.dist.diffs)[3]/1000
  min.dist <-
    min(aggregated.dist.diffs[aggregated.dist.diffs != 0])/1000
  
  #' The `binwidth`, `cutoff`, and `bin number` values are
  #' now determined on the fly.
  #' binwidth is 3 times the minimum non-zero distance.
  #' the cutoff length is the median distance between soundings.
  vgm.binwidth <- ceiling(3*min.dist)
  vgm.cutoff <- median.dist
  vgm.bins <- ceiling(vgm.cutoff/vgm.binwidth)
  
  #build a distance matrix and convert to km instead of meters.
  X_R <- geodist::geodist(x = reordered.obs[,1:2],
                          measure = 'haversine')
  X_R <- X_R/1000
  
  #prepare the dataframe for the variogram calculataion
  coordinates(reordered.obs) <-  ~lon+lat
  proj4string(reordered.obs) <- CRS('+proj=longlat +datum=WGS84')
  
  if(nrow(obs) > 3) { #this loop is for testing purposes
    #calculate the variogram, fit the model
    vgm <- variogram(diff~lon+lat, data = reordered.obs,
                     width = vgm.binwidth, 
                     cutoff = vgm.cutoff, cressie = F)
    model.fit <- fit.variogram(vgm, vgm('Exp'))
    
    #pull the range value from the variogram fit and solve for L.
    #3*L=range
    #' If the calculated L is larger than the median distance, use the
    #' median distance instead.
    if(model.fit[2,3]/3 <= vgm.cutoff) {
      R_ls <- model.fit[2,3]/3
    } else {
      R_ls <- vgm.binwidth
    }
  } else {
    ### All of this can be removed after testing!
    #perform the spatial decay calculation
    X_R <- exp(-X_R/7)
    
    #get the original diagonal values (summed RMSE) from the R matrix
    diag.value <- unique(diag(R))
    if(length(diag.value) > 1) stop('non-unique values on the R diagonal!')
    #make the new R matrix
    R_Xs <- diag.value*X_R
    return(R_Xs)
  }

  #perform the spatial decay calculation
  X_R <- exp(-X_R/R_ls)
  
  #get the original diagonal values (summed RMSE) from the R matrix
  diag.value <- unique(diag(R))
  if(length(diag.value) > 1) stop('non-unique values on the R diagonal!')
  #make the new R matrix
  R_Xs <- diag.value*X_R
  
  #plot output values
  plot_spatial_off_diags(R = R_Xs, R_ls, variogram.model = vgm,
                         variogram.fit = model.fit,
                         vgm.binwidth, vgm.cutoff, output.path = 'out')
  
  return(R_Xs)
}

# ~~~~~~~~~~~~~~~~~~~~ Begin main script ~~~~~~~~~~~~~~~~~~~~ #
# load in receptor info
recep_filepath <- paste0(out_path, "receptors.rds")
receptors <- readRDS(recep_filepath)
recep_times <- as.numeric(receptors[, 1])
recep_names <- receptors[, 2]
nobs <- nrow(receptors)

# convert receptor times from seconds-since-epoch to POSIX
class(recep_times) <- c("POSIXt", "POSIXct")
attributes(recep_times)$tzone <- "UTC"

# get each site's indices in receptor list
isites <- lapply(paste0(sites, '/'),
                        FUN = function(x) grep(x, recep_names))

# load grid info
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# load in z file - anthropogenic enhancements
z_file <- paste0(out_path, "z.rds")
z <- readRDS(z_file)

# Model Data Mismatch is made up of multiple parts R = R_part + R_aggr + R_eddy +
# R_bkgd + R_transPBL + R_transWIND + R_instr + R_bio + R_other

# ~~~~~~~~~~~~~~~ R_part ~~~~~~~~~~~~~~~#
R_part <- rep(rmse_part^2, nobs)

# ~~~~~~~~~~~~~~~ R_aggr ~~~~~~~~~~~~~~~#
R_aggr <- rep(rmse_aggr^2, nobs)

# ~~~~~~~~~~~~~~~ R_eddy ~~~~~~~~~~~~~~~#
R_eddy <- rep(rmse_eddy^2, nobs)

# ~~~~~~~~~~~~~~~ R_bkgd ~~~~~~~~~~~~~~~#
R_bkgd <- rep(rmse_bkgd^2, nobs)

# ~~~~~~~~~~~~~~~ R_transPBL ~~~~~~~~~~~~~~~#
R_transPBL <- rep(rmse_transPBL^2, nobs)

# ~~~~~~~~~~~~~~~ R_transWIND ~~~~~~~~~~~~~~~#
R_transWIND <- rep(rmse_transWIND^2, nobs)

# ~~~~~~~~~~~~~~~ R_instr ~~~~~~~~~~~~~~~#
R_instr <- rep(rmse_instr^2, nobs)

# if extra instrument error for specific sites
for(ii in 1:nsites){
  extra_site_rmse <- extra_sites_err_rmse[[ii]]
  isite <- isites[[ii]]
  R_instr[isite] = extra_site_rmse^2 #R_instr[isite] + extra_site_rmse^2
}

# ~~~~~~~~~~~~~~~ R_bio ~~~~~~~~~~~~~~~#
R_bio <- rep(rmse_bio^2, nobs)

# ~~~~~~~~~~~~~~~ R_other ~~~~~~~~~~~~~~~#
R_other <- rep(rmse_other^2, nobs)

# ~~~~~~~~~~~~~~~ Add up all components of R ~~~~~~~~~~~~~~~#

R <- diag(R_part + R_eddy + R_aggr + R_bkgd + R_transPBL +
            R_transWIND + R_instr + R_bio + R_other)

# ~~~~~~ Now add covariance (off-diagonals) ~~~~~~~#
ncomponents <- nrow(R_arr)

#add temporally related off diagonals
for(ii in 1:ncomponents) {
  
  rmse_in <- as.numeric(R_arr[ii, "RMSE"])
  cor_scale_in <- as.numeric(R_arr[ii, "correlation_scale"])
  non_decaying_in <- as.logical(R_arr[ii, "non_decaying?"])
  cor_btwn_days_in <- as.logical(R_arr[ii, "correlate_btwn_days?"])
  
  if(is.na(cor_scale_in) & !non_decaying_in)
    next
  
  for(site in sites){
    R <- add_off_diags(R, recep_names, recep_times, site, rmse_in, cor_scale_in,
                       non_decaying_in, cor_btwn_days_in)
  }
}

#' If dense XCO2 observations are being used (XCO2), apply spatial correlations
#' to the R matrix as well
if(file.exists('X_config.RData')) {
  
  #' In later iterations of this code, the following parameters can be
  #' included in the main script and passed to the inversion code via
  #' the `setup.config()` function. Now, it is assumed these file paths
  #' will not change.
  obs.path <- '../include/xco2_observations.csv'
  z.path <- 'out/z.rds'
  receptors.path <- 'out/receptors.rds'
  vgm.cutoff <- 30; vgm.bins <- 15
  output.path <- 'out'
  
  #' Create and plot the variogram and new R matrix using the nested
  #' functions from above.
  new.R <- add_spatial_off_diags(R, obs.path, z.path, receptors.path,
                                 output.path)
  
  #overwrite the original R matrix with the new one
  R <- new.R
}


# now aggregate into daily observational errors if aggregate_obs is toggled
if (aggregate_obs) {
  
  recep_aggr_filepath <- paste0(out_path, "receptors_aggr.rds")
  receps_aggr <- readRDS(recep_aggr_filepath)
  recep_sites_aggr <- receps_aggr[,1]
  recep_times_aggr <- receps_aggr[,2]
  class(recep_times_aggr) <- c("POSIXt", "POSIXct")
  attributes(recep_times_aggr)$tzone <- "UTC"
  
  # turn into diagonal n x n matrix
  R_aggregated <-
    aggregate_R(R, recep_names, recep_times, recep_sites_aggr, recep_times_aggr)
  
  # save R to file
  print("saving model data mismatch matrix to R.rds")
  filepath <- paste0(out_path, "R.rds")
  saveRDS(R_aggregated, filepath)
  
} else {
  
  # save R to file
  print("saving model data mismatch matrix to R.rds")
  filepath <- paste0(out_path, "R.rds")
  saveRDS(R, filepath)
  
}
