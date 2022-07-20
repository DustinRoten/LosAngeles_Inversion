#load essential libraries
library(ggplot2); library(plotly)
library(mcr); library(Metrics)

#setup the input file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/OCO3_LosAngeles/Bayesian_Inversion'

#city center
city.center <- c(-118.138, 34.004)
local.tz <- 'America/Los_Angeles'

#number of iterations in the regression analysis
iterations <- 1000

#initialize functions
setwd(file.path(home.dir, work.ext))
source('r/dependencies.r')

############################
### Setup the Dataframes ###
############################
#list the available SAMs
SAM.path <- list.files(file.path(home.dir, data.ext),
                       pattern = 'out_',
                       full.names = TRUE)

#this loop will skip any missing data without crashing the inversion
for(i in 1:length(SAM.path)) {
  
  #obtain the TCCON signal
  TCCON.signal.path <- file.path(SAM.path[i], 'include',
                                 'TCCON_signal.csv')
  if(file.exists(TCCON.signal.path)) {
    TCCON.signal <- read.csv(TCCON.signal.path)
  } else{next}
  
  #obtain the TCCON background
  TCCON.background.path <- file.path(SAM.path[i], 'include',
                                     'TCCON_background.csv')
  if(file.exists(TCCON.background.path)) {
    TCCON.background <- read.csv(TCCON.background.path)
    names(TCCON.background) <- c('TCCON.bkg', 'TCCON.uncert')
  } else{next}
  
  #obtain the OCO3 background (remove bio background)
  OCO3.background.path <- file.path(SAM.path[i], 'include',
                                    'OCO3_background.csv')
  if(file.exists(OCO3.background.path)) {
    OCO3.background <- read.csv(OCO3.background.path)
    OCO3.background$background <-
      OCO3.background$background -
      OCO3.background$bio.background
    OCO3.background <- OCO3.background[,1:2]
    names(OCO3.background) <- c('OCO3.bkg', 'OCO3.uncert')
  } else{next}
  
  xco2.path <- file.path(SAM.path[i], 'include',
                         'xco2.csv')
  if(file.exists(xco2.path)) {
    xco2 <- read.csv(xco2.path)
  } else{next}
  
  #format background values
  background <- data.frame(TCCON.background, OCO3.background)
  
  if(i == 1) {
    TCCON.signal_df <- TCCON.signal
    background_df <- background
    xco2_df <- cbind(xco2, OCO3.background)
  } else if(i > 1) {
    TCCON.signal_df <- rbind(TCCON.signal_df, TCCON.signal)
    background_df <- rbind(background_df, background)
    xco2_df <- rbind(xco2_df,
                     cbind(xco2, OCO3.background))
  }
}
############################
############################



############################
### Calculate Extra Shit ###
############################
xco2_df <- subset(xco2_df, !is.na(OCO3.bkg))
xco2_df <- subset(xco2_df, (xco2 - (bio + OCO3.bkg)) > 0)
xco2_df$date.time <- as.POSIXct(xco2_df$date.time,
                                tz = 'UTC')
attr(xco2_df$date.time, 'tzone') <- local.tz

#determine sounding distance from "city center"
xco2_df$dist <- 10^-3*distGeo(p1 = city.center,
                              p2 = xco2_df[,c('lon', 'lat')])

SAM.list <- unique(xco2_df$date.time)

#determine bulk properties
bulk.properties <- data.frame(matrix(NA, nrow = 0, ncol = 4))
names(bulk.properties) <- c('SAM', 'month', 'hour', 'dist_to_city')
for(i in 1:length(SAM.list)) {
  
  SAM <- subset(xco2_df, date.time == SAM.list[i])
  
  add.line <- data.frame(SAM = SAM.list[i],
                         month = month(SAM.list[i]),
                         hour = hour(SAM.list[i]),
                         dist_to_city = mean(SAM$dist))
  
  bulk.properties <- rbind(bulk.properties, add.line)
}
############################
############################

#obtain the list of prior sectors
priors <- read.csv('ext/defined_vulcan_sectors.csv')[1:5,]


#############################
### Calculate the Control ###
#############################
SAMs <- xco2_df

#make lambda_p
lambda_p <- as.matrix(1)

#make S_p
S_p <- as.matrix(runif(1, min = 0, max = 0.15))

#make K
K.all <- as.vector(rowSums(SAMs[,gsub(' ', '', priors[,1])]))

#make z
z.all <- as.vector(SAMs$xco2 - (SAMs$bio + SAMs$OCO3.bkg))

#make R
xco2.errors_df <- z.all - K.all

spatial.error_df <-
  data.frame(time = SAMs$date.time,
             lon = SAMs$lon,
             lat = SAMs$lat,
             error = abs(xco2.errors_df))

R.all <- suppressMessages(spatial.correlation(spatial.error = spatial.error_df,
                                              vgm.binwidth = 3, vgm.cutoff = 20,
                                              included.xco2.errors =
                                                'ext/included.errors.csv',
                                              plot.output.path = 'Out/R'))

#do fancy math
term.1 <- S_p %*% t(K.all)
term.2 <- solve((K.all %*% S_p %*% t(K.all)) + R.all)
term.3 <- (z.all - (K.all %*% lambda_p))

lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
S_error <- solve(t(K.all) %*% solve(R.all) %*% K.all + solve(S_p))
#############################
#############################



collect.SAM.data <- data.frame(matrix(NA, nrow = 0, ncol = 47))
col.names <- c('posterior.mean.err', 'posterior.sd.err',
               'rel.prior.mean.err', 'rel.posterior.mean.err',
               'rel.prior.sd.err', 'rel.posterior.sd.err',
               'lambda_hat', 'S_p', 'S_error',
               'num_of_SAMs', 'num_of_soundings',
               'min.month', 'q1.month', 'median.month',
               'mean.month', 'q3.month', 'max.month',
               'range.month',
               'min.hour', 'q1.hour', 'median.hour',
               'mean.hour', 'q3.hour', 'max.hour',
               'range.hour',
               'min.dist', 'q1.dist', 'median.dist',
               'mean.dist', 'q3.dist', 'max.dist',
               'range.dist',
               'min.dist_to_city', 'q1.dist_to_city',
               'median.dist_to_city', 'mean.dist_to_city',
               'q3.dist_to_city', 'max.dist_to_city',
               'range.dist_to_city',
               'min.signal', 'q1.signal', 'median.signal',
               'mean.signal', 'q3.signal', 'max.signal',
               'range.signal', 'file.path')
#,
#'min.inter_dist', 'q1.inter_dist',
#'median.inter_dist', 'mean.inter_dist',
#'q3.inter_dist', 'max.inter_dist',
#'range.inter_dist')

names(collect.SAM.data) <- col.names

#let the shitshow begin!
for(i in 1:iterations) {
  
  num_of_SAMs <- round(runif(1, min = 5, max = length(SAM.list)))
  
  SAMs_to_use <- sample(SAM.list, num_of_SAMs)
  
  which.soundings <- xco2_df$date.time %in% SAMs_to_use
  SAMs <- xco2_df[which.soundings,]

  #save the SAMs used in the iteration
  iter.file.path <- paste0('Out/SAMs/iter_', i, '.csv')
  write.csv(SAMs, file = iter.file.path,
            row.names = FALSE)
  
  #################
  ### Inversion ###
  #################
  #make lambda_p
  lambda_p <- as.matrix(1)
  
  #make S_p
  S_p <- as.matrix(runif(1, min = 0, max = 0.25))
  
  #make K
  K <- as.vector(rowSums(SAMs[,gsub(' ', '', priors[,1])]))
  
  #make z
  z <- as.vector(SAMs$xco2 - (SAMs$bio + SAMs$OCO3.bkg))
  
  #make R
  xco2.errors_df <- z - K
  
  spatial.error_df <-
    data.frame(time = SAMs$date.time,
               lon = SAMs$lon,
               lat = SAMs$lat,
               error = abs(xco2.errors_df))
  
  R <- suppressMessages(spatial.correlation(spatial.error = spatial.error_df,
                                            vgm.binwidth = 3, vgm.cutoff = 20,
                                            included.xco2.errors =
                                              'ext/included.errors.csv',
                                            plot.output.path = 'Out/R'))
  
  #do fancy math
  term.1 <- S_p %*% t(K)
  term.2 <- solve((K %*% S_p %*% t(K)) + R)
  term.3 <- (z - (K %*% lambda_p))
  
  lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
  S_error <- solve(t(K) %*% solve(R) %*% K + solve(S_p))
  
  posterior.mean.err <- mean(z.all - K.all %*% lambda_hat)
  posterior.sd.err <- sd(z.all - K.all %*% lambda_hat)
  
  rel.prior.mean.err <- mean(term.3)
  rel.posterior.mean.err <- mean(z - (K %*% lambda_hat))
  
  rel.prior.sd.err <- sd(term.3)
  rel.posterior.sd.err <- sd(z - (K %*% lambda_hat))
  
  #determine bulk statistics
  bulk.data <- bulk.properties[bulk.properties$SAM %in% SAMs_to_use,]
  hour.summary <- summary(bulk.data$hour)
  month.summary <- summary(bulk.data$month)
  dist_to_city.summary <- summary(bulk.data$dist_to_city)
  
  dist.summary <- summary(SAMs$dist)
  signal.summary <- summary(SAMs$xco2 - (SAMs$bio + SAMs$OCO3.bkg))
  
  add.line <- data.frame(posterior.mean.err, posterior.sd.err,
                         rel.prior.mean.err, rel.posterior.mean.err,
                         rel.prior.sd.err, rel.posterior.sd.err,
                         lambda_hat[1,1], S_p[1,1], S_error[1,1],
                         num_of_SAMs,
                         num_of_soundings = nrow(SAMs),
                         
                         min.month = as.numeric(month.summary[1]),
                         q1.month = as.numeric(month.summary[2]),
                         median.month = as.numeric(month.summary[3]),
                         mean.month = as.numeric(month.summary[4]),
                         q3.month = as.numeric(month.summary[5]),
                         max.month = as.numeric(month.summary[6]),
                         range.month = as.numeric(month.summary[6] -
                                                  month.summary[1]),
                         
                         min.hour = as.numeric(hour.summary[1]),
                         q1.hour = as.numeric(hour.summary[2]),
                         median.hour = as.numeric(hour.summary[3]),
                         mean.hour = as.numeric(hour.summary[4]),
                         q3.hour = as.numeric(hour.summary[5]),
                         max.hour = as.numeric(hour.summary[6]),
                         range.hour = as.numeric(hour.summary[6] - 
                                                 hour.summary[1]),
                         
                         min.dist = as.numeric(dist.summary[1]),
                         q1.dist = as.numeric(dist.summary[2]),
                         median.dist = as.numeric(dist.summary[3]),
                         mean.dist = as.numeric(dist.summary[4]),
                         q3.dist = as.numeric(dist.summary[5]),
                         max.dist = as.numeric(dist.summary[6]),
                         range.dist = as.numeric(dist.summary[6] -
                                                 dist.summary[1]),
                         
                         min.dist_to_city = as.numeric(dist_to_city.summary[1]),
                         q1.dist_to_city = as.numeric(dist_to_city.summary[2]),
                         median.dist_to_city = as.numeric(dist_to_city.summary[3]),
                         mean.dist_to_city = as.numeric(dist_to_city.summary[4]),
                         q3.dist_to_city = as.numeric(dist_to_city.summary[5]),
                         max.dist_to_city = as.numeric(dist_to_city.summary[6]),
                         range.dist_to_city = as.numeric(dist_to_city.summary[6] -
                                                         dist_to_city.summary[1]),
                         
                         min.signal = as.numeric(signal.summary[1]),
                         q1.signal = as.numeric(signal.summary[2]),
                         median.signal = as.numeric(signal.summary[3]),
                         mean.signal = as.numeric(signal.summary[4]),
                         q3.signal = as.numeric(signal.summary[5]),
                         max.signal = as.numeric(signal.summary[6]),
                         range.signal = as.numeric(signal.summary[6] -
                                                   signal.summary[1]),
                         file.path = iter.file.path)
  
  collect.SAM.data <- rbind(collect.SAM.data, add.line)
  
  #message to user
  cat(paste0(round(100*i/iterations, 2), '% complete.     ', '\r'))
  
}; names(collect.SAM.data) <- col.names

#save the file for use later
write.csv(collect.SAM.data, file = 'MC_Output.csv',
          row.names = FALSE)
#################
#################