#setup the control data
#load essential libraries
library(ggplot2); library(plotly); library(ggmap)
library(mcr); library(Metrics); library(ape)
library(qcc); library(patchwork)

#setup the input file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/OCO3_LosAngeles/Bayesian_Inversion'

#city center
city.center <- c(-118.138, 34.004)
local.tz <- 'America/Los_Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#background threshold from OCO-3 method
#(corresponding to values in first column of Vulcan sectors)
#Power Industry, Industry, Cement, OnRoad, NonRoad,
#Rail, Aviation, Marine, Commercial, Residential
OCO3.bkg.threshold <- 0.01
lambda_p <- rep(1, 10)
S_p <- c(0.05^2, 0.38^2, 0.10^2, 0.50^2, 0.06^2,
         0.30^2, 1^2, 0.25^2, 0.03^2, 0.001^2)

#read in the prior list
priors <- read.csv('ext/defined_vulcan_sectors.csv')[,1]
priors <- grep(priors, pattern = '_', invert = TRUE, value = TRUE)

#get domain map
#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)
ggmap <- get_map(location = city.center,
                 maptype = 'satellite',
                 color = 'bw',
                 zoom = 9)

#get list of filtered SAMs
# SAMs <- read.csv('ext/filtered_SAMs.csv')
# SAMs <- as.POSIXct(SAMs[,1], format = '%Y-%m-%d %H:%M:%H',
#                    tz = local.tz)
# attr(SAMs, 'tzone') <- 'UTC'

#initialize functions
setwd(file.path(home.dir, work.ext))
source('r/dependencies.r')

############################
### Setup the Dataframes ###
############################
#list the available SAMs
SAM.path <- list.files(file.path(home.dir, data.ext),
                       pattern = 'out_', full.names = TRUE)
# filtered.SAM.path <- grep(SAM.path,
#                           pattern = paste(strftime(SAMs, format = '%Y%m%d',
#                                                    tz = 'UTC'),
#                                           collapse = '|'),
#                           value = TRUE)
filtered.SAM.path <- SAM.path
filtered.SAM.path <- grep(filtered.SAM.path, pattern = 'out_2019',
                          value = TRUE, invert = TRUE)
remove('SAM.path')

#this loop will skip any missing data without crashing the inversion
for(i in 1:length(filtered.SAM.path)) {
  
  #obtain the TCCON signal
  TCCON.signal.path <- file.path(filtered.SAM.path[i], 'include',
                                 'TCCON_signal.csv')
  if(file.exists(TCCON.signal.path)) {
    TCCON.signal <- read.csv(TCCON.signal.path)
  } else{next}
  
  #obtain the TCCON background
  TCCON.background.path <- file.path(filtered.SAM.path[i], 'include',
                                     'TCCON_background.csv')
  if(file.exists(TCCON.background.path)) {
    TCCON.background <- read.csv(TCCON.background.path)
    names(TCCON.background) <- c('TCCON.bkg', 'TCCON.uncert')
  } else{next}
  
  #obtain the OCO3 background (remove bio background)
  OCO3.background.path <- file.path(filtered.SAM.path[i], 'include',
                                    'OCO3_background.csv')
  if(file.exists(OCO3.background.path)) {
    OCO3.background <- read.csv(OCO3.background.path)
    if(is.na(OCO3.background[1,1])) {next}
    OCO3.background$background <-
      OCO3.background$background -
      OCO3.background$bio.background
    OCO3.background <- OCO3.background[,1:2]
    names(OCO3.background) <- c('OCO3.bkg', 'OCO3.uncert')
  } else{next}
  
  xco2.path <- file.path(filtered.SAM.path[i], 'include',
                         'xco2.csv')
  
  if(file.exists(xco2.path)) {
    xco2 <- read.csv(xco2.path)
    if(any(is.na(xco2$xco2))) {print(i); stop()}
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
filtered.SAM.times <- unique(xco2_df$date.time)

########################
### Do New Inversion ###
########################
lambda_hat_df <- data.frame(matrix(NA, nrow = 0,
                                   ncol = (3+length(priors))))
names(lambda_hat_df) <- c('Iteration', 'time', priors, 'mean.error')

error_hat_df <- data.frame(matrix(NA, nrow = 0,
                                  ncol = (2+length(priors))))
names(error_hat_df) <- c('Iteration', 'time', priors)

#make lambda_p
lambda_p <- as.matrix(lambda_p)
remove('lambda_hat'); remove('S_error')
lambda_initial <- lambda_p
for(i in 1:length(filtered.SAM.times)) {
  
  print(paste0('Starting iteration ', i))
  
  sub.xco2_df <-
    subset(xco2_df, date.time == filtered.SAM.times[i])
  
  if(exists('lambda_hat')) {lambda_p <- lambda_hat}
  if(exists('S_error')) {S_p <- as.matrix(S_error)}
  
  #make S_p
  if(!exists('S_error')) {
    tmp.matrix <- matrix(0, nrow = length(S_p), ncol = length(S_p))
    diag(tmp.matrix) <- S_p
    S_p <- tmp.matrix
  }
  
  #scale with new lambda here!
  K <- as.matrix(sub.xco2_df[,gsub(' ', '', priors)])
  scaled.K <- K %*% lambda_p
  
  #make R
  xco2.errors_df <-
    #observed
    (sub.xco2_df$xco2 - (sub.xco2_df$bio + sub.xco2_df$OCO3.bkg)) -
    #modeled (MODIFIED HERE!)
    rowSums(K)
  
  spatial.error_df <- data.frame(time = sub.xco2_df$date.time,
                                 lon = sub.xco2_df$lon,
                                 lat = sub.xco2_df$lat,
                                 error = abs(xco2.errors_df))
  
  R <- spatial.correlation(spatial.error = spatial.error_df,
                           vgm.binwidth = 3, vgm.cutoff = 20,
                           included.xco2.errors =
                             'ext/included.errors.csv',
                           plot.output.path = 'Out/new.R/')
  
  z <- sub.xco2_df$xco2 - (sub.xco2_df$bio + sub.xco2_df$OCO3.bkg)
  
  term.1 <- S_p %*% t(K)
  term.2 <- solve((K %*% S_p %*% t(K)) + R)
  term.3 <- (z - (K %*% lambda_p))
  
  lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
  S_error <- solve(t(K) %*% solve(R) %*% K + solve(S_p))
  
  #add the iterations to the lambda data frame
  add.line <- data.frame(i, filtered.SAM.times[i],
                         as.data.frame(t(lambda_hat)),
                         mean(z - (K %*% lambda_hat)))
  names(add.line) <- names(lambda_hat_df)
  lambda_hat_df <- rbind(lambda_hat_df, add.line)
  
  #add the errors to the lambda error data frame
  add.line <- data.frame(i, filtered.SAM.times[i],
                         t(sqrt(diag(S_error))))
  error_hat_df <- rbind(error_hat_df, add.line)
  
}

plot.data <- data.frame(lambda_hat_df$Iteration,
                        melt(lambda_hat_df[,2:(length(priors)+2)]))
names(plot.data) <- c('SAM', 'time', 'Sector', 'ScalingFactor')

plot.data.error <- data.frame(error_hat_df$i,
                              melt(error_hat_df[,2:(length(priors)+2)]))
names(plot.data.error) <- c('SAM', 'time', 'Sector', 'Error')

plot.data$error <- plot.data.error$Error

write.csv(plot.data, file = 'Out/Results/Bayesian_Results.csv',
          row.names = FALSE)

ggplot() +
  geom_ribbon(data = plot.data,
              aes(x = SAM,
                  ymin = ScalingFactor - error,
                  ymax = ScalingFactor + error),
              fill = 'lightblue', alpha = 0.5) +
  geom_line(data = plot.data,
            aes(x = SAM, y = ScalingFactor)) +
  theme_classic() +
  facet_wrap(Sector ~ ., ncol = 2, scales = 'free_y')

#save the S_error
write.csv(S_error, file = 'Out/Results/S_posterior.csv')

#build correlation matrix
diag_vals <- sqrt(diag(S_error))
diag_matrix <- matrix(0, nrow = length(diag_vals), ncol = length(diag_vals))
diag(diag_matrix) <- diag_vals
r_matrix <- solve(diag_matrix) %*% S_error %*% solve(diag_matrix)
write.csv(r_matrix, file = 'Out/Results/r_matrix.csv')
