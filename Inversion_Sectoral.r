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
#(corresponding to values in first column of vulcan sectors)
OCO3.bkg.threshold <- 0.01
lambda_p <- c(0.95, 0.93, 0.93, 0.9, 1)
S_p <- c(0.03, 0.05, 0.03, 0.03, 0.03)
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
SAMs <- read.csv('ext/filtered_SAMs.csv')
SAMs <- as.POSIXct(SAMs[,1], format = '%Y-%m-%d %H:%M:%H',
                   tz = local.tz)
attr(SAMs, 'tzone') <- 'UTC'

#initialize functions
setwd(file.path(home.dir, work.ext))
source('r/dependencies.r')

############################
### Setup the Dataframes ###
############################
#list the available SAMs
SAM.path <- list.files(file.path(home.dir, data.ext),
                       pattern = 'out_', full.names = TRUE)
filtered.SAM.path <- grep(SAM.path,
                          pattern = paste(strftime(SAMs, format = '%Y%m%d',
                                                   tz = 'UTC'),
                                          collapse = '|'),
                          value = TRUE)
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

########################
### Do New Inversion ###
########################
#make lambda_p
lambda_p <- as.matrix(lambda_p)

#make S_p
tmp.matrix <- matrix(0, nrow = length(S_p), ncol = length(S_p))
diag(tmp.matrix) <- S_p
S_p <- tmp.matrix

#make R
xco2.errors_df <-
  #observed
  (xco2_df$xco2 - (xco2_df$bio + xco2_df$OCO3.bkg)) -
  #modeled
  (rowSums(xco2_df[,gsub(' ', '', priors)]))

spatial.error_df <- data.frame(time = xco2_df$date.time,
                               lon = xco2_df$lon,
                               lat = xco2_df$lat,
                               error = abs(xco2.errors_df))

R <- spatial.correlation(spatial.error = spatial.error_df,
                         vgm.binwidth = 3, vgm.cutoff = 20,
                         included.xco2.errors =
                           'ext/included.errors.csv',
                         plot.output.path = 'Out/new.R/')

K <- as.matrix(xco2_df[,gsub(' ', '', priors)])
z <- xco2_df$xco2 - (xco2_df$bio + xco2_df$OCO3.bkg)

term.1 <- S_p %*% t(K)
term.2 <- solve((K %*% S_p %*% t(K)) + R)
term.3 <- (z - (K %*% lambda_p))

lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
S_error <- solve(t(K) %*% solve(R) %*% K + solve(S_p))

mean(term.3)
mean(z - (K %*% lambda_hat))
