#load essential libraries
library(ggplot2); library(mcr); library(Metrics)

#setup the input file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/OCO3_LosAngeles/Bayesian_Inversion'

inversion.breaks <- c('20190814', '20200401', '20200505', '20201231')
#inversion.breaks <- c('20190814', '20211231')

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



#########################
### Preliminary Plots ###
#########################
#plot the TCCON/OCO-3 background comparison
sub.background_df <- background_df[complete.cases(background_df),]
lm.fit.1 <- mcreg(x = sub.background_df$TCCON.bkg,
                  y = sub.background_df$OCO3.bkg,
                  method.reg = 'Deming')
lm.fit.1.sum <- getCoefficients(lm.fit.1)
ggplot() +
  ggtitle('OCO-3 and TCCON Background Estimates') +
  xlab('TCCON [ppm]') +
  ylab('OCO-3 [ppm]') +
  geom_label(aes(x = -Inf, y = Inf),
             hjust = 0, vjust = 1,
             label = paste0('m=', round(lm.fit.1.sum[2], 2), '\n',
                            'b=', round(lm.fit.1.sum[1], 2), 'ppm', '\n',
                            'RMSE=', round(rmse(sub.background_df$TCCON.bkg,
                                                sub.background_df$OCO3.bkg), 2),
                            'ppm')) +
               geom_abline(slope = 1,
              intercept = 0,
              linetype = 'dashed') +
  geom_abline(slope = lm.fit.1.sum[2],
              intercept = lm.fit.1.sum[1],
              linetype = 'solid', color = 'red') +
  geom_errorbar(data = sub.background_df,
                aes(x = TCCON.bkg,
                    ymax = OCO3.bkg + OCO3.uncert,
                    ymin = OCO3.bkg - OCO3.uncert)) +
  geom_errorbarh(data = sub.background_df,
                 aes(y = OCO3.bkg,
                     xmax = TCCON.bkg + TCCON.uncert,
                     xmin = TCCON.bkg - TCCON.uncert)) +
  geom_point(data = sub.background_df,
             aes(x = TCCON.bkg, y = OCO3.bkg)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#plot the TCCON/OCO-3 signal comparison
sub.signal_df <- TCCON.signal_df[complete.cases(TCCON.signal_df),]
lm.fit.2 <- mcreg(x = sub.signal_df$signal,
                  y = sub.signal_df$oco3.signal,
                  method.reg = 'Deming')
lm.fit.2.sum <- getCoefficients(lm.fit.2)
ggplot() +
  ggtitle('OCO-3 and TCCON Urban Estimates') +
  xlab('TCCON [ppm]') +
  ylab('OCO-3 [ppm]') +
  geom_label(aes(x = -Inf, y = Inf),
             hjust = 0, vjust = 1,
             label = paste0('m=', round(lm.fit.2.sum[2], 2), '\n',
                            'b=', round(lm.fit.2.sum[1], 2), 'ppm', '\n',
                            'RMSE=', round(rmse(sub.signal_df$signal,
                                                sub.signal_df$oco3.signal), 2),
                            'ppm')) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 'dashed') +
  geom_abline(slope = lm.fit.2.sum[2],
              intercept = lm.fit.2.sum[1],
              linetype = 'solid', color = 'red') +
  geom_errorbar(data = sub.signal_df,
                aes(x = signal,
                    ymax = oco3.signal + oco3.uncert,
                    ymin = oco3.signal - oco3.uncert)) +
  geom_errorbarh(data = sub.signal_df,
                 aes(y = oco3.signal,
                     xmax = signal + uncert,
                     xmin = signal - uncert)) +
  geom_point(data = sub.signal_df,
             aes(x = signal, y = oco3.signal)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#plot the model-obs errors (both backgrounds)
all.backgrounds_df <- xco2_df[complete.cases(xco2_df[,c('background',
                                                        'OCO3.bkg')]),]
#read in the vulcan sectors
priors <- read.csv('ext/defined_vulcan_sectors.csv')
which.p <- grep(priors$category, pattern = '_', invert = TRUE)
priors <- priors[which.p,]

ggplot() +
  geom_histogram(data = all.backgrounds_df, color = 'black',
                 aes(((xco2 - bio) - (OCO3.bkg - 0.3)) -
                       rowSums(all.backgrounds_df[,gsub(' ', '', priors[,1])])),
                 fill = 'red', alpha = 0.5) +

  geom_histogram(data = all.backgrounds_df, color = 'black',
                 aes(((xco2 - bio) - background) -
                   rowSums(all.backgrounds_df[,gsub(' ', '', priors[,1])])),
                 fill = 'white', alpha = 0.5)

#########################
#########################



#############################
### Running the Inversion ###
#############################
#ensure that every sounding has a background value
sub.xco2_df <- xco2_df[complete.cases(xco2_df$background),]

#subset signals above the background value
sub.xco2_df <- subset(sub.xco2_df,
                      (xco2 - (bio + OCO3.bkg) > 0))

#obtain the sector and error information
sectors <- read.csv('ext/defined_vulcan_sectors.csv')
xco2.errors <- read.csv('ext/included.errors.csv')

#loop through the time breaks
inversion.breaks <- as.POSIXct(inversion.breaks,
                               format = '%Y%m%d',
                               tz = 'UTC')

#change the date.time values in sub.xco2_df to POSIX
sub.xco2_df$date.time <- as.POSIXct(sub.xco2_df$date.time,
                                    origin = '1970-01-01',
                                    tz = 'UTC')

for(i in 1:(length(inversion.breaks) - 1)) {
  
  #read in the corresponding prior file
  priors <- read.csv(paste0('ext/S_prior.', i, '.csv'))
  
  #check that the sector names match
  if(any(sectors[,1] != priors[,1]))
    stop('Check category names in /ext/ files!')
  
  #currently, the large point source category isn't needed.
  which.p <- grep(priors$category, pattern = '_', invert = TRUE)
  priors <- priors[which.p,]
  
  #make lambda_p
  lambda_p <- 1
  
  #make S_p
  S_p <- 1-0.75
  
  #get all of the soundings within the timeframe
  time.sub.xco2_df <- subset(sub.xco2_df,
                             date.time >= inversion.breaks[i] &
                             date.time < inversion.breaks[i+1])
  
  #make K
  K <- rowSums(time.sub.xco2_df[,gsub(' ', '', priors[,1])])
  
  #make z
  z <-
    (time.sub.xco2_df$xco2 - time.sub.xco2_df$bio) -
    (time.sub.xco2_df$background)
  
  #make R
  xco2.errors_df <-
    #observed
    ((time.sub.xco2_df$xco2 - time.sub.xco2_df$bio) -
    (time.sub.xco2_df$background)) -
    #modeled
    (rowSums(time.sub.xco2_df[,gsub(' ', '', priors[,1])]))
  
  spatial.error_df <-
    data.frame(time = time.sub.xco2_df$date.time,
               lon = time.sub.xco2_df$lon,
               lat = time.sub.xco2_df$lat,
               error = abs(xco2.errors_df))
  
  R <- spatial.correlation(spatial.error = spatial.error_df,
                           vgm.binwidth = 3, vgm.cutoff = 20,
                           included.xco2.errors =
                             'ext/included.errors.csv',
                           plot.output.path = 'Out/R')
  
  #do fancy math
  term.1 <- as.matrix(S_p) %*% t(K)
  term.2 <- solve((K %*% as.matrix(S_p) %*% t(K)) + R)
  term.3 <- (z - (K %*% as.matrix(lambda_p)))
  
  lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
  S_error <- solve(t(K) %*% solve(R) %*% K + solve(S_p))
  
  errors <- data.frame(prior.errors = term.3,
                       posterior.errors = (z - (K %*% lambda_hat)),
                       period = i)
  lambdas <- data.frame(category= 'Total',
                        lambda = lambda_hat,
                        period = i)
  S_errors <- data.frame(category = 'Total',
                         S = diag(S_error),
                         period = i)
  if(i == 1) {
    errors_df <- errors
    lambda_df <- lambdas
    S_errors_df <- S_errors
  }
  if(i > 1) {
    errors_df <- rbind(errors_df, errors)
    lambda_df <- rbind(lambda_df, lambdas)
    S_errors_df <- rbind(S_errors_df, S_errors)
  }
}


stop()

#scaling factors plot
ggplot() +
  ggtitle('Posterior Scaling Factors') +
  xlab('Period') +
  ylab(expression(hat(lambda))) +
  geom_bar(data = lambda_df, stat = 'identity',
           position = 'dodge', color = 'black',
           aes(x = as.character(period), y = lambda,
               fill = category)) +
  geom_errorbar(data = S_errors_df, position = 'dodge',
                size = 0.5,
                aes(x = as.character(period),
                    ymax = lambda_df$lambda + S,
                    ymin = lambda_df$lambda - S,
                    group = category)) +
  scale_fill_discrete(name = 'Sectors') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

#posterior percent difference
ggplot() +
  ggtitle('Posterior Percent Difference') +
  xlab('Time Period') +
  ylab('Percent Difference [%]') +
  geom_bar(data = lambda_df, stat = 'identity',
           position = 'dodge', color = 'black',
           aes(x = as.character(period), y = 100*(lambda-1),
               fill = category)) +
  geom_errorbar(data = S_errors_df, position = 'dodge',
                size = 0.5,
                aes(x = as.character(period),
                    ymax = 100*((lambda_df$lambda-1) + S),
                    ymin = 100*((lambda_df$lambda-1) - S),
                    group = category)) +
  scale_fill_discrete(name = 'Sectors:') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

#error plot
ggplot() +
  geom_histogram(data = errors_df,
                 aes(prior.errors),
                 color = 'black',
                 fill = 'white',
                 alpha = 1) +
  geom_histogram(data = errors_df,
                 aes(posterior.errors),
                 color = 'black',
                 fill = 'red',
                 alpha = 0.5) +
  scale_y_continuous(trans = 'log10') +
  theme_classic() +
  facet_wrap(. ~ period, ncol = 1)


  