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
OCO3.bkg.threshold <- 0.01
lambda_p <- 1
S_p <- 0.12

#get domain map
#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)
ggmap <- get_map(location = city.center,
                 maptype = 'satellite',
                 color = 'bw',
                 zoom = 9)

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
############################
############################



#######################################
### Analysis of Monte Carlo Results ###
#######################################
priors <- read.csv('ext/defined_vulcan_sectors.csv')[1:5,]

#make K
K.all <- as.vector(rowSums(xco2_df[,gsub(' ', '', priors[,1])]))

#make z
z.all <- as.vector(xco2_df$xco2 - (xco2_df$bio + xco2_df$OCO3.bkg))

#calculate difference
prior_diff <- mean(z.all - K.all)

#read in Monte Carlo data
collect.SAM.data <- read.csv('MC_Output.csv')

#quick formatting change
collect.SAM.data$is.eff <- NA
collect.SAM.data$is.eff[collect.SAM.data$posterior.mean.err <
                          prior_diff] <- 'Yes'
collect.SAM.data$is.eff[collect.SAM.data$posterior.mean.err >=
                          prior_diff] <- 'No'

#lm.fit.1
lm.1 <- lm(posterior.mean.err ~ mean.signal, data = collect.SAM.data)
intercept.1 <- lm.1$coefficients[1]
slope.1 <- lm.1$coefficients[2]
r_val.1 <- sqrt(summary(lm.1)$r.squared)

#regression
signal.plot <- ggplot() +
  ggtitle(paste0('Posterior Error Correction vs.', '\n',
                 'Mean Enhancement')) +
  labs(subtitle = paste0('y=', round(slope.1, 2),
                         'x', if(intercept.1 >= 0) {'+'} else{'-'},
                         round(intercept.1, 2), '; r=',
                         round(r_val.1, 2)),
       caption = paste0('Prior mean error (black dashed line): ',
                        round(prior_diff, 2), 'ppm.')) +
  xlab(expression(bar('z')[i]~'[ppm]')) +
  ylab(expression(bar(z[tot] - bold(K)[tot]*hat(lambda))~'[ppm]')) +
  geom_abline(slope = 0, intercept = prior_diff,
              linetype = 'dashed', color = 'black') +
  geom_point(data = collect.SAM.data, shape = 21,
             aes(x = mean.signal, y = posterior.mean.err,
                 fill = as.character(is.eff),
                 size = num_of_SAMs),
             color = 'black', alpha = 0.75) +
  geom_abline(slope = slope.1, intercept = intercept.1,
              linetype = 'dashed', color = 'blue',
              size = 1) +
  scale_fill_discrete(name = 'Reduced Error:') +
  scale_size_continuous(name = 'Number of SAMs:') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

#lm.fit.2
sounding.density <-
  collect.SAM.data$num_of_soundings/collect.SAM.data$num_of_SAMs

lm.2 <- lm(collect.SAM.data$posterior.mean.err ~ sounding.density)
intercept.2 <- lm.2$coefficients[1]
slope.2 <- lm.2$coefficients[2]
r_val.2 <- sqrt(summary(lm.2)$r.squared)

#regression
density.plot <- ggplot() +
  ggtitle(paste0('Posterior Error Correction vs.', '\n',
                 'Sounding Density')) +
  labs(subtitle = paste0('y=', round(slope.2, 5),
                         'x', if(intercept.2 >= 0) {'+'} else{'-'},
                         round(intercept.2, 2), '; r=',
                         round(r_val.2, 2)),
       caption = paste0('Prior mean error (black dashed line): ',
                        round(prior_diff, 2), 'ppm.')) +
  xlab('Sounding Density') +
  ylab(expression(bar(z[tot] - bold(K)[tot]*hat(lambda))~'[ppm]')) +
  geom_abline(slope = 0, intercept = prior_diff,
              linetype = 'dashed', color = 'black') +
  geom_point(data = collect.SAM.data, shape = 21,
             aes(x = num_of_soundings/num_of_SAMs,
                 y = posterior.mean.err,
                 fill = as.character(is.eff),
                 size = num_of_SAMs),
             color = 'black', alpha = 0.75) +
  geom_abline(slope = slope.2, intercept = intercept.2,
              linetype = 'dashed', color = 'blue',
              size = 1) +
  scale_fill_discrete(name = 'Reduced Error:') +
  scale_size_continuous(name = 'Number of SAMs:') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

multi.regression <-
  lm(collect.SAM.data$posterior.mean.err ~
       collect.SAM.data$mean.signal + sounding.density)
#######################################
#######################################



##############################################
### Collect Effective and Ineffective SAMs ###
##############################################
#bottom five percent
bottom.val <- quantile(collect.SAM.data$posterior.mean.err, 0.05)

#collect iterations of effective SAMs
eff.SAM.list <- subset(collect.SAM.data,
                       posterior.mean.err <= bottom.val)$file.path
SAM.list <- NULL
for(i in 1:length(eff.SAM.list)) {
  iter <- read.csv(eff.SAM.list[i])
  SAM.list <- c(SAM.list, unique(iter$date.time))
}

#top five percent
top.val <- quantile(collect.SAM.data$posterior.mean.err, 0.95)


#collect iterations of ineffective SAMs
ineff.SAM.list <- subset(collect.SAM.data,
                         posterior.mean.err >= top.val)$file.path
SAM.list.2 <- NULL
for(i in 1:length(ineff.SAM.list)) {
  iter <- read.csv(ineff.SAM.list[i])
  SAM.list.2 <- c(SAM.list.2, unique(iter$date.time))
}
##############################################
##############################################



#################
### Frequency ###
#################
#count the frequency of occurrence of (in)effective SAMs
SAM.table <- table(SAM.list)
SAM.table <- SAM.table[order(SAM.table, decreasing = TRUE)]
SAM.table_df <- as.data.frame(SAM.table)

SAM.table.2 <- table(SAM.list.2)
SAM.table.2 <- SAM.table.2[order(SAM.table.2, decreasing = TRUE)]
SAM.table.2_df <- as.data.frame(SAM.table.2)

#match the values from SAM.table.2 to SAM.table
SAM.table_df$Freq.2 <- NA
for(i in 1:nrow(SAM.table_df)) {
  SAM.table_df$Freq.2[i] <-
    SAM.table.2_df$Freq[which(as.character(SAM.table.2_df$SAM.list.2) ==
                                as.character(SAM.table_df$SAM.list[i]))]
}
SAM.table_df$SAM.list <- as.character(SAM.table_df$SAM.list)
SAM.table_df$SAM.list <- factor(SAM.table_df$SAM.list,
                                levels = SAM.table_df$SAM.list)

#Pareto chart of effective SAMs
ggplot() +
  ggtitle('SAMs Contributing to Error Reduction') +
  xlab('SAM') +
  ylab('Count') +
  geom_bar(data = SAM.table_df, stat = 'identity',
           aes(x = SAM.list, y = Freq.2),
           color = 'black', fill = 'black') +
  geom_bar(data = SAM.table_df, stat = 'identity',
           aes(x = SAM.list, y = Freq),
           color = 'black', fill = 'green', alpha = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90))

#effective SAMs
SAM.table_df$SAM.list <- as.POSIXct(SAM.table_df$SAM.list,
                                    tz = local.tz)
eff.SAMs <- xco2_df$date.time %in% SAM.table_df$SAM.list[1:30]
eff.SAMs <- xco2_df[eff.SAMs,]

#ineffective SAMs
SAM.table.2_df$SAM.list.2 <- as.POSIXct(SAM.table.2_df$SAM.list,
                                        tz = local.tz)
ineff.SAMs <- xco2_df$date.time %in% SAM.table.2_df$SAM.list.2[1:30]
ineff.SAMs <- xco2_df[ineff.SAMs,]

#get unique (in)effective SAMs
unique.eff.SAMs <- eff.SAMs$date.time %in% ineff.SAMs$date.time
unique.eff.SAMs <- eff.SAMs[!unique.eff.SAMs,]

unique.ineff.SAMs <- ineff.SAMs$date.time %in% eff.SAMs$date.time
unique.ineff.SAMs <- ineff.SAMs[!unique.ineff.SAMs,]

ineff.plot <- ggmap(ggmap) +
  ggtitle('SAMs Unique to Ineffective Error Reduction') +
  ylab('Ineffective') +
  geom_point(data = unique.ineff.SAMs, shape = 23,
             aes(x = lon, y = lat,
                 fill = xco2- (bio + OCO3.bkg))) +
  facet_wrap(. ~ strftime(date.time, tz = local.tz),
             ncol = 5) +
  scale_fill_gradientn(colors = c('green', 'yellow',
                                  'orange', 'red'),
                       name = 'Enhancement [ppm]',
                       limits = c(0,10)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

eff.plot <- ggmap(ggmap) +
  ylab('Effective') +
  geom_point(data = unique.eff.SAMs, shape = 23,
             aes(x = lon, y = lat,
                 fill = xco2- (bio + OCO3.bkg))) +
  facet_wrap(. ~ strftime(date.time, tz = local.tz),
             ncol = 5) +
  scale_fill_gradientn(colors = c('green', 'yellow',
                                  'orange', 'red'),
                       name = 'Enhancement [ppm]',
                       limits = c(0,10)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1,
                                frame.colour = "black",
                                frame.linewidth = 1)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = 'bottom',
        legend.key.width = unit(0.9, 'in'))

#create the combo plot
combo.plot <- ineff.plot / eff.plot
combo.plot <- combo.plot +
  plot_annotation(title = 'SAMs Unique to (In)effective Error Reduction',
  theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave(combo.plot, filename = 'Out/Results/Unique_SAMs.jpg',
       device = 'jpg', height = 9, width = 7, units = 'in')
####################################
####################################



eff.SAMs.agg_df <- data.frame(matrix(NA, nrow = 0, ncol = 7))
names(eff.SAMs.agg_df) <- c('date.time', 'soundings', 'bkg.soundings',
                            'mean.signal', 'Moran.i', 'mean.dist',
                            'in.domain')
ineff.SAMs.agg_df <- data.frame(matrix(NA, nrow = 0, ncol = 7))
names(ineff.SAMs.agg_df) <- names(eff.SAMs.agg_df)
for(i in 1:10) {
  
  #obtain each effective SAM
  eff.SAM <- subset(unique.eff.SAMs,
                    date.time == unique(unique.eff.SAMs$date.time)[i])
  
  #over domain
  in.domain <- subset(eff.SAM,
                      (lon >= xmin & lon <= xmax) &
                      (lat >= ymin & lat <= ymax))
  perc_in_dom <- 100*nrow(in.domain)/nrow(eff.SAM)
  
  #calculate background
  eff.SAM.bkg <- subset(eff.SAM,
                        rowSums(eff.SAM[,gsub(' ', '', priors[,1])]) <
                        OCO3.bkg.threshold)
  eff.num_bkg <- nrow(eff.SAM.bkg)
  
  #Moran's i
  xco2.dists <- as.matrix(dist(cbind(eff.SAM$lon, eff.SAM$lat)))
  xco2.dists.inv <- 1/xco2.dists
  diag(xco2.dists.inv) <- 0
  M.i <- Moran.I(eff.SAM$xco2, xco2.dists.inv)
  
  #calculate density
  dist.vals <- as.matrix(dist(cbind(eff.SAM$lon, eff.SAM$lat)))
  dist.vals <- unique(as.vector(dist.vals))
  
  #build add line
  eff.add.line <- data.frame(date.time = unique(eff.SAMs$date.time)[i],
                             soundings = nrow(eff.SAM),
                             bkg.soundings = eff.num_bkg,
                             mean.signal = mean(eff.SAM$xco2 -
                                               (eff.SAM$bio +
                                                eff.SAM$OCO3.bkg)),
                             Moran.i = as.numeric(M.i[1]),
                             mean.dist = mean(dist.vals),
                             in.domain = perc_in_dom)
  eff.SAMs.agg_df <- rbind(eff.SAMs.agg_df, eff.add.line)
  
  #obtain each effective SAM
  ineff.SAM <- subset(unique.ineff.SAMs,
                      date.time == unique(unique.ineff.SAMs$date.time)[i])
  
  #over domain
  in.domain <- subset(ineff.SAM,
                      (lon >= xmin & lon <= xmax) &
                        (lat >= ymin & lat <= ymax))
  perc_in_dom <- 100*nrow(in.domain)/nrow(ineff.SAM)
  
  #calculate background
  ineff.SAM.bkg <- subset(ineff.SAM,
                          rowSums(ineff.SAM[,gsub(' ', '', priors[,1])]) <
                          OCO3.bkg.threshold)
  ineff.num_bkg <- nrow(ineff.SAM.bkg)
  
  #Moran's i
  xco2.dists <- as.matrix(dist(cbind(ineff.SAM$lon, ineff.SAM$lat)))
  xco2.dists.inv <- 1/xco2.dists
  diag(xco2.dists.inv) <- 0
  M.i <- Moran.I(ineff.SAM$xco2, xco2.dists.inv)
  
  #calculate density
  dist.vals <- as.matrix(dist(cbind(ineff.SAM$lon, ineff.SAM$lat)))
  dist.vals <- unique(as.vector(dist.vals))
  
  ineff.add.line <- data.frame(date.time = unique(ineff.SAMs$date.time)[i],
                               soundings = nrow(ineff.SAM),
                               bkg.soundings = ineff.num_bkg,
                               mean.signal = mean(ineff.SAM$xco2 -
                                                 (ineff.SAM$bio +
                                                  ineff.SAM$OCO3.bkg)),
                               Moran.i = as.numeric(M.i[1]),
                               mean.dist = mean(dist.vals),
                               in.domain = perc_in_dom)
  ineff.SAMs.agg_df <- rbind(ineff.SAMs.agg_df, ineff.add.line)
  
}

combo_df <- rbind(data.frame(eff.SAMs.agg_df,
                             Inversion = 'Effective'),
                  data.frame(ineff.SAMs.agg_df,
                             Inversion = 'Not Effective'))

#testing Moran's i
test_1 <- wilcox.test(eff.SAMs.agg_df$Moran.i,
                      ineff.SAMs.agg_df$Moran.i,
                      'greater')
test_1.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_1$p.value, 3))) +
  ylab("Moran's i") +
  geom_boxplot(data = combo_df,
               aes(x = Inversion, y = Moran.i)) +
  geom_point(data = combo_df,
             aes(x = Inversion, y = Moran.i)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#testing soundings
test_2 <- wilcox.test(eff.SAMs.agg_df$soundings,
                      ineff.SAMs.agg_df$soundings,
                      'greater')
test_2.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_2$p.value, 3))) +
  ylab('Number of Soundings') +
  geom_boxplot(data = combo_df,
               aes(x = Inversion, y = soundings)) +
  geom_point(data = combo_df,
             aes(x = Inversion, y = soundings)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#testing number of background soundings
test_3 <- wilcox.test(eff.SAMs.agg_df$bkg.soundings,
                      ineff.SAMs.agg_df$bkg.soundings,
                      alternative = 'less')
test_3.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_3$p.value, 3))) +
  ylab(paste0('Number of' , '\n', 'Background Soundings')) +
  geom_boxplot(data = combo_df,
               aes(x = Inversion, y = bkg.soundings)) +
  geom_point(data = combo_df,
             aes(x = Inversion, y = bkg.soundings)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#testing ratio of background soundings
test_4 <- wilcox.test(100*eff.SAMs.agg_df$bkg.soundings/
                        eff.SAMs.agg_df$soundings,
                      100*ineff.SAMs.agg_df$bkg.soundings/
                        ineff.SAMs.agg_df$soundings,
                      alternative = 'less')
test_4.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_4$p.value, 3))) +
  ylab(paste0(paste0('Percent of ', '\n',
                     'Background Soundings [%]'))) +
  geom_boxplot(data = combo_df,
               aes(x = Inversion,
                   y = 100*bkg.soundings/soundings)) +
  geom_point(data = combo_df,
             aes(x = Inversion,
                 y = 100*bkg.soundings/soundings)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#testing mean signal
test_5 <- wilcox.test(eff.SAMs.agg_df$mean.signal,
                      ineff.SAMs.agg_df$mean.signal,
                      'greater')
test_5.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_5$p.value, 3))) +
  ylab('Mean Enhancement [ppm]') +
  geom_boxplot(data = combo_df,
               aes(x = Inversion, y = mean.signal)) +
  geom_point(data = combo_df,
             aes(x = Inversion, y = mean.signal)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#testing domain coverage
test_6 <- wilcox.test(eff.SAMs.agg_df$in.domain,
                      ineff.SAMs.agg_df$in.domain,
                      'greater')
test_6.plot <- ggplot() +
  labs(subtitle = paste0('p=', round(test_6$p.value, 3))) +
  ylab('Domain Coverage [%]') +
  geom_boxplot(data = combo_df,
               aes(x = Inversion, y = in.domain)) +
  geom_point(data = combo_df,
             aes(x = Inversion, y = in.domain)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

t.test_plot <- (test_1.plot + test_2.plot + test_3.plot)/
  (test_4.plot + test_5.plot + test_6.plot)
t.test_plot <- t.test_plot +
  plot_annotation(title = 'T-test Results of SAMs',
                  theme = 
                    theme(plot.title = element_text(hjust = 0.5)))
ggsave(t.test_plot, filename = 'Out/Results/t_tests.jpg',
       device = 'jpg', height = 5, width = 7, unit = 'in')


###################
### Filter SAMs ###
###################
m.intercept <- multi.regression$coefficients[1]
m.slope.1 <- multi.regression$coefficients[2]
m.slope.2 <- multi.regression$coefficients[3]

#equation
linear.filter <- function(a, x, b, y, c) {
  z = a*x + b*y + c
  return(z)
}

SAM.times <- aggregate((xco2 - (bio + OCO3.bkg)) ~ date.time,
                        data = xco2_df, mean)
SAM.times$num_of_soundings <- aggregate(xco2 ~ date.time,
                                        data = xco2_df, NROW)[,2]
names(SAM.times) <- c('date.time', 'mean.xco2', 'num_of_soundings')

#outcomes from linear filtering
SAM.times$outcomes <- linear.filter(m.slope.1, SAM.times$mean.xco2,
                                    m.slope.2, SAM.times$num_of_soundings,
                                    m.intercept)

new.SAMs <-
  SAM.times$date.time[which(SAM.times$num_of_soundings >= 275)]
new.soundings <- xco2_df$date.time %in% new.SAMs
new.soundings <- xco2_df[new.soundings,]

#make R
xco2.errors_df <-
  #observed
  (new.soundings$xco2 - (new.soundings$bio +
                           new.soundings$OCO3.bkg)) -
  #modeled
  (rowSums(new.soundings[,gsub(' ', '', priors[,1])]))

spatial.error_df <-
  data.frame(time = new.soundings$date.time,
             lon = new.soundings$lon,
             lat = new.soundings$lat,
             error = abs(xco2.errors_df))

R <- spatial.correlation(spatial.error = spatial.error_df,
                         vgm.binwidth = 3, vgm.cutoff = 20,
                         included.xco2.errors =
                           'ext/included.errors.csv',
                         plot.output.path = 'Out/new.R/')

########################
### Do New Inversion ###
########################
K <- rowSums(new.soundings[,gsub(' ', '', priors[,1])])
z <- new.soundings$xco2 -
  (new.soundings$bio + new.soundings$OCO3.bkg)

term.1 <- as.matrix(S_p) %*% t(K)
term.2 <- solve((K %*% as.matrix(S_p) %*% t(K)) + R)
term.3 <- (z - (K %*% as.matrix(lambda_p)))

lambda_hat <- lambda_p + (term.1 %*% (term.2 %*% term.3))
S_error <- solve(t(K) %*% solve(R) %*% K + solve(S_p))

write.csv(new.SAMs, file = paste0('ext/filtered_SAMs.csv'),
          row.names = FALSE)
