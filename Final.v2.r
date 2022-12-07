library(ggplot2); library(ggmap); library(scales); library(patchwork)
library(lubridate); library(stringr)

home.dir <- '/uufs/chpc.utah.edu/common/home'
work.dir <- 'u1211790/LosAngeles_Inversion'
data.dir <- 'lin-group14/DDR/Roten_InputData'
output.dir <- 'Out/Results'
setwd(file.path(home.dir, work.dir))
source('r/dependencies.r')

#time filter
base.time <- '201501010000-201512311159'
target.time <- '202001010000-202106301159'

#base year
comparison.year <- 2015
tzone <- 'America/Los_Angeles'

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

#obtain map
map <- get_map(location = 'Los Angeles', maptype = 'terrain',
               color = 'bw', zoom = 10)

##########
#convert times
agg.base.year <- 2015
base.time <- str_split_fixed(base.time, pattern = '-', n = 2)
base.time <- as.POSIXct(base.time, format = '%Y%m%d%H%M',
                        tz = 'UTC')
target.time <- str_split_fixed(target.time, pattern = '-', n = 2)
target.time <- as.POSIXct(target.time, format = '%Y%m%d%H%M',
                          tz = 'UTC')

### Identify all of the input data paths ###
#inversion results
bayesian.path <- file.path(home.dir, work.dir, 'Out/Results')
bayesian.data <- format_bayesian.data(bayesian.path, tzone = tzone)

#traffic
traffic.path <- file.path(home.dir, data.dir, 'PeMS')
scaling.factors.path <- file.path(home.dir, work.dir, 'ext')

traffic.data <-
  format_traffic.data(traffic.path, scaling.factors.path,
                      tzone = tzone, comparison.year = comparison.year)
traffic.data <- subset(traffic.data, year(time) > 2019)

daily.traffic.data <- format_traffic.data.daily(traffic.data, tzone)
daily.traffic.data <- subset(daily.traffic.data, year(time) > 2019)

traffic.data$Sector <- 'OnRoad'
daily.traffic.data$Sector <- 'OnRoad'

#port
port.path <- file.path(home.dir, data.dir, 'Ports')
port.data <- read.csv(list.files(port.path, pattern = 'Port_Data.csv',
                                 full.names = TRUE))
port.data <- format_port.data(port.data, comparison.year)
port.data$Sector <- 'Marine'

#aircraft
aircraft.path <- file.path(home.dir, data.dir, 'LAWA')
aircraft.data <- format_aircraft.data(aircraft.path, comparison.year)
aircraft.data$Sector <- 'Aviation'

#manufacturing
manual.manufacturing <- data.frame(Value = -65, Sector = 'Industry') # (%) from GHGRP 
manufacturing.path <- file.path(home.dir, data.dir, 'GHGRP')
manufacturing.data <- read.csv(file.path(manufacturing.path,
                                         'GHGRP_Comparison.csv'))
mean.manufacturing.data <-
  data.frame(Value = mean(manufacturing.data$perc.diff),
             Sector = 'Industry')

#power plants
powerplant.path <- file.path(home.dir, data.dir, 'eGRID')
manual.powerplant <- data.frame(Value = -5, Sector = 'Power Industry') # (%) from eGRID

### Plot the Data ###
all.sectors.plots <- ggplot() +
  
  ggtitle('Optimized Vulcan 3.0 Sectors and Proxy Data') +
  xlab('Time') +
  ylab('Percent Change [%]') +
  geom_hline(yintercept = 0, linetype = 'dotted',
             color = 'gray') +
  
  #add traffic
  # geom_line(data = traffic.data,
  #           aes(x = time, y = perc_change)) +
  geom_line(data = daily.traffic.data, color = 'red',
            aes(x = time, y = perc_change)) +
  
  #add ports
  geom_line(data = port.data, linetype = 'dashed',
            aes(x = timestamp, y = Rel.TEU,
                color = Port, group = Port)) +
  
  #add aircraft
  geom_line(data = aircraft.data,
            color = 'blue', linetype = 'dashed',
            aes(x = Time, y = Domestic)) +
  geom_line(data = aircraft.data,
            color = 'green', linetype = 'dashed',
            aes(x = Time, International)) +
  geom_line(data = aircraft.data,
            color = 'red', linetype = 'dashed',
            aes(x = Time, Cargo)) +
  
  #add manufacturing
  geom_hline(data = mean.manufacturing.data,
             aes(yintercept = Value),
             color = 'blue', linetype = 'dashed') +
  geom_hline(data = manual.manufacturing,
             aes(yintercept = Value),
             color = 'red', linetype = 'dashed') +
  
  #add power plant
  geom_hline(data = manual.powerplant,
             aes(yintercept = Value),
             color = 'red', linetype = 'dashed') +
  
  geom_ribbon(data = bayesian.data, color = 'gray',
              alpha = 0.25,
              aes(x = time,
                  ymin = ScalingFactor - 100*error,
                  ymax = ScalingFactor + 100*error)) +
  
  geom_line(data = bayesian.data,
            aes(x = time, y = ScalingFactor)) +
  
  geom_point(data = bayesian.data, color = 'black',
             size = 0.75,
             aes(x = time, y = ScalingFactor)) +
  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  
  facet_wrap(. ~ Sector, ncol = 2, scales = 'free')
ggsave(all.sectors.plots, filename = 'Out/Results/All_Sectors.jpg',
       device = 'jpg', height = 8, width = 8, units = 'in')

#Covariance plots
cov.data.power <- subset(bayesian.data, Sector == 'Power Industry')
cov.data.industry <- subset(bayesian.data, Sector == 'Industry')
cov.data.road <- subset(bayesian.data, Sector == 'OnRoad')
cov.data.marine <- subset(bayesian.data, Sector == 'Marine')

ggplot() +
  geom_abline(slope = 1, intercept = 0, color = 'black',
            linetype = 'dashed') +
  geom_point(aes(x = cov.data.road$ScalingFactor,
                 y = cov.data.power$ScalingFactor))

ggplot() +
  geom_abline(slope = 1, intercept = 0, color = 'black',
              linetype = 'dashed') +
  geom_point(aes(x = cov.data.power$ScalingFactor,
                 y = cov.data.industry$ScalingFactor))
