#########################
### Plot Traffic Data ###
#########################
#load necessary packages
library(lubridate); library(ggplot2); library(scales)

#setup the input file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/Roten_InputData'

setwd(file.path(home.dir, work.ext))
dir.create('Out', showWarnings = FALSE)

base.year <- 2015
tz <- 'America/Los_Angeles'

#data breaks (UTC)
data.breaks <- c('20200101', '20210630')

#construct the data directory
traffic.data.dir <- file.path(home.dir, data.ext, 'PeMS')
traffic.data <- read.csv(file.path(traffic.data.dir,
                                   'interstate_data.csv'))

#POSIX time format
traffic.data$Hour <- as.POSIXct(traffic.data$Hour,
                                format = '%Y-%m-%d %H:%M:%S',
                                tz = tz)
#add week number
traffic.data$Week <- as.numeric(strftime(traffic.data$Hour,
                                         format = '%W',
                                         tz = tz))

#add day number
traffic.data$Day <- wday(traffic.data$Hour)

#read in the traffic co2 data
traffic.co2 <- read.csv('ext/traffic_co2.csv')
year.list <- unique(year(traffic.data$Hour))

traffic.data$co2 <- NA
for(yr in year.list) {
  
  which.yr <- which(year(traffic.data$Hour) == yr)
  
  scaling.factor <- traffic.co2$grams_per_mi[traffic.co2$Year == yr]
  
  traffic.data$co2[which.yr] <-
    (10^-6)*(traffic.data$VMT[which.yr]*scaling.factor)
  
}

base.traffic.data <- subset(traffic.data, year(Hour) == base.year)
sub.traffic.data <- subset(traffic.data, year(Hour) != base.year)

#set up matrix for relative data
relative.data <- data.frame(matrix(NA, nrow = 0, ncol = 3))
names(relative.data) <- c('Hour', 'Relative', 'Interstate')
for(i in 1:nrow(sub.traffic.data)) {
  
  #grab the corresponding entry from the base data
  base.line <- subset(base.traffic.data,
                      Interstate == sub.traffic.data$Interstate[i] &
                      Week == sub.traffic.data$Week[i] &
                      Day == sub.traffic.data$Day[i] &
                      hour(Hour) == hour(sub.traffic.data$Hour[i]))
  
  if(nrow(base.line) != 0) {
    #grab the entry of interest
    add.line <- data.frame(Hour =
                             sub.traffic.data$Hour[i],
                           Relative =
                             100*(sub.traffic.data$co2[i]/base.line$co2 - 1),
                           Interstate =
                             sub.traffic.data$Interstate[i])
    relative.data <- rbind(relative.data, add.line)
  }
  
  #message to user
  msg <- paste0(round(100*i/nrow(sub.traffic.data), 2),
                '% complete.     ', '\r')
  if(!exists('old.msg')) {old.msg <- 'none'}
  if(old.msg != msg) {
    cat(msg)
    old.msg <- msg
  }
  
}; cat(paste('\n'))

#convert the data breaks into
data.breaks <- as.POSIXct(data.breaks, format = '%Y%m%d',
                          tz = 'UTC')
attr(data.breaks, 'tzone') <- tz

#divide the data into periods
relative.data$Period <- NA
for(i in 1:(length(data.breaks) - 1)) {
  
  #select the values within the timeframe
  idx.values <- which(relative.data$Hour >= data.breaks[i] &
                        relative.data$Hour < data.breaks[i+1])
  
  #add the period number to these entries
  relative.data$Period[idx.values] <- i
  
}

aggregated.data <- aggregate(Relative ~ Interstate + Period,
                             data = relative.data, mean)
ggplot() +
  geom_bar(data = aggregated.data, stat = 'identity',
           aes(x = Interstate, y = Relative,
               fill = Interstate), color = 'black')
