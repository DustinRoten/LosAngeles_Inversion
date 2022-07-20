library(lubridate); library(ggplot2)

home.path <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
epa.ext <- 'lin-group14/DDR/Roten_InputData/EPA'
data.tz <- 'UTC'
#indicate which year to use as background
base.year <- 2015
target.years <- c(2019, 2020, 2021)

#data breaks (UTC)
data.breaks <- c('20190814', '20200401', '20200505', '20201231')
base.breaks <- c('20140814', '20150401', '20150505', '20151231')

#construct the epa.path
epa.path <- file.path(home.path, epa.ext)

#list the available PeMS files
epa.files <- list.files(epa.path, full.names = TRUE,
                        pattern = '\\.csv$')

#get base year
base.year_df <- grep(epa.files, pattern = paste0(base.year),
                     value = TRUE)
base.year_df <- read.csv(base.year_df)[,1:10]
names(base.year_df) <- c('State', 'Name', 'ID', 'Year', 'Date',
                         'Hour', 'SO2_lbs', 'NOX_lbs', 'CO2_ShortTons',
                         'Heat_MMBTU')

#get target year
target.years_df <- grep(epa.files, pattern = paste0(base.year),
                        invert = TRUE, value = TRUE)
target.years_df <- read.csv(target.years_df)[,1:10]
names(target.years_df) <- names(base.year_df)

#add a POSIX column
base.year_df$date.time <- as.POSIXct(paste0(base.year_df$Date, ' ',
                                            base.year_df$Hour),
                                     format = '%Y-%m-%d %H',
                                     origin = '1970-01-01',
                                     tz = 'UTC')
target.years_df$date.time <- as.POSIXct(paste0(target.years_df$Date, ' ',
                                               target.years_df$Hour),
                                        format = '%Y-%m-%d %H',
                                        origin = '1970-01-01',
                                        tz = 'UTC')

#combine all of the available data
all.years_df <- rbind(base.year_df, target.years_df)

#add the Week number
base.year_df$Week <- strftime(base.year_df$date.time,
                              format = '%W')
target.years_df$Week <- strftime(target.years_df$date.time,
                                 format = '%W')
#add the Day number
base.year_df$Day <- wday(base.year_df$date.time)
target.years_df$Day <- wday(target.years_df$date.time)

#modify breaks to subset the base year
data.breaks <- as.POSIXct(data.breaks, format = '%Y%m%d',
                          tz = data.tz)
base.breaks <- as.POSIXct(base.breaks, format = '%Y%m%d',
                          tz = data.tz)

#divide the base year data into periods
base.year_df$Period <- NA
for(i in 1:(length(base.breaks) - 1)) {
  #select the values within the timeframe
  idx.values <- which(base.year_df$date.time >= base.breaks[i] &
                        base.year_df$date.time < base.breaks[i+1])
  #add the period number to these entries
  base.year_df$Period[idx.values] <- i
}; base.year_df <- subset(base.year_df,
                          !is.na(CO2_ShortTons) & !is.na(Period))

#divide the base year data into periods
target.years_df$Period <- NA
for(i in 1:(length(data.breaks) - 1)) {
  #select the values within the timeframe
  idx.values <- which(target.years_df$date.time >= data.breaks[i] &
                        target.years_df$date.time < data.breaks[i+1])
  #add the period number to these entries
  target.years_df$Period[idx.values] <- i
}; target.years_df <- subset(target.years_df,
                             !is.na(CO2_ShortTons) & !is.na(Period))

#build a standard error function
sd.err <- function(x) {sd(x)/sqrt(length(x))}

#aggregate the base data
agg.base <- aggregate(CO2_ShortTons ~ Name + Period,
                      data = base.year_df, mean)
agg.base$sd.err <- aggregate(CO2_ShortTons ~ Name + Period,
                             data = base.year_df, sd.err)[,3]

#aggregate the target data
agg.target <- aggregate(CO2_ShortTons ~ Name + Period,
                        data = target.years_df, mean)
agg.target$sd.err <- aggregate(CO2_ShortTons ~ Name + Period,
                               data = target.years_df, sd.err)[,3]

#construct a dataframe of relative values
agg.relative <- data.frame(agg.target[,1:2],
                           Relative = 100*(agg.target[,3]/agg.base[,3] - 1))
agg.relative$sd.err <- sqrt(agg.base$sd.err^2 + agg.target$sd.err^2)

all.plants <- aggregate(Relative ~ Period,
                        data = agg.relative, mean)

ggplot() +
  ggtitle(expression(paste('Percent Change in CO'[2], ' Emissions'))) +
  labs(subtitle = 'Power Plants') +
  ylab('Percent Difference [%]') +
  xlab('Time Period') +
  geom_hline(yintercept = 0) +
  geom_bar(data = agg.relative, stat = 'identity',
           aes(x = Period, y = Relative,
               fill = Name),
           color = 'black',
           position = 'dodge') +
  # geom_errorbar(data = agg.relative,
  #               aes(x = Period,
  #                   ymin = Relative - sd.err,
  #                   ymax = Relative + sd.err),
  #               position = position_dodge(0.9)) +
  geom_bar(data = all.plants, stat = 'identity',
           aes(x = Period, y = Relative),
           color = 'black', fill = 'white',
           linetype = 'dashed', alpha = 0.75) +
  scale_fill_discrete(name = 'Plant Name:') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
