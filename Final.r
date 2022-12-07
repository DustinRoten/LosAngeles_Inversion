library(ggplot2); library(ggmap); library(scales); library(patchwork)
library(lubridate); library(stringr)

home.dir <- '/uufs/chpc.utah.edu/common/home'
work.dir <- 'u1211790/LosAngeles_Inversion'
data.dir <- 'lin-group14/DDR/Roten_InputData'
output.dir <- 'Out/Results'
setwd(file.path(home.dir, work.dir))

#time filter
base.time <- '201501010000-201512311159'
target.time <- '202001010000-202106301159'

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

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

#data paths
traffic.path <- file.path(home.dir, data.dir, 'PeMS')
aircraft.path <- file.path(home.dir, data.dir, 'LAWA')
manufacturing.path <- file.path(home.dir, data.dir, 'GHGRP')
power.path <- -3.06

#read in data
# power industry, manufacturing, on-road, non-road, buildings
# lambda_hat <- c(0.97, 0.82, 1.07, 1.25, 1.09)
# lambda_err <- c(0.02, 0.01, 0.01, 0.03, 0.03)

#datasets
traffic.data <- read.csv(file.path(traffic.path,
                                   'interstate_data.csv'))
traffic.scaling.factors <- read.csv(file.path(home.dir,
                                              work.dir,
                                              'ext',
                                              'traffic_co2.csv'))
aircraft.data <- read.csv(file.path(aircraft.path,
                                    'LAX.csv'))
manufacturing.data <- read.csv(file.path(manufacturing.path,
                                         'GHGRP_Comparison.csv'))
####################
### Traffic Data ###
####################
#scale the traffic data based on CO2/mi
traffic.years <- unique(year(traffic.data$Hour))
for(i in 1:length(traffic.years)) {
  
  #match years and scale them to g/mi CO2
  traffic.data$VMT[year(traffic.data$Hour) == traffic.years[i]] <-
    traffic.data$VMT[year(traffic.data$Hour) == traffic.years[i]]*
    traffic.scaling.factors$Year[traffic.scaling.factors$Year ==
                                   traffic.years[i]]
  
}

#subset by prior and posterior data
traffic.data$time.period <- NA

#convert time to posix
traffic.data$Hour <- as.POSIXct(traffic.data$Hour,
                                format = '%Y-%m-%d %H:%M:%S',
                                tz = 'UTC')

traffic.data$time.period[traffic.data$Hour >= base.time[1] &
                         traffic.data$Hour <= base.time[2]] <- 'Prior'
traffic.data$time.period[traffic.data$Hour >= target.time[1] &
                           traffic.data$Hour <= target.time[2]] <- 'Posterior'
traffic.data <- subset(traffic.data, !is.na(time.period))

#finish converting to mtCO2/mi
traffic.data$VMT <- ((traffic.data$VMT/1000)/1000)
names(traffic.data) <- c('Hour', 'co2', 'VHT', 'Interstate')

#subset time periods
sd.err <- function(x) {sd(x, na.rm = TRUE)/length(x)}
agg.traffic.data <-
  aggregate(co2 ~ year(Hour) + Interstate,
            data = traffic.data, mean)
agg.traffic.data$error <-
  aggregate(co2 ~ year(Hour) + Interstate,
            data = traffic.data, sd.err)$co2
names(agg.traffic.data) <- c('Year', 'Interstate', 'co2', 'error')

#create a new dataframe to store 
rel.traffic.data <- agg.traffic.data
for(i in 1:nrow(agg.traffic.data)) {
  
  agg.traffic.data_line <- agg.traffic.data[i,]
  agg.base.year.value <-
    subset(agg.traffic.data,
           Year == agg.base.year &
             Interstate == agg.traffic.data_line$Interstate)
  
  #add the relative CO2 value
  rel.traffic.data$co2[i] <-
    agg.traffic.data_line$co2/agg.base.year.value$co2
  
  #add the associated error
  rel.traffic.data$error[i] <-
    100*sqrt(agg.traffic.data_line$error^2 +
               agg.base.year.value$error^2)
  
}

#one final aggregation
traffic.plot.vals <- subset(rel.traffic.data, Year != agg.base.year)

#a small function for adding errors in quadrature
in.quadrature <- function(x) {sqrt(sum(x^2))}

tmp <- aggregate(co2 ~ Interstate, data = traffic.plot.vals, mean)
tmp$error <- aggregate(error ~ Interstate, data = traffic.plot.vals,
                       in.quadrature)$error
traffic.plot.vals <- tmp
traffic.plot.vals$co2 <- 100*(traffic.plot.vals$co2 - 1)
names(traffic.plot.vals) <- c('Interstate', 'Perc.Change', 'Error')
####################
####################
####################


##############################
### Nonroad Transportation ###
##############################
passenger.plot <- ggplot() +
  ggtitle('Passenger') +
  xlab('Month') +
  ylab('Passengers (in millions)') +
  geom_line(data = aircraft.data,
            aes(x = as.integer(Month),
                y = (Domestic.Passengers + International.Passengers)/1e6,
                color = as.character(Year),
                group = as.character(Year))) +
  geom_point(data = aircraft.data,
             aes(x = as.integer(Month),
                 y = (Domestic.Passengers + International.Passengers)/1e6,
                 color = as.character(Year),
                 group = as.character(Year))) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_discrete(name = 'Year') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

cargo.plot <- ggplot() +
  ggtitle('Cargo') +
  xlab('Month') +
  ylab('Tons of Cargo (in thousands)') +
  geom_line(data = aircraft.data,
            aes(x = Month,
                y = Cargo..Tons./1e3,
                color = as.character(Year),
                group = as.character(Year))) +
  geom_point(data = aircraft.data,
             aes(x = Month,
                 y = Cargo..Tons./1e3,
                 color = as.character(Year),
                 group = as.character(Year))) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_discrete(name = 'Year') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

aircraft.plots <- passenger.plot | cargo.plot
aircraft.plots <- aircraft.plots +
  plot_annotation(title = 'Aircraft Transport (LAX)',
                  theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave(aircraft.plots,
       filename = file.path(output.dir, 'Aircraft_Data.jpg'),
       device = 'jpg', unit = 'in',
       height = 4, width = 6)

#selected aircraft data
agg.aircraft.data <- aggregate(Domestic.Passengers +
                               International.Passengers ~ Year,
                               data = aircraft.data, sum)

agg.aircraft.data$Cargo <- aggregate(Cargo..Tons. ~ Year,
                                     data = aircraft.data, sum)$Cargo..Tons.
names(agg.aircraft.data) <- c('Year', 'Passengers', 'Cargo')

#create a new dataframe to store relative aircraft information
rel.aircraft.data <- agg.aircraft.data
for(i in 1:nrow(agg.aircraft.data)) {
  
  agg.aircraft.data_line <- agg.aircraft.data[i,]
  agg.base.year.value <- subset(agg.aircraft.data, Year == agg.base.year)
  
  #add the relative passenger value
  rel.aircraft.data$Passengers[i] <-
    agg.aircraft.data_line$Passengers/agg.base.year.value$Passengers
  
  #add the relative cargo value
  rel.aircraft.data$Cargo[i] <-
    agg.aircraft.data_line$Cargo/agg.base.year.value$Cargo
  
}

#add final values for plotting
aircraft.plot.vals <- subset(rel.aircraft.data, Year != agg.base.year)
aircraft.passengers <- 100*(mean(aircraft.plot.vals$Passengers) - 1)
aircraft.cargo <- 100*(mean(aircraft.plot.vals$Cargo) - 1)
aircraft.plot.vals <- data.frame(Passengers = aircraft.passengers,
                                 Cargo = aircraft.cargo)
##############################
##############################
##############################



#####################
### Manufacturing ###
#####################
#calculate the mean percent change
mean.perc.manufacturing <- round(mean(manufacturing.data$perc.diff), 2)
mean.perc.manufacturing <- as.character(mean.perc.manufacturing)

#generate spatial plot
spatial.manufacturing <- ggmap(map) +
  ggtitle(paste0('Percent Change: 2015 to 2020', '\n',
                 'Manufacturing Emissions')) +
  labs(subtitle = as.expression(bquote(bar(Delta)['%'] == .(mean.perc.manufacturing)*'%'))) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_point(data = manufacturing.data, shape = 21,
             aes(x = lon, y = lat, fill = perc.diff,
                 size = co2.2020)) +
  scale_fill_gradient2(name = 'Percent Change [%]',
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 0, trans = pseudo_log_trans(),
                       limits = c(-100, 100),
                       breaks = c(-100, -10, 0, 10, 100)) +
  scale_size_continuous(name = expression(paste('2020 CO'[2], ' Emissions [tons]')),
                        trans = 'log10') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.direction = 'horizontal', legend.box = 'vertical',
        legend.key.width = unit(0.75, 'in'))
ggsave(spatial.manufacturing,
       filename = file.path(output.dir, 'Manufacturing_Data.jpg'),
       device = 'jpg', unit = 'in',
       height = 9, width = 9)
#####################
#####################
#####################



########################
### Plot All Results ###
########################
converted.lambdas <- lambda_hat
converted.lambdas <- 100*(converted.lambdas-1)
converted.errors <- 100*lambda_err

#lambda values
lambda_df <- data.frame(value = converted.lambdas,
                        error = converted.errors,
                        source = 'OCO-3',
                        category = c('Power Industry',
                                     'Manufacturing',
                                     'On-Road',
                                     'Non-Road',
                                     'Buildiings'))

#traffic data
traffic_df <- data.frame(value = traffic.plot.vals$Perc.Change,
                         error = traffic.plot.vals$Error,
                         source = traffic.plot.vals$Interstate,
                         category = 'On-Road')

#aircraft data
aircraft_df.1 <- data.frame(value = aircraft.plot.vals$Passengers,
                            error = NA,
                            source = 'Aircraft Passenger',
                            category = 'Non-Road')
aircraft_df.2 <- data.frame(value = aircraft.plot.vals$Cargo,
                            error = NA,
                            source = 'Aircraft Cargo',
                            category = 'Non-Road')
aircraft_df.3 <- data.frame(value = aircraft.plot.vals$Cargo,
                            error = NA,
                            source = 'Aircraft Cargo',
                            category = 'On-Road')

#manufacturing data
manufacturing_df <- data.frame(value = as.numeric(mean.perc.manufacturing),
                               error = NA,
                               source = 'GHGRP',
                               category = 'Manufacturing')

combined_df <- rbind(lambda_df,
                     traffic_df,
                     aircraft_df.1,
                     aircraft_df.2,
                     aircraft_df.3,
                     manufacturing_df)
combined_df$source <- factor(combined_df$source,
                             levels = c('OCO-3',
                                        'I10', 'I105', 'I210', 'I405', 'I605',
                                        'Aircraft Passenger', 'Aircraft Cargo',
                                        'GHGRP'))

#temporary - add power plants
add.line <- data.frame(value = 100-3.06,
                       error = NA,
                       source = 'eGRID',
                       category = 'Power Industry')
combined_df <- rbind(combined_df, add.line)

#testing
combined_df$value <- combined_df$value + 100

all.plot.data <- ggplot(data = combined_df,
                        aes(x = category, y = value, group = source)) +
  ggtitle('Posterior Scaling Factors') +
  xlab('Emission Categories') +
  ylab('[%]') +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  geom_col(aes(fill = source),
           color = 'black', position = position_dodge(width = 0.8),
           width = 0.8) +
  geom_text(aes(label = source,
                y = ifelse(value <= 0, 1, value - 1),
                hjust = ifelse(value >= 0, 1, 0)),
            angle = 90, vjust = 0.5,
            position = position_dodge(width = 0.8)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
ggsave(all.plot.data,
       filename = file.path(output.dir, 'All_Data.jpg'),
       device = 'jpg', unit = 'in',
       height = 4, width = 6)

#quick melt data for plotting
melt_df <- melt(lambda_hat_df)
melt_df$Iteration <- rep(1:10, 7)
melt_df$variable <- as.character(melt_df$variable)
melt_df <-
  melt_df[!(melt_df$variable %in% c('Iteration', 'mean.error')),]


lambda_iterations <- ggplot() +
  ggtitle('Iterative Scaling Factors') +
  xlab('Iteration') +
  ylab(expression(hat(lambda))) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_line(data = melt_df,
            aes(x = Iteration, y = value,
                group = variable, color = variable)) +
  geom_point(data = melt_df,
             aes(x = Iteration, y = value,
                 group = variable, color = variable)) +
  scale_color_discrete(name = 'Category') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
ggsave(lambda_iterations,
       filename = 'Out/Results/Lambda_Iterations.jpg',
       device = 'jpg', unit = 'in',
       height = 4, width = 6.5)
