library(ggplot2); library(lubridate)

TCCON.data <- read.csv(file.path('/uufs/chpc.utah.edu/common/home',
                                 'lin-group14/DDR/Roten_InputData/TCCON/Caltech',
                                 'ci20120920_20201229.public.csv'))

bkg.data <- read.csv(file.path('/uufs/chpc.utah.edu/common/home',
                               'lin-group14/DDR/Roten_InputData/TCCON/Edwards',
                               'df20130720_20201231.public.csv'))

TCCON.data$times <- as.POSIXct(TCCON.data$times,
                               origin = '1970-01-01',
                               tz = 'UTC')
bkg.data$times <- as.POSIXct(bkg.data$times,
                             origin = '1970-01-01',
                             tz = 'UTC')

#aggregate the TCCON data to hourly
TCCON_hourly <- aggregate(xco2 ~ year(times) + month(times) + day(times) + hour(times),
                          data = TCCON.data, mean)
bkg_hourly <- aggregate(xco2 ~ year(times) + month(times) + day(times) + hour(times),
                        data = bkg.data, mean)

names(TCCON_hourly) <- c('year', 'month', 'day', 'hour', 'xco2')
names(bkg_hourly) <- c('year', 'month', 'day', 'hour', 'xco2')

enhancement_hourly <- data.frame(matrix(NA, nrow = 0, ncol = 5))
names(enhancement_hourly) <- names(TCCON_hourly)
for(i in 1:nrow(TCCON_hourly)) {
  
  bkg.val <- subset(bkg_hourly,
                    year == TCCON_hourly$year[i] &
                    month == TCCON_hourly$month[i] &
                    day == TCCON_hourly$day[i] &
                    hour == TCCON_hourly$hour[i])
  
  if(nrow(bkg.val) != 0) {
    add.line <- TCCON_hourly[i,]
    add.line$xco2 <- add.line$xco2 - bkg.val$xco2
  } else if(nrow(bkg.val) == 0) {
    add.line <- TCCON_hourly[i,]
    add.line$xco2 <- NA
  }
    
  enhancement_hourly <- rbind(enhancement_hourly, add.line)
  
  print(i)
}

#format the time
char_add.zero <- function(x, zeros = 1) {
  
  x_out <- NULL
  for(i in 1:length(x)) {
    if(x[i] <= 9) {x_out[i] <- paste0(rep('0', zeros), x[i])}
    if(x[i] > 9) {x_out[i] <- x[i]}
  }
  return(x_out)
}

#plot time series
time.series.plot <- data.frame(subset(TCCON.data, year(times) >= 2015),
                               site = 'Urban')
time.series.plot <- rbind(time.series.plot,
                          data.frame(subset(bkg.data, year(times) >= 2015),
                                     site = 'Background'))
time.series.plot$Display <- 'Total'

#convert enhancement to formatted hourly steps
hourly.times <- paste0(enhancement_hourly$year,
                       char_add.zero(enhancement_hourly$month),
                       char_add.zero(enhancement_hourly$day),
                       char_add.zero(enhancement_hourly$hour))
hourly.times <- as.POSIXct(hourly.times,
                           format = '%Y%m%d%H',
                           origin = '1970-01-01',
                           tz = 'UTC')
enhancement_hourly$timestamp <- hourly.times
enhancement_hourly <- subset(enhancement_hourly,
                             year(timestamp) >= 2015)

#properly format
add.time.series.plot <- data.frame(times = enhancement_hourly$timestamp,
                                   xco2 = enhancement_hourly$xco2,
                                   xco2.uncert = NA,
                                   site = 'Urban',
                                   Display = 'Enhancement')

time.series.plot <- rbind(time.series.plot,
                          add.time.series.plot)

time.series.plot$Display <- factor(time.series.plot$Display,
                                   levels = c('Total', 'Enhancement'))

#get average values
sig.2015 <- subset(enhancement_hourly, year == 2015)
sig.2015$TCCON <- 'Base Year'

sig.OTHER <- subset(enhancement_hourly, year == 2019 | year == 2020)
sig.OTHER$TCCON <- 'Target Years'

mean.2015 <- mean(subset(enhancement_hourly,
                         year == 2015)$xco2, na.rm = TRUE)
err.2015 <- sd(subset(enhancement_hourly,
                      year == 2015)$xco2, na.rm = TRUE)/
  sqrt(nrow(subset(enhancement_hourly, year == 2015)))

mean.other <- mean(subset(enhancement_hourly,
                          year == 2019 | year == 2020)$xco2, na.rm = TRUE)
err.other <- sd(subset(enhancement_hourly,
                       year == 2019 | year == 2020)$xco2, na.rm = TRUE)/
  sqrt(nrow(subset(enhancement_hourly, year == 2019 | year == 2020)))

#TCCON plot
ggplot() +
  ggtitle(expression(paste('XCO'[2], ' Measurements from TCCON Sites'))) +
  labs(subtitle = as.expression(bquote(bar(Delta)[1] == .(round(mean.2015, 3))*'ppm;'~
                                       bar(Delta)[2] == .(round(mean.other, 3))*'ppm;'~
                                       Delta['%'] == .(round(100*(mean.other/mean.2015 - 1), 3))*
                                         '%'))) +
  xlab('Time') +
  ylab(expression(paste('XCO'[2], ' [ppm]'))) +
  geom_rect(data = time.series.plot,
            aes(xmin = as.POSIXct('2015-01-01', tz = 'UTC'),
                xmax = as.POSIXct('2016-01-01', tz = 'UTC'),
                ymin = -Inf, ymax = Inf),
            fill = 'lightgray', color = 'black', alpha = 0.25) +
  geom_rect(data = time.series.plot,
            aes(xmin = as.POSIXct('2020-01-01', tz = 'UTC'),
                xmax = as.POSIXct('2021-01-01', tz = 'UTC'),
                ymin = -Inf, ymax = Inf),
            fill = 'lightgray', color = 'black', alpha = 0.25) +
  geom_line(data = time.series.plot,
            aes(x = times, y = xco2, color = site,
                group = site)) +
  facet_wrap(Display ~ ., ncol = 1, scales = 'free_y') +
  scale_color_discrete(name = 'Location: ') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

#make boxplots
sig.2015 <- subset(enhancement_hourly, year == 2015)
sig.2015$TCCON <- 'Base Year'

sig.OTHER <- subset(enhancement_hourly, year == 2019 | year == 2020)
sig.OTHER$TCCON <- 'Target Years'

mean.2015 <- mean(subset(enhancement_hourly,
            year == 2015)$xco2, na.rm = TRUE)
err.2015 <- sd(subset(enhancement_hourly,
                      year == 2015)$xco2, na.rm = TRUE)/
  sqrt(nrow(subset(enhancement_hourly, year == 2015)))

mean.other <- mean(subset(enhancement_hourly,
                          year == 2019 | year == 2020)$xco2, na.rm = TRUE)
err.other <- sd(subset(enhancement_hourly,
                       year == 2019 | year == 2020)$xco2, na.rm = TRUE)/
  sqrt(nrow(subset(enhancement_hourly, year == 2019 | year == 2020)))

year_info <- data.frame(TCCON = c('Base Year',
                                  'Target Years'),
                        mean.data = c(mean.2015,
                                      mean.other),
                        err.data = c(err.2015,
                                     err.other))

t.test(subset(enhancement_hourly,
              year == 2015)$xco2,
       subset(enhancement_hourly,
              year == 2019 | year == 2020)$xco2,
       alternative = 'less')

print(paste0(round(100*mean.other/mean.2015, 4), '%'))


