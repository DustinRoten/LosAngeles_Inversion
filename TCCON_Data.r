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

mean.2015 <- mean(subset(enhancement_hourly,
            year == 2015)$xco2, na.rm = TRUE)
err.2015 <- sd(subset(enhancement_hourly,
                      year == 2015)$xco2, na.rm = TRUE)/
  nrow(subset(enhancement_hourly, year == 2015))

mean.other <- mean(subset(enhancement_hourly,
                          year == 2019 | year == 2020)$xco2, na.rm = TRUE)
err.2015 <- sd(subset(enhancement_hourly,
                      year == 2015)$xco2, na.rm = TRUE)/
  nrow(subset(enhancement_hourly, year == 2015))

t.test(subset(enhancement_hourly,
              year == 2015)$xco2,
       subset(enhancement_hourly,
              year == 2019 | year == 2020)$xco2,
       alternative = 'less')
