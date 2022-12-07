format_traffic.data <- function(traffic.path = NULL, scaling.factors.path = NULL,
                                traffic.name = 'interstate_data.csv',
                                scaling.factors.name = 'traffic_co2.csv',
                                tzone = 'America/Los_Angeles',
                                comparison.year = NULL) {
  
  #read in the traffic dataset (VMT and VHT)
  traffic.data <- read.csv(file.path(traffic.path, traffic.name))
  
  #convert the traffic date/time to POSIX format
  traffic.data$Hour <- tryCatch(expr = as.POSIXct(traffic.data$Hour,
                                                  tz = tzone),
                                error = function() {
                                  stop('Unable to format traffic timestamps.')
                                  return(NA)
                                },
                                finally = as.POSIXct(traffic.data$Hour,
                                                     tz = tzone))
  
  #read in the appropriate scaling factors
  traffic.scaling.factors <- read.csv(file.path(scaling.factors.path,
                                                scaling.factors.name))
  
  #total VMT up across all interstates
  total.traffic <- aggregate(VMT ~ Hour, data = traffic.data, sum)
  
  #use a simple for-loop to convert VMT to CO2
  #scale the traffic data based on CO2/mi
  traffic.years <- unique(year(total.traffic$Hour))
  for(i in 1:length(traffic.years)) {
    
    #match years and scale them to g/mi CO2
    total.traffic$VMT[year(total.traffic$Hour) == traffic.years[i]] <-
      total.traffic$VMT[year(total.traffic$Hour) == traffic.years[i]]*
      traffic.scaling.factors$Year[traffic.scaling.factors$Year ==
                                     traffic.years[i]]
    
  }; total.traffic$VMT <- total.traffic$VMT/1e6
  fix.name <- which(names(total.traffic) == 'VMT')
  names(total.traffic)[fix.name] <- 'mtCO2'
  
  #identify the week number and day for each hour
  total.traffic$week <- week(total.traffic$Hour)
  total.traffic$day <- as.character(wday(total.traffic$Hour, label = TRUE))
  
  #identify the comparison year
  compare.year <- subset(total.traffic, year(Hour) == comparison.year)
  other.years <- subset(total.traffic, year(Hour) != comparison.year)
  
  #create a dataframe of relative values (%)
  rel.values <- data.frame(matrix(NA, nrow = 0, ncol = 2))
  names(rel.values) <- c('time', 'perc_change')
  for(i in 1:nrow(other.years)) {
    
    entry.other.years <- other.years[i,]
    
    #get the corresponding entry from the base year
    entry.compare.year <- subset(compare.year,
                                 week == other.years$week[i] &
                                 day == other.years$day[i] &
                                 hour(Hour) == hour(other.years$Hour[i]))
    
    #calculate the ratio between the year values
    hr.ratio <- entry.other.years$mtCO2/entry.compare.year$mtCO2
    if(length(hr.ratio) == 0) {hr.ratio <- NA}
    
    #add the line to the relative dataframe
    rel.values <- rbind(rel.values,
                        data.frame(time = other.years$Hour[i],
                                   perc_change = 100*(hr.ratio - 1)))
    
  } #closes calculations for relative values (i)
  
  return(rel.values)
}