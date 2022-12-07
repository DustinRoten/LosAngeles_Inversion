format_traffic.data.daily <- function(traffic.data = NULL, tzone = NULL) {
  
  require(stats)
  
  daily.traffic.data <-
    aggregate(perc_change ~ year(time) + month(time) + day(time),
              data = traffic.data, mean)
  
  dates <- paste0(daily.traffic.data$`year(time)`, '-',
                  daily.traffic.data$`month(time)`, '-',
                  daily.traffic.data$`day(time)`, ' ',
                  '12:00')
  
  new.dates <- as.POSIXct(dates,
                          format = '%Y-%m-%d %H:%M',
                          tz = tzone)
  
  new.traffic.data <-
    data.frame(time = new.dates,
               perc_change = daily.traffic.data$perc_change)
  
  return(new.traffic.data)
  
}