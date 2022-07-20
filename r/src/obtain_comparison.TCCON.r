obtain_comparison.TCCON <- function(TCCON.comparison.path = NULL, xco2.path = NULL,
                                    timestr = NULL, output.path = NULL, tz = 'UTC',
                                    hr.min = 0.5, hr.max = 1.0, TCCON.lon = NULL,
                                    TCCON.lat = NULL, dlon = NULL, dlat = NULL) {
  
  #change the timestr to a POSIX object
  timestr.posix <- as.POSIXct(timestr, format = '%Y%m%d%H%M',
                              origin = '1970-01-01', tz = tz)
  
  #read in the data from the TCCON background site
  TCCON <- read.csv(TCCON.comparison.path)
  TCCON$times <- as.POSIXct(TCCON$times, origin = '1970-01-01',
                            tz = tz)
  
  #subset the TCCON data to the appropriate time window
  TCCON.vals <- subset(TCCON,
                       times <= timestr.posix + hr.min*(60*60) &
                         times >= timestr.posix - hr.min*(60*60))
  
  #if there isn't enough data in the first window, extend it
  if(nrow(TCCON.vals) < 20)
    TCCON.vals <- subset(TCCON,
                         times <= timestr.posix + hr.max*(60*60) &
                           times >= timestr.posix - hr.max*(60*60))
  
  #determine the signal calculations
  TCCON_signal <-
    data.frame(signal = mean(TCCON.vals$xco2),
               uncert = sd(TCCON.vals$xco2)/sqrt(nrow(TCCON.vals)))
  
  #obtain the OCO-3 sounding values
  xco2 <- read.csv(xco2.path)
  xco2 <- subset(xco2,
                 (lon >= TCCON.lon - dlon & lon <= TCCON.lon + dlon) &
                 (lat >= TCCON.lat - dlat & lat <= TCCON.lat + dlat))
  oco3_signal <-
    data.frame(oco3.signal = mean(xco2$xco2),
               oco3.uncert = sd(xco2$xco2)/sqrt(nrow(xco2)))
  
  
  #add the oco3 data to the TCCON data
  TCCON_signal <- data.frame(TCCON_signal, oco3_signal)
  
  #save the information in a background file
  write.csv(TCCON_signal,
            file = file.path(output.path, 'TCCON_signal.csv'),
            row.names = FALSE)
  
}