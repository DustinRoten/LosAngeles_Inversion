obtain_background.TCCON <- function(TCCON.background.path = NULL, timestr = NULL,
                                    output.path = NULL, tz = 'UTC', hr.min = 0.5,
                                    hr.max = 1.0) {
  
  #change the timestr to a POSIX object
  timestr.posix <- as.POSIXct(timestr, format = '%Y%m%d%H%M',
                              origin = '1970-01-01', tz = tz)
  
  #read in the data from the TCCON background site
  TCCON <- read.csv(TCCON.background.path)
  TCCON$times <- as.POSIXct(TCCON$times, origin = '1970-01-01',
                            tz = tz)
  
  #subset the TCCON data to the appropriate time window
  TCCON.bkg.vals <- subset(TCCON,
                           times <= timestr.posix + hr.min*(60*60) &
                           times >= timestr.posix - hr.min*(60*60))
  
  #if there isn't enough data in the first window, extend it
  if(nrow(TCCON.bkg.vals) < 20)
    TCCON.bkg.vals <- subset(TCCON,
                             times <= timestr.posix + hr.max*(60*60) &
                               times >= timestr.posix - hr.max*(60*60))
  
  #determine the background calculations
  TCCON_background <-
    data.frame(background = mean(TCCON.bkg.vals$xco2),
               uncert = sd(TCCON.bkg.vals$xco2)/sqrt(nrow(TCCON.bkg.vals)))
  
  #save the information in a background file
  write.csv(TCCON_background,
            file = file.path(output.path, 'TCCON_background.csv'),
            row.names = FALSE)
  
}