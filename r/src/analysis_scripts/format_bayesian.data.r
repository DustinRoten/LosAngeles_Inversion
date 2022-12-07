format_bayesian.data <- function(bayesian.path = NULL,
                                 bayesian.name = 'Bayesian_Results.csv',
                                 tzone = NULL) {
  
  #read in the data
  bayesian.data <- read.csv(file.path(bayesian.path, bayesian.name))
  
  #convert to local time
  bayesian.data$time <- as.POSIXct(bayesian.data$time, tz = 'UTC')
  attr(bayesian.data$time, 'tzone') <- tzone
  
  #convert to relative values
  bayesian.data$ScalingFactor <- 100*(bayesian.data$ScalingFactor - 1)
 
  return(bayesian.data)
   
}