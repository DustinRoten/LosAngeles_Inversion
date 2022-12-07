format_aircraft.data <- function(aircraft.path, comparison.year) {
  
  raw.aircraft <- list.files(aircraft.path, full.names = TRUE)
  raw.aircraft <- read.csv(raw.aircraft)
  
  base.year <- subset(raw.aircraft, Year == comparison.year)
  output.year <- subset(raw.aircraft, Year != comparison.year)
  
  rel.year <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(rel.year) <- c('Time', 'Domestic', 'International', 'Cargo')
  for(i in 1:nrow(output.year)) {
    
    #determine the time
    timestamp <- as.POSIXct(paste0(output.year$Year[i], '-',
                                   output.year$Month[i], '-',
                                   '15'),
                            format = '%Y-%m-%d',
                            tz = 'UTC')
    
    tmp <- subset(base.year, Month == output.year$Month[i])
    
    #construc the add line for the rel.year dataframe
    add.line <- data.frame(Time = timestamp,
                           Domestic =
                             100*(output.year$Domestic.Passengers[i]/tmp$Domestic.Passengers - 1),
                           International =
                             100*(output.year$International.Passengers[i]/tmp$International.Passengers - 1),
                           Cargo =
                             100*(output.year$Cargo..Tons.[i]/tmp$Cargo..Tons. - 1))
    
    rel.year <- rbind(rel.year, add.line)
    
  }
  
  return(rel.year)
  
}