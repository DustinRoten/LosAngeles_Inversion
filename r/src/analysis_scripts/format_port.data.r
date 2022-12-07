format_port.data <- function(port.data = NULL, comparison.year = NULL) {
  
  comparison.data <- subset(port.data, Year == comparison.year)
  monthly.data <- subset(port.data, Year != comparison.year)
  
  relative.data <- data.frame(matrix(NA,
                                     nrow = nrow(monthly.data),
                                     ncol = 4))
  names(relative.data) <- c('Port', 'Year', 'Month', 'Rel.TEU')
  for(i in 1:nrow(monthly.data)) {
    
    tmp <- subset(comparison.data,
                  Port == monthly.data$Port[i] &
                  Month == monthly.data$Month[i])
    
    relative.data$Port[i] <- monthly.data$Port[i]
    relative.data$Year[i] <- monthly.data$Year[i]
    relative.data$Month[i] <- monthly.data$Month[i]
    relative.data$Rel.TEU[i] <- 100*(monthly.data$TEU[i]/tmp$TEU - 1)
    
  }
  
  #convert to numeric dates
  relative.data$Month[relative.data$Month == 'January'] <- 1
  relative.data$Month[relative.data$Month == 'February'] <- 2
  relative.data$Month[relative.data$Month == 'March'] <- 3
  relative.data$Month[relative.data$Month == 'April'] <- 4
  relative.data$Month[relative.data$Month == 'May'] <- 5
  relative.data$Month[relative.data$Month == 'June'] <- 6
  relative.data$Month[relative.data$Month == 'July'] <- 7
  relative.data$Month[relative.data$Month == 'August'] <- 8
  relative.data$Month[relative.data$Month == 'September'] <- 9
  relative.data$Month[relative.data$Month == 'October'] <- 10
  relative.data$Month[relative.data$Month == 'November'] <- 11
  relative.data$Month[relative.data$Month == 'December'] <- 12
  
  #convert to posix time
  relative.data$timestamp <-
    as.POSIXct(paste0(relative.data$Year, '-',
                      relative.data$Month, '-15'),
               format = '%Y-%m-%d', tz = 'UTC')
  
  return(relative.data)
  
}