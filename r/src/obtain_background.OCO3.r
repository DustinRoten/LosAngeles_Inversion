obtain_background.OCO3 <- function(xco2.path = NULL, sector.list = NULL,
                                   bkg.tolerance = 0, sounding.tolerance = 20,
                                   output.path = NULL) {
  
  #read in the previously generated xco2 file
  if(file.exists(xco2.path)) {
    xco2 <- read.csv(xco2.path)
  } else{return()}
  
  #get the names from the sector list
  sector.names <- sector.list[,1]
  sector.names <- grep(sector.names, pattern = '_',
                       value = TRUE, invert = TRUE)
  sector.names <- gsub(' ', '', sector.names)
  
  #identify and select the columns of interest
  cols.select <- which(names(xco2) %in% sector.names)
  sub.xco2 <- xco2[,cols.select]
  
  #sum the rows
  xco2.totals <- rowSums(sub.xco2)
  
  #identify background values
  xco2.bkg.select <- which(xco2.totals <= bkg.tolerance)
  xco2.bkg <- xco2[xco2.bkg.select,]
  
  #calculate the mean background and error
  oco3.bkg <- data.frame(background = mean(xco2.bkg$xco2),
                         uncert = sd(xco2.bkg$xco2)/nrow(xco2.bkg),
                         bio.background = mean(xco2.bkg$bio),
                         bio.uncertainty = sd(xco2.bkg$bio)/nrow(xco2.bkg))
  
  #write background csv file
  write.csv(oco3.bkg, file = file.path(output.path, 'OCO3_background.csv'),
            row.names = FALSE)
  
}