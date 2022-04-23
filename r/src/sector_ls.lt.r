sector_ls.lt <- function(truth.path = NULL, prior.path = NULL, sector.path = NULL,
                         sam.timestamp = NULL) {
  
  #list all of the sectors/categories in the directory
  sector.list <- list.files(sector.path, full.names = TRUE)
  
  #grab the two total inventories
  prior <- brick(prior.path)
  truth <- brick(truth.path)
  
  for(i in 1:length(sector.list)) {
    
    #get the sector and its name
    sector <- brick(sector.list[i]); sector[is.na(sector)] <- 0
    sector.name <- separate.camelback(gsub('.nc', '',
                                           basename(sector.list[i])))
    
    #determine the fraction of the total that this sector makes up
    fraction <- (sector/prior); fraction[is.na(fraction)] <- 0
    r2 <- fraction*truth; r2[is.na(r2)] <- 0
   
    ls <-
      variogram_flux.err.ls(r1 = sector, r2 = r2,
                            plot.label = paste0('SAM: ',
                                                strftime(sam.timestamp,
                                                         tz = 'UTC')),
                            output = file.path(obs.dir,
                                               paste0(sector.name, '_',
                                                      'flux_variogram')))
    
  }
  
}
