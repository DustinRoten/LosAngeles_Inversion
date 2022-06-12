make_vulcan.prior_ncdf4 <- function(site, vulcan.path, sector.list, times,
                                    inner.extent, obs.dir, inner.name,
                                    lps.names) {
  
  #remove other indicators beginning with '_' in the sector list
  new.sector.list <- grep(x = sector.list$category,
                          pattern = '_', invert = TRUE)
  new.vulcan.sector_df <- sector.list[new.sector.list,]
  
  #get the list of sectors making up the large point sources
  .LPS <- grep(x = sector.list$category, pattern = '_LPS')
  .LPS <- sector.list[.LPS,]
  
  #check the format of the input times
  if(!any(class(times) == 'POSIXt'))
    stop('Input times must be POSIXt format!')
  
  #obtain an arbitrary footprint to set up the reference raster
  ref.path <- file.path(obs.dir, '../footprints')
  ref.path <- list.files(ref.path, pattern = 'footprint_1$',
                         full.names = TRUE)
  ref.raster <- brick(list.files(ref.path, full.names = TRUE))[[1]]
  
  make_vulcan.sectors(vulcan_emiss = vulcan.path,
                      sector_df = new.vulcan.sector_df,
                      date.times = times,
                      ref.raster = ref.raster,
                      inner.extent = inner.extent,
                      output.path = obs.dir)
  ##############################
  ##############################
  ##############################
  
} #closes function