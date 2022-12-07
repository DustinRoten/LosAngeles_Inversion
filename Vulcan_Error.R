vulcan_emiss <- list.files(file.path(home.dir, data.dir, 'Vulcan3.0/data'),
                           full.names = TRUE)

categories <- c('airport', 'cement', 'cmv', 'commercial', 'elec_prod',
                'industrial', 'nonroad', 'onroad', 'rail', 'residential')

year <- '2015'

nc.extent <- 

get_vulcan_emiss_error <- function(vulcan_emiss = NULL, categories = NA,
                                   nc.extent = NA, year = NA) {
  
  #list all the mean values
  mean.emissions_path <- grep(vulcan_emiss, pattern = '_Emissions',
                              value = TRUE)
  mean.rasters <- grep(list.files(mean.emissions_path, full.names = TRUE),
                       pattern = 'total|AK',
                       value = TRUE, invert = TRUE)
  
  #list all the high values
  hi.emissions_path <- grep(vulcan_emiss, pattern = '_Hi',
                              value = TRUE)
  hi.rasters <- grep(list.files(hi.emissions_path, full.names = TRUE),
                     pattern = 'total|AK',
                     value = TRUE, invert = TRUE)
  
  #list all the low values
  lo.emissions_path <- grep(vulcan_emiss, pattern = '_Lo',
                            value = TRUE)
  lo.rasters <- grep(list.files(lo.emissions_path, full.names = TRUE),
                     pattern = 'total|AK',
                     value = TRUE, invert = TRUE)
  
  for(cat in categories) {
    
    #mean raster
    mn.raster <- grep(mean.rasters, pattern = cat,
                      value = TRUE)
    mn.raster <- brick(mn.raster)
    layer.idx <- grep(names(mn.raster), pattern = as.character(year))
    mn.raster <- mn.raster[[layer.idx]]
    mn.total <- cellStats(mn.raster, sum)
    remove('mn.raster')

    #high raster
    hi.raster <- grep(hi.rasters, pattern = cat,
                      value = TRUE)
    hi.raster <- brick(hi.raster)
    layer.idx <- grep(names(hi.raster), pattern = as.character(year))
    hi.raster <- hi.raster[[layer.idx]]
    hi.total <- cellStats(hi.raster, sum)
    remove('hi.raster')
    
    #low raster
    lo.raster <- grep(lo.rasters, pattern = cat,
                      value = TRUE)
    lo.raster <- brick(lo.raster)
    layer.idx <- grep(names(lo.raster), pattern = as.character(year))
    lo.raster <- lo.raster[[layer.idx]]
    lo.total <- cellStats(lo.raster, sum)
    remove('lo.raster')
    
    hi.diff <- 100*(hi.total - mn.total)/mn.total
    lo.diff <- 100*(lo.total - mn.total)/mn.total
    
    print(paste0(cat, ' - ', round(max(c(abs(hi.diff), abs(lo.diff))),3), '%'))
    
  }
  
}

get_vulcan_emiss_error(vulcan_emiss, categories, nc.extent, year)
