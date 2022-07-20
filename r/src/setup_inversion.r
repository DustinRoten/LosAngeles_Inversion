#' Run the Bayesian Inversion
setup_inversion <- function(api.key, homedir, input.data, workdir, site,
                            xmin, xmax, ymin, ymax, lon_res, lat_res,
                            domain.inventory, domain.path, background.path, 
                            sector.definitions.path, downscaling, SMUrF.path,
                            SMUrF.output, use.year, flux_units, output.directory,
                            included.errors, prior.emissions, truth.emissions,
                            background.emissions, prior.uncertainty,
                            bio.emissions, prior.lonlat, background.lonlat,
                            observation.values, background.values, footprints.dirs,
                            OCO_3, TCCON.comparison.path, TCCON.lon, TCCON.lat,
                            dlon, dlat, TCCON.background.path, tmp.path) {
  
  #set working directory, load dependenices, load functions
  setwd(workdir); source('r/dependencies.r')
  library(patchwork); library(geodist); #library(RNetCDF)
  library(ggmap)
  
  #set the tmp file storage path
  rasterOptions(tmpdir = tmp.path)
  
  #if needed, add api.key after workdir in the arguments
  register_google(key = api.key)
  
  #read in external data (from /ext)
  sector_df <- read.csv(sector.definitions.path)
  
  #generate prior extent (inner domain)
  nc.extent <- extent(xmin, xmax, ymin, ymax)
  
  #begin the process here...
  for(i in 1:length(footprints.dirs)) {
    
    #' First, the appropriate directories need to be constructed
    #' and the Bayesian Inversion scripts added. This is done using
    #' the `output.directory` variable. In each directory, three more
    #' directories must be included:
    #' *footprints*, *include*, *inversion*
    out.dirs <- file.path(output.directory, basename(footprints.dirs[i]))
    necessary.dirs <- file.path(out.dirs,
                                rep(c('footprints', 'include'),
                                    each = length(out.dirs)))
    
    #create the directories
    for(dir in 1:length(necessary.dirs))
      dir.create(necessary.dirs[dir], recursive = T, showWarnings = F)
    
    #create the operator.dir (footprints)
    operator.dir <- file.path(out.dirs, 'footprints')
    obs.dir <- file.path(out.dirs, 'include')
    
    #' Compile the list of footprints here. The standard output of
    #' X-STILT includes a `footprints` folder of symbolic links to
    #' the actual footprint files. If this directory is present, 
    #' the file path will be automatically updated. If it is NOT 
    #' present, the script will search for footprints within the 
    #' original specified directory.
    listed.dir <- list.files(footprints.dirs[i], full.names = TRUE)
    idx <- which(basename(listed.dir) == 'footprints')
    
    #update directory as needed
    if(length(idx) != 0)
      footprint.list <- list.files(listed.dir[idx],
                                   full.names = TRUE, pattern = '.nc')
    if(length(idx) == 0)
      footprint.list <- list.files(footprints.dirs[i],
                                   full.names = TRUE, pattern = '.nc')
    
    ######################################################
    #' The maximum extent required for the background
    #' emission inventory file must be determined. The
    #' low-tech approach here simply reads in each 
    #' footprint file and records the maximum lons/lats.
    #' Additionally, averaged timesteps will be generated.
    ######################################################
    min.lon <- NULL; max.lon <- NULL
    min.lat <- NULL; max.lat <- NULL
    start.times <- NULL; timesteps <- NULL
    for(j in 1:length(footprint.list)) {
      
      #read in the raster, remove all zero weights
      footprint <- brick(footprint.list[j])
      
      #convert the footprint to a dataframe
      footprint_df <- raster::as.data.frame(footprint, xy = TRUE)
      
      #some footprints may only have one timestep
      if(ncol(footprint_df) == 3) {
        summed.flux.values <- footprint_df[,3:ncol(footprint_df)]
      } else if(ncol(footprint_df) > 3) {
        summed.flux.values <- rowSums(footprint_df[,3:ncol(footprint_df)])
      }
      tmp_df <- data.frame(footprint_df[,1:2], summed.flux.values)
      tmp_df <- subset(tmp_df, summed.flux.values > 0)
      
      #record the number of layers
      timesteps[j] <- nlayers(footprint)
      
      #Determine the min/max values here
      min.lon[j] <- min(tmp_df$x); max.lon[j] <- max(tmp_df$x)
      min.lat[j] <- min(tmp_df$y); max.lat[j] <- max(tmp_df$y)
      
      #' Record the start time of each footprint.
      #' THE START TIME MUST BE IN THE NAME OF THE FILE!
      start.times[j] <-
        unlist(strsplit(basename(footprint.list[j]), split = '_'))[1]
      
      cat(paste0('Constraining footprint domain: ',
                 round(100*j/length(footprint.list), 2), '%     '),
          '\r')
      
    }; remove('footprint'); remove('tmp_df')
    
    #Overwrite the unnecessary vectors and replace with a single
    #min/max value.
    #' Testing Note: For SAM_2020022419, the values are as follows:
    #' min.lon <- -121.6875; max.lon <- -115.4708
    #' min.lat <- 32.79583; max.lat <- 41.10417
    min.lon <- min(min.lon); max.lon <- max(max.lon)
    min.lat <- min(min.lat); max.lat <- max(max.lat)
    outer.extent <- extent(min.lon, max.lon, min.lat, max.lat)
    max.layers <- max(timesteps)
    
    #construct the averaged timesteps here
    #HOURLY steps are assumed!
    #strip off seconds!!
    pos.start.times <- as.POSIXlt(start.times,
                                  format = '%Y%m%d%H%M',
                                  tz = 'UTC')
    avg.start.time <- mean(pos.start.times)
    avg.start.time <- avg.start.time - second(avg.start.time)
    
    #format the average start time into a YYYYMMDDHHMM timestamp
    timestr <- strftime(avg.start.time,
                        format = '%Y%m%d%H%M',
                        tz = 'UTC')
    
    #build the sequence here
    step.in.seconds <- 3600
    avg.timesteps <-
      seq(as.numeric(avg.start.time),
          as.numeric(avg.start.time) - step.in.seconds*max.layers,
          -1*step.in.seconds)
    
    #remove the start time
    #(its not included in the footprint time steps)
    #reorder the list to be chronologicial
    timestep.list <- as.POSIXct(avg.timesteps,
                                origin = '1970-01-01',
                                tz = 'UTC')[2:length(avg.timesteps)]
    timestep.list <- sort(timestep.list)
    ######################################################
    
    
    
    #########################################
    #' Construct the modified footprint files
    #########################################
    #' When a footprint is generated, the timesteps will be
    #' mapped to the averaged timesteps from above. The footprint
    #' will also be cropped according to the outer domain that
    #' was constrained above.
    make_foots_ncdf4(footprint.list = footprint.list,
                     timestr = timestr,
                     averaged.timesteps = timestep.list,
                     operator.dir = operator.dir,
                     nc.extent = outer.extent)
    gc() #garbage collection
    #########################################
    
    
    
    ############################
    #' Make the `prior_emiss.nc`
    ############################
    make_vulcan.prior_ncdf4(site, vulcan.path = domain.path,
                            sector.list = read.csv(sector.definitions.path),
                            times = timestep.list,
                            inner.extent = nc.extent,
                            obs.dir, inner.name = prior.emissions,
                            lps.names = c('_LPS'))
    gc() #garbage collection
    ############################


    
    ####################################################
    #' Make the `bio_flux.nc` file using `SMUrF` output.
    ####################################################
    #make the layered (NetCDF) file for the biospheric flux (from SMUrF)
    outer.domain.path <- file.path(obs.dir, '../footprints/footprint_1')
    outer.domain <- raster(list.files(outer.domain.path,
                                      full.names = TRUE))
    make_bio_ncdf4(SMUrF.path, SMUrF.output,
                   times = timestep.list,
                   match.year = use.year,
                   match.raster = outer.domain,
                   bio.output = file.path(obs.dir,
                                          bio.emissions))
    gc() #garbage collection
    ####################################################
    
    
    
    #######################################
    ### Determine Background from TCCON ###
    #######################################
    #save the background value to a csv file in the 'include' directory
    obtain_background.TCCON(TCCON.background.path, timestr,
                            output.path = obs.dir)
    #######################################
    
    
    
    ###############################
    ### Incorporate XCO2 Values ###
    ###############################
    obtain_XCO2(obs.dir,
                TCCON.bkg.path = file.path(obs.dir,
                                           'TCCON_background.csv'),
                OCO3.path = OCO_3,
                SMUrF.path = file.path(obs.dir,
                                       'bio_flux.nc'),
                sector.list = read.csv(sector.definitions.path))
    gc() #garbage collection
    ###############################
    ###############################
    
    
    
    ################################################
    ### Determine the TCCON Value for Comparison ###
    ################################################
    xco2.path <- file.path(obs.dir, 'xco2.csv')
    obtain_comparison.TCCON(TCCON.comparison.path = TCCON.comparison.path,
                            xco2.path = xco2.path, timestr = timestr,
                            output.path = obs.dir, TCCON.lon = TCCON.lon,
                            TCCON.lat = TCCON.lat, dlon = dlon,
                            dlat = dlat)
    ################################################
    ################################################
    
    
    
    #############################################
    ### Attempt to Calculate OCO-3 Background ###
    #############################################
    obtain_background.OCO3(xco2.path,
                           sector.list = read.csv(sector.definitions.path),
                           bkg.tolerance = 0.01, sounding.tolerance = 20,
                           output.path = obs.dir)
    #############################################
    #############################################
    
  } #closes footprint dirs 
} #close function
  