##########################################################
### Run a Bayesian Inversion for the Los Angeles Basin ###
##########################################################

#set up the environment
options(stringsAsFactors = FALSE)

#################################
##### Setup User Parameters #####
#################################
homedir <- '/uufs/chpc.utah.edu/common/home'
input.data <- 'lin-group11/Roten_InputData'
work.ext <- 'u1211790/LosAngeles_Inversion'
workdir <- file.path(homedir, work.ext)
setwd(workdir); source('r/dependencies.r')
library(patchwork); library(geodist); library(RNetCDF)

#domain of interest
site <- 'Los Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599
#################################
#################################



###################################
### Domain Emissions Parameters ###
###################################
#original input resolution
lon_res <- 1/120; lat_res <- 1/120

#specify the emission inventory to use for the domain
domain.inventory <- 'Vulcan'

#direct R to the inventory directory
domain.path <- file.path(homedir, input.extension, 'Vulcan3.0/hourly')
  
background.path <- 

#sector definitions if Vulcan/EDGAR is selected
sector.definitions.path <- 'ext/defined_vulcan_sectors.csv'
  
#files used for downscaling if ODIAC/EDGAR is selected
downscaling <- NA
###################################
###################################


#include background and biospheric components from observations
include_outer <- FALSE
include_bio <- FALSE

#flux units
flux_units <- 'umol/(m2 s)'

#identify the sector of interest here
#use NA for the total XCO2 amount
#use 'Biosphere' for biospheric flux
#must be in both EDGAR and Vulcan sectors!
emissions.category <- NA

#register the api key for google maps to work
api.key <-
  read.csv('insert_ggAPI.csv', header = FALSE, col.names = 'api.key')
register_google(key = api.key)

#identify the location of the output if 'case 1' is used.
if(case == 1) {
  output.directory <- file.path(homedir,
                                paste0('lin-group11/Bayesian_Inversion/',
                                       'F', aggregation.factor))
} else if(case == 2) {
  output.directory <- file.path(homedir,
                                paste0('lin-group11/Bayesian_Inversion/',
                                       'Custom_F', aggregation.factor))
} else if(case == 3) {
  if(ignore.lps) lps <- 'noLPS'
  if(!ignore.lps) lps <- 'LPS'
  output.directory <- file.path(homedir,
                                paste0('lin-group11/Bayesian_Inversion/',
                                       bias, '_', lps))
} else if(case == 4) {
  output.directory <- file.path(homedir,
                                paste0('lin-group11/Bayesian_Inversion/',
                                       'cellwise_corrections'))
}

#simulated XCO2 observation errors
#construct a dataframe with the name of each error source
#second column is he mean value of the error (or NA if it is a %)
#third column is the standard deviation of the eorror (or NA if it is a %)
#fourth column is the percent error (other numeric columns must be NA)
included.errors <- 'ext/included.errors.csv'

#set up file names for required NetCDF files
prior.emissions <- 'prior_emiss.nc'
truth.emissions <- 'truth_emiss.nc'
background.emissions <- 'outer_emiss.nc'
prior.uncertainty <- 'prior_uncert.nc'
bio.emissions <- 'bio_flux.nc'

#set up file names for required RDS files
prior.lonlat <- 'lonlat_domain.rds'
background.lonlat <- 'lonlat_outer.rds'
observation.values <- 'obs.rds'
background.values <- 'background.rds'

#include SLURM options
slurm_options <- list(time = '48:00:00',
                      account = 'lin-np',
                      partition = 'lin-np')
#################################
#################################
#################################

#create a dataframe to be sent to SLURM
p.table <-
  prepare.parallelization(homedir = homedir,
                          input.extension = input.extension,
                          workdir = workdir,
                          api.key = api.key,
                          footprint.dirs = footprint.dirs,
                          odiac.path = odiac.path,
                          odiac.background.path = odiac.background.path,
                          odiac.downscaling = odiac.downscaling,
                          vulcan.path = vulcan.path,
                          vulcan.sector.path = vulcan.sector.path,
                          smurf.path = smurf.path,
                          use.year = use.year,
                          carma.path = carma.path,
                          edgar.path = edgar.path,
                          edgar.downscaling = edgar.downscaling,
                          edgar.sector.path = edgar.sector.path,
                          output.directory = output.directory,
                          site = site,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          emissions.category = emissions.category,
                          included.errors = included.errors,
                          prior.emissions = prior.emissions,
                          truth.emissions = truth.emissions,
                          background.emissions = background.emissions,
                          prior.uncertainty = prior.uncertainty,
                          bio.emissions = bio.emissions,
                          prior.lonlat = prior.lonlat,
                          background.lonlat = background.lonlat,
                          observation.values = observation.values,
                          background.values = background.values,
                          lon_res = lon_res, lat_res = lat_res,
                          aggregation.factor = aggregation.factor,
                          ls = ls, lt = lt,
                          obs.error, include_outer, include_bio,
                          flux_units = flux_units,
                          case = case, bias = bias, ignore.lps = ignore.lps,
                          TEST = TEST, test.which = 1)

####################################################
##### Begin the Inversion Process for each SAM #####
####################################################
#run the Bayesian Inversion in parallel!
if(!TEST) {
  
  #job names
  if(case == 1)
    job <- paste0('Inversion_F', aggregation.factor)
  if(case == 2)
    job <- paste0('Inversion_Custom_F', aggregation.factor)
  if(case == 3) {
    if(ignore.lps) string <- 'noLPS'
    if(!ignore.lps) string <- 'LPS'
    job <- paste0('Inversion_', bias, '_', string)
  }
  if(case == 4)
    job <- paste0('Cellwise_Corrections')
  
  # run jobs in parallel with SLURM
  rslurm::slurm_apply(f = BayesianInversion,
                      params = p.table,
                      jobname = job,
                      nodes = 1,
                      cpus_per_node = nrow(p.table),
                      pkgs = 'base',
                      slurm_options = slurm_options)
  
} else if (TEST) {
  
  p.table <- p.table[1,]
  
  # run one job on an interactive node
  BayesianInversion(site = p.table$site,
                    workdir = p.table$workdir,
                    api.key = p.table$api.key,
                    included.errors = p.table$included.errors,
                    footprint.dirs = p.table$footprint.dirs,
                    odiac.path = p.table$odiac.path,
                    odiac.background.path = p.table$odiac.background.path,
                    odiac.downscaling = p.table$odiac.downscaling,
                    vulcan.path = p.table$vulcan.path,
                    vulcan.sector.path = p.table$vulcan.sector.path,
                    edgar.path = p.table$edgar.path,
                    edgar.downscaling = p.table$edgar.downscaling,
                    carma.path = p.table$carma.path,
                    smurf.path = p.table$smurf.path,
                    use.year = p.table$use.year,
                    output.directory = p.table$output.directory,
                    emissions.category = p.table$emissions.category,
                    edgar.sector.path = p.table$edgar.sector.path,
                    xmin = p.table$xmin,
                    xmax = p.table$xmax,
                    ymin = p.table$ymin,
                    ymax = p.table$ymax,
                    prior.emissions = p.table$prior.emissions,
                    truth.emissions = p.table$truth.emissions,
                    prior.lonlat = p.table$prior.lonlat,
                    background.emissions = p.table$background.emissions,
                    background.lonlat = p.table$background.lonlat,
                    bio.emissions = p.table$bio.emissions,
                    prior.uncertainty = p.table$prior.uncertainty,
                    observation.values = p.table$observation.values,
                    background.values = p.table$background.values,
                    lon_res = p.table$lon_res,
                    lat_res = p.table$lat_res,
                    aggregation.factor = p.table$aggregation.factor,
                    ls = p.table$ls,
                    lt = p.table$lt,
                    obs.error = p.table$obs.error,
                    include_outer = p.table$include_outer,
                    include_bio = p.table$include_bio,
                    flux_units = p.table$flux_units,
                    TEST = TRUE, case = p.table$case,
                    bias = p.table$bias,
                    ignore.lps = p.table$ignore.lps)
  
}