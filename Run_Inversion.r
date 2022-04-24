##########################################################
### Run a Bayesian Inversion for the Los Angeles Basin ###
##########################################################
#set up the environment
options(stringsAsFactors = FALSE)

#register the api key for google maps to work
api.key <-
  read.csv('insert_ggAPI.csv', header = FALSE, col.names = 'api.key')
register_google(key = api.key)
##########################################################
##########################################################



#################################
##### Setup User Parameters #####

#################################
#set up working directories
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
#input resolution
lon_res <- 1/120; lat_res <- 1/120

### FFCO2 Emissions ###
#specify the emission inventory to use for the domain
domain.inventory <- 'Vulcan'

#direct R to the inventory directory
domain.path <- file.path(homedir, input.data, 'Vulcan3.0/hourly')
  
background.path <- domain.path #for now, same as domain

#sector definitions if Vulcan/EDGAR is selected
sector.definitions.path <- 'ext/defined_vulcan_sectors.csv'
  
#files used for downscaling if ODIAC/EDGAR is selected
downscaling <- NA

### SMUrF Emissions ###
SMUrF.path <- file.path(homedir, input.data, 'SMUrF')
SMUrF.output <- 'NEE' #GPP, Reco, NEE

#flux units
flux_units <- 'umol/(m2 s)'
###################################
###################################



#######################################
### Set up File and Directory Names ###
#######################################
#identify the output directory
output.directory <- file.path(homedir, 'lin-group11',
                              'OCO3_LosAngeles',
                              'Bayesian_Inversion')

#include error values for the R matrix
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
job <- paste0('LA_Inversion')
#################################
#################################
#################################



##################
### Input Data ###
##################
#original X-STILT footprints (from X-STILT_UrBAnFlux)
footprint.dirs <- list.files(file.path(homedir, 'lin-group11',
                                       'OCO3_LosAngeles',
                                       'raw_footprints'),
                             pattern = 'out_',
                             full.names = TRUE)

#OCO intrument information
instrument <- 'OCO-3';
data.level <- 'L2'
data.version <- 'V10.4r'

#construct the file path from the input parameters above
OCO_3 <- list.files(file.path(homedir, input.data,
                              instrument, data.level),
                    pattern = data.version,
                    full.names = TRUE)

#TCCON data (compare can be NA)
#provide path for coincident TCCON site
TCCON.compare.path <- list.files(file.path(homedir, input.data,
                                           'TCCON', 'Caltech'),
                                 pattern = '*.public.csv',
                                 full.names = TRUE)

#provide path for background TCCON data
TCCON.background.path <- list.files(file.path(homedir, input.data,
                                              'TCCON', 'Edwards'),
                                    pattern = '*.public.csv',
                                    full.names = TRUE)
##################
##################



####################
### Set up SLURM ###
####################

#create dataframe of function arguments
p.table <- data.frame(api.key, homedir, input.data, workdir,
                      site, xmin, xmax, ymin, ymax,
                      lon_res, lat_res, domain.inventory,
                      domain.path, background.path, 
                      sector.definitions.path, downscaling,
                      SMUrF.path, SMUrF.output, flux_units,
                      output.directory, included.errors,
                      prior.emissions, truth.emissions,
                      background.emissions, prior.uncertainty,
                      bio.emissions, prior.lonlat,
                      background.lonlat, observation.values,
                      background.values, footprint.dirs,
                      OCO_3, TCCON.compare.path,
                      TCCON.background.path)

# run jobs in parallel with SLURM
rslurm::slurm_apply(f = BayesianInversion,
                    params = p.table,
                    jobname = job,
                    nodes = 1,
                    cpus_per_node = nrow(p.table),
                    pkgs = 'base',
                    slurm_options = slurm_options)