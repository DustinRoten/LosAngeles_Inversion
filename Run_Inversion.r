##########################################################
### Run a Bayesian Inversion for the Los Angeles Basin ###
##########################################################
#set up the environment
options(stringsAsFactors = FALSE)

#delegate space for potential tmp files from the process
tmp.path <- '/scratch/general/vast/u1211790'
if(!dir.exists(tmp.path)) dir.create(tmp.path)
##########################################################
##########################################################



#################################
##### Setup User Parameters #####
#################################
#set up working directories
homedir <- '/uufs/chpc.utah.edu/common/home'
input.data <- 'lin-group14/DDR/Roten_InputData'
work.ext <- 'u1211790/LosAngeles_Inversion'
workdir <- file.path(homedir, work.ext)
setwd(workdir); source('r/dependencies.r')

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')

#domain of interest
site <- 'Los Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#Comparison TCCON
TCCON.lon <- -118.127; TCCON.lat <- 34.136
dlon <- 0.05; dlat <- 0.05
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
SMUrF.output <- 'NEE_mean' #GPP, Reco, NEE
use.year <- 2019

#flux units
flux_units <- 'umol/(m2 s)'
###################################
###################################



#######################################
### Set up File and Directory Names ###
#######################################
#identify the output directory
output.directory <- file.path(homedir, 'lin-group14',
                              'DDR', 'OCO3_LosAngeles',
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
slurm_options <- list(time = '72:00:00',
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
footprints.path <- file.path(homedir, input.data,
                             'XSTILT_output',
                             gsub(' ', '', site))
footprints.dirs <- list.files(footprints.path,
                              full.names = TRUE,
                              pattern = 'out_')

#OCO instrument information
instrument <- 'OCO-3';
data.version <- 'V10.4r'

#construct the file path from the input parameters above
OCO_3 <- list.files(file.path(homedir, input.data, instrument),
                    pattern = data.version,
                    full.names = TRUE)

#TCCON data (compare can be NA)
#provide path for coincident TCCON site
TCCON.comparison.path <- list.files(file.path(homedir, input.data,
                                              'TCCON', 'Caltech'),
                                    pattern = '*.public.qc.csv',
                                    full.names = TRUE)

#provide path for background TCCON data
TCCON.background.path <- list.files(file.path(homedir, input.data,
                                              'TCCON', 'Edwards'),
                                    pattern = '*.public.qc.csv',
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
                      SMUrF.path, SMUrF.output, use.year,
                      flux_units, output.directory, included.errors,
                      prior.emissions, truth.emissions,
                      background.emissions, prior.uncertainty,
                      bio.emissions, prior.lonlat,
                      background.lonlat, observation.values,
                      background.values, footprints.dirs,
                      OCO_3, TCCON.comparison.path,
                      TCCON.lon, TCCON.lat, dlon, dlat,
                      TCCON.background.path, tmp.path)
p.table <- p.table[37,]

#run jobs in parallel with SLURM
# rslurm::slurm_apply(f = setup_inversion,
#                     params = p.table,
#                     jobname = job,
#                     nodes = 3,
#                     cpus_per_node = 25,
#                     pkgs = 'base',
#                     slurm_options = slurm_options)

setup_inversion(api.key = p.table$api.key,
                homedir = p.table$homedir,
                input.data = p.table$input.data,
                workdir = p.table$workdir,
                site = p.table$site,
                xmin = p.table$xmin, xmax = p.table$xmax,
                ymin = p.table$ymin, ymax = p.table$ymax,
                lon_res = p.table$lon_res,
                lat_res = p.table$lat_res,
                domain.inventory = p.table$domain.inventory,
                domain.path = p.table$domain.path,
                background.path = p.table$background.path, 
                sector.definitions.path = p.table$sector.definitions.path,
                downscaling = p.table$downscaling,
                SMUrF.path = p.table$SMUrF.path,
                SMUrF.output = p.table$SMUrF.output,
                use.year = p.table$use.year,
                flux_units = p.table$flux_units,
                output.directory = p.table$output.directory,
                included.errors = p.table$included.errors,
                prior.emissions = p.table$prior.emissions,
                truth.emissions = p.table$truth.emissions,
                background.emissions = p.table$background.emissions,
                prior.uncertainty = p.table$prior.uncertainty,
                bio.emissions = p.table$bio.emissions,
                prior.lonlat = p.table$prior.lonlat,
                background.lonlat = p.table$background.lonlat,
                observation.values = p.table$observation.values,
                background.values = p.table$background.values,
                footprints.dirs = p.table$footprints.dirs,
                OCO_3 = p.table$OCO_3,
                TCCON.comparison.path = p.table$TCCON.comparison.path,
                TCCON.lon = p.table$TCCON.lon,
                TCCON.lat = p.table$TCCON.lat,
                dlon = p.table$dlon, dlat = p.table$dlat,
                TCCON.background.path = p.table$TCCON.background.path,
                tmp.path = p.table$tmp.path)
