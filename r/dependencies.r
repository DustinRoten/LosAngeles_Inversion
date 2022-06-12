### Source R scripts for Dien's X-STILT, 05/23/2018, DW
rsc <- dir('r/src', pattern = '.*\\.r$',
           full.names = T, recursive = F)
invisible(lapply(rsc, source))

#additional scripts for plotting
rsc.plots <- dir('r/src/analysis_scripts', pattern = '.*\\.r$',
           full.names = T, recursive = F)
invisible(lapply(rsc.plots, source))

rsc.nested <- dir('r/src/analysis_scripts/nested_functions',
                  pattern = '.r$', full.names = T,
                  recursive = F)
invisible(lapply(rsc.nested, source))

#scripts for the cellwise corrections
cellwise <- dir('r/src/inversion_scripts',
                pattern = '.r$', full.names = T)
invisible(lapply(cellwise, source))

#scripts for the oco functions
oco_fun <- dir('r/src/oco_functions',
               pattern = '.r$', full.names = T)
invisible(lapply(oco_fun, source))

# Other relevant packages
library(ggplot2); library(viridis); library(patchwork)
library(geosphere); library(geodist); library(raster)
library(ggmap); library(ncdf4); library(dplyr)
library(lutz); library(lubridate); library(stringr)
library(reshape2); library(patchwork)
