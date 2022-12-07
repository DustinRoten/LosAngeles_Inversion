#plot the sectors
library(raster); library(ggplot2); library(viridis)
library(patchwork); library(stringr)

#initial file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/OCO3_LosAngeles/Bayesian_Inversion'
setwd(file.path(home.dir, work.ext))
source('r/dependencies.r')

#domain
site <- 'Los Angeles'
local.tz <- 'America/Los_Angeles'
TCCON.site <- c(-118.12702, 34.13634)
TCCON.bkg <- c(-117.8811, 34.9599)
dlon <- dlat <- 0.05

xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#list the available footprint files
list.foots <- list.files(file.path(home.dir, data.ext),
                         full.names = TRUE,
                         pattern = 'out_')

#list the directory names in vector form
foot_df <- str_split_fixed(basename(list.foots),
                           pattern = '_', n = 4)[,2]

#convert to local time
foot_df <- as.POSIXct(foot_df, format = '%Y%m%d%H',
                      tz = 'UTC')
attr(foot_df, 'tzone') <- local.tz

selected.foots <- which(hour(foot_df) >= 11 & hour(foot_df) <= 13)
selected.foots <- list.foots[selected.foots]

for(i in 1:length(selected.foots)) {

  #get the list of footprints
  foot.list <- list.files(file.path(selected.foots[i], 'footprints'),
                          full.names = TRUE)
  foot.list <- list.files(foot.list, full.names = TRUE)
  
  foot.geolocations <- as.data.frame(str_split_fixed(basename(foot.list),
                                                     pattern = '_',
                                                     n = 3)[,2:3])
  foot.geolocations[,2] <- gsub('.nc', '', foot.geolocations[,2])
  foot.geolocations[,1] <- as.numeric(foot.geolocations[,1])
  foot.geolocations[,2] <- as.numeric(foot.geolocations[,2])
  names(foot.geolocations) <- c('lon', 'lat')
  
  #which footprints are over the TCCON site?
  tccon.soundings <- which(foot.geolocations$lon >= TCCON.site[1] - dlon &
                           foot.geolocations$lon <= TCCON.site[1] + dlon &
                           foot.geolocations$lat >= TCCON.site[2] - dlat &
                           foot.geolocations$lat <= TCCON.site[2] + dlat)
  
  if(i == 1)
    all.foots <- foot.list[tccon.soundings]
  if(i > 1)
    all.foots <- c(all.foots, foot.list[tccon.soundings])
  
} #close the subsetting loop

#average the selected footprints
for(i in 1:length(all.foots)) {
  
  tmp <- calc(brick(all.foots[i]), mean)
  
  if(i == 1)
    summed.footprints <- tmp
  if(i > 1)
    summed.footprints <- summed.footprints + tmp
  
}; averaged.footprints <- summed.footprints/length(all.foots)
removeTmpFiles(h=0)

#raster to df
avg.foot_df <- raster::as.data.frame(averaged.footprints,
                                     xy = TRUE)
avg.foot_df <- subset(avg.foot_df, layer != 0)

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

#obtain map
map <- get_map(location = c(mean(c(xmin,xmax)),
                            mean(c(ymin,ymax))),
               maptype = 'terrain',
               color = 'bw', zoom = 9)

influence.plot <- ggmap(map) +
  ggtitle(expression(paste('CO'[2], ' Emission Domain'))) +
  xlab('Longitude') +
  ylab('Latitude') +
  # geom_tile(data = avg.foot_df, alpha = 0.75,
  #           aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(name = expression(paste('Footprint [ppm/', mu, 'mol/',
                                             'm'^2, '/s]')),
                     trans = 'log10') +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = 'black', fill = NA, linetype = 'dashed') +
  # geom_point(shape = 24, fill = 'red', size = 2,
  #            aes(x = TCCON.site[1], y = TCCON.site[2])) +
  # geom_point(shape = 24, fill = 'red', size = 2,
  #            aes(x = TCCON.bkg[1], y = TCCON.bkg[2])) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1,
                                frame.colour = "black",
                                frame.linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.text = element_text(angle = 45, hjust = 1))
ggsave(influence.plot,
       filename = file.path('Out/Results/Influence_Domain.jpg'),
       device = 'jpg', height = 5, width = 5, units = 'in')
