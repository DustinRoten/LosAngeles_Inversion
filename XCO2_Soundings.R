library(raster); library(ggmap)

#directories
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.dir <- 'u1211790/LosAngeles_Inversion'
data.dir <- 'lin-group14/KaiW/XSTILT_OCO/OCO_obs'
output.dir <- 'Out/Results'
setwd(file.path(home.dir, work.dir))
source('r/dependencies.r')

#city information
city.center <- c(-118.138, 34.004)
local.tz <- 'America/Los_Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599
nc.extent <- c(xmin, xmax, ymin, ymax)

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

#obtain map
map <- get_map(location = 'Los Angeles', maptype = 'satellite',
               color = 'bw', zoom = 10)

#list the overpasses
oco2 <- list.files(file.path(home.dir, data.dir, 'overpass_obs'),
                   pattern = 'LosAngeles', full.names = TRUE)
oco2.idx <- 59

oco2.transect <- read.delim(oco2[oco2.idx], sep = '')[,c(2,1,3)]
names(oco2.transect) <- c('lon', 'lat', 'xco2')
oco2.transect$Source <- 'OCO-2'


oco3 <- list.files(file.path(home.dir, 'lin-group11',
                             'Bayesian_Inversion', 'F1'),
                   pattern = '20200224',
                   full.names = TRUE)
oco3 <- read.csv(file.path(oco3, 'include',
                           'all_xco2_observations.csv'))[,c(2:4)]
names(oco3) <- c('lon', 'lat', 'xco2')
oco3$Source <- 'OCO-3'

ggmap(map) +
  ggtitle(expression(paste('XCO'[2], ' Soundings from an Overpass'))) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_point(data = rbind(oco2.transect, oco3),
             shape = 5, color = 'white',
             aes(x = lon, y = lat)) +
  scale_fill_viridis() +
  facet_wrap(. ~ Source, ncol = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
