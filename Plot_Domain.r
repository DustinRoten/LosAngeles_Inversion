#plot the domain
library(ggplot2); library(ggmap); library(scales); library(patchwork)
library(lubridate); library(stringr)

home.dir <- '/uufs/chpc.utah.edu/common/home'
work.dir <- 'u1211790/LosAngeles_Inversion'
data.dir <- 'lin-group14/DDR/Roten_InputData'
output.dir <- 'Out/Results'
setwd(file.path(home.dir, work.dir))

#time filter
base.time <- '201501010000-201512311159'
target.time <- '202001010000-202106301159'

#base year
comparison.year <- 2015
tzone <- 'America/Los_Angeles'

#domain
site <- 'Los Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

#obtain map
map <- get_map(location = c(mean(c(xmin,xmax)),
                            mean(c(ymin,ymax))),
               maptype = 'satellite',
               color = 'bw', zoom = 8)

domain.plot <- ggmap(map) +
  ggtitle(site) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = 'yellow', fill = NA, linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(domain.plot, filename = 'Out/Results/Domain.jpg',
       device = 'jpg', height = 3, width = 3, units = 'in')
