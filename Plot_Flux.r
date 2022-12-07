#plot the sectors
library(raster); library(ggplot2); library(viridis)
library(patchwork)

#initial file paths
home.dir <- '/uufs/chpc.utah.edu/common/home'
work.ext <- 'u1211790/LosAngeles_Inversion'
data.ext <- 'lin-group14/DDR/OCO3_LosAngeles/Bayesian_Inversion'
setwd(file.path(home.dir, work.ext))
source('r/dependencies.r')

#domain
site <- 'Los Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#register the api key for google maps to work
api.key <- read.csv('insert_ggAPI.csv', header = FALSE,
                    col.names = 'api.key')
register_google(key = api.key)

#obtain map
map <- get_map(location = site,
               maptype = 'satellite',
               color = 'bw', zoom = 10)

SAM.dirs <- list.files(file.path(home.dir, data.ext),
                       full.names = TRUE,
                       pattern = 'out_')

#get the list of sectors from the original csv fil
sector.file <- file.path(home.dir, work.ext, 'ext',
                         'defined_vulcan_sectors.csv')
sectors <- gsub(' ', '', read.csv(sector.file)[,1])
sectors <- grep(sectors, pattern = '_LPS', invert = TRUE,
                value = TRUE)

sectors_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
names(sectors_df) <- c('lon', 'lat', 'flux', 'sector')
for(i in 1:length(sectors)) {
  
  #identify the sector
  sector <- sectors[i]
  
  for(j in 1:length(SAM.dirs)) {
    
    #get the SAM path
    SAM <- SAM.dirs[j]
    sector.raster <- list.files(file.path(SAM, 'include', 'sectors'),
                                pattern = sector,
                                full.names = TRUE)
    sector.raster <- brick(sector.raster)
    sector.raster <- calc(sector.raster, mean)
    
    if(j == 1)
      combined.sectors <- sector.raster
    
    if(j > 1)
      combined.sectors <- combined.sectors + sector.raster
    
  }; combined.sectors <- combined.sectors/length(SAM.dirs)
  
  combined.sectors_df <- raster::as.data.frame(combined.sectors,
                                               xy = TRUE)
  combined.sectors_df$sector <- sector
  names(combined.sectors_df) <- names(sectors_df)
  
  #add to the aggregated dataframe
  sectors_df <- rbind(sectors_df, combined.sectors_df)
  
}

sectors_df <- subset(sectors_df, flux != 0)
sectors_df$sector <- separate.camelback(sectors_df$sector)

#create a sector-specific flux map
flux.map.plot <- ggmap(map) +
  ggtitle('Sector-Specific Flux') +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_tile(data = sectors_df,
            aes(x = lon, y = lat, fill = flux)) +
  scale_fill_viridis(name = expression(paste('Mean Flux [', mu, 'mol m'^-2, ' s'^-1, ']')),
                     trans = 'log10') +
  theme_classic() +
  facet_wrap(. ~ sector) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.key.width = unit(0.5, 'in'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1,
                                frame.colour = "black",
                                frame.linewidth = 1))
ggsave(flux.map.plot, filename = 'Out/Flux_Map.jpg',
       device = 'jpg', height = 7, width = 7, units = 'in')

summed_sectors <- aggregate(flux ~ sector, data = sectors_df, sum)
summed_sectors$flux <- 100*summed_sectors$flux/sum(summed_sectors$flux)

#reorder by decreasing values
summed_sectors <-
  summed_sectors[order(summed_sectors$flux, decreasing = TRUE),]

summed_sectors$sector <- factor(summed_sectors$sector,
                              levels = summed_sectors$sector)

contributions.plot <- ggplot() +
  ggtitle('Sector Contributions') +
  xlab('Sector') +
  ylab('Contribution [%]') +
  geom_bar(data = summed_sectors, fill = 'blue', color = 'black',
           aes(x = sector, y = flux), stat = 'identity') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(contributions.plot, filename = 'Out/Flux_Contributions.jpg',
       device = 'jpg', height = 5, width = 5, units = 'in')


