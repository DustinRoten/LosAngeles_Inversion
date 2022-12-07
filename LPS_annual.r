#############
### GHGRP ###
#############
GHGRP.2015 <- list.files(manufacturing.path,
                         pattern = '2015', full.names = TRUE)
GHGRP.2015 <- read.csv(GHGRP.2015)

GHGRP.2020 <- list.files(manufacturing.path,
                         pattern = '2020', full.names = TRUE)
GHGRP.2020 <- read.csv(GHGRP.2020)

#Industry
GHGRP.2015_industry <- subset(GHGRP.2015,
                              Industry.Type..sectors. != 'Power Plants')
GHGRP.2015_total <- sum(as.numeric(GHGRP.2015_industry$CO2.emissions..non.biogenic.),
                        na.rm = TRUE)

GHGRP.2020_industry <- subset(GHGRP.2020,
                              Industry.Type..sectors. != 'Power Plants')
GHGRP.2020_total <- sum(as.numeric(GHGRP.2020_industry$CO2.emissions..non.biogenic.),
                        na.rm = TRUE)

nrow(GHGRP.2015_industry); nrow(GHGRP.2020_industry)
100*(GHGRP.2020_total/GHGRP.2015_total - 1)
#############
#############


#############
### eGRID ###
#############
#2014
eGRID.2014 <- list.files(powerplant.path, full.names = TRUE,
                         pattern = '2014')
eGRID.2014 <- read.csv(eGRID.2014)
eGRID.2014_total <- sum(as.numeric(eGRID.2014$PLCO2AN),
                        na.rm = TRUE)

#2016
eGRID.2016 <- list.files(powerplant.path, full.names = TRUE,
                         pattern = '2016')
eGRID.2016 <- read.csv(eGRID.2016)
eGRID.2016_total <- sum(as.numeric(eGRID.2016$PLCO2AN),
                        na.rm = TRUE)

eGRID.2015_total <- mean(c(eGRID.2014_total, eGRID.2016_total))

#2020
eGRID.2020 <- list.files(powerplant.path, full.names = TRUE,
                         pattern = '2020')
eGRID.2020 <- read.csv(eGRID.2020)
eGRID.2020_total <- sum(as.numeric(eGRID.2020$PLCO2AN),
                        na.rm = TRUE)

eGRID.change <- 100*(eGRID.2020_total/eGRID.2015_total - 1)
#############
#############
