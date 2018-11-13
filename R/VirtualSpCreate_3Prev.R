###   Creation of virtualspecies for aucPR simulations
##  Following vignette: http://borisleroy.com/files/virtualspeciestutorial.html
##  Used virtualspecies package for simulation of response curves and true suitability raster (not thresholding/sampling)

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  17 Nov 2017 
##  13 Feb 2018: updated to export all bioclim in addition to 3 I'm using to simulate sp; no change to sp
##  1 March 2018: updated to fix prevalence
##  16 May 2018: masked Great Lakes in continuous output
##  13 Nov 2018: add USGS disclaimer

library(raster) # 2.6-7
library(virtualspecies) # 1.4-2
library(tidyverse) # 1.2.1
library(maptools) # 0.9-2
# R 3.4.3

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Import climate covariates    ####
###

####  Worldclim
## Worldclim version 1.4; downloaded 11/17/2017; 5 minute resolution
## Note that getData isn't pulling from worldclim 2.0 yet, even though it's available
## Units: temp in C*10; precip in mm
curr.bio = getData('worldclim', var = 'bio', res = 5)
proj4string(curr.bio) # wgs84
names(curr.bio)


##  Subset to USA:
USA = getData("GADM", country = "USA", level = 1) 
# drop AK, HI
USA = USA[!USA$NAME_1 %in% c("Alaska", "Hawaii"), ]

curr.bio = crop(curr.bio, USA)
curr.bio = mask(curr.bio, USA)
plot(curr.bio, 1)

##  I can't handle including the great lakes in the US outline...
lakes = maps::map('lakes', regions = "Great Lakes")
IDs <- sapply(strsplit(lakes$names, ":"), function(x) x[2])
lakes = map2SpatialPolygons(lakes, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
curr.bio = mask(curr.bio, lakes, inverse = TRUE)


##  Fix units:
# temp vars: divide by 10; no changes needed for precip (mm)
bioclim.divide.vec = c(10, 10, 1, 10, 10, 10, 10, 10, 10, 10, 10, rep(1, 8))
curr.bio = curr.bio / bioclim.divide.vec

plot(curr.bio)

##  Export all 19:
write_rds(curr.bio, "curr.bio.all19.rds")

##  Subset to variables I'll use:
curr.bio = raster::subset(curr.bio, c("bio1", "bio5", "bio12"))

##  Export:
write_rds(curr.bio, "curr.bio.rds")

###
#####   Generate suitability relationships    ####
###
##  Want 3 species that vary in prevalence; increase SD to get higher prev
##  Targeting prevalence of .5, .25, .05

##  Will want to calculate prevalence based on west of the 100th parallel:
WestExtent = extent(curr.bio)
WestExtent@xmax <- -100
West = crop(curr.bio[[1]], WestExtent)
nCellWest = length(na.omit(values(West)))
nCellUS = length(na.omit(values(curr.bio[[1]])))

###
#####   Sp1: target prevalence = 0.05     #####
###

###   Generate sp based on environmental relationships:
##  Using multiplicative form:
par.sp1 = formatFunctions(bio1 = c(fun = 'dnorm', mean = 20, sd = 8),
                          bio5 = c(fun = 'dnorm', mean = 45, sd = 8),
                          bio12 = c(fun = 'dnorm', mean = 150, sd = 77))

sp1 = generateSpFromFun(raster.stack = curr.bio,
                        parameters = par.sp1,
                        species.type = "multiplicative",
                        plot = TRUE)
plotResponse(sp1)
names(sp1$suitab.raster) <- "TrueContSuit.sp1"


# check if these two ways are the same
set.seed(123)
sp1.probprev <- calc(sp1$suitab.raster, fun = function(x) rbinom(1, 1, x))
plot(sp1.probprev)
prob.v2 <- raster(sp1$suitab.raster)
set.seed(123)
values(prob.v2) <- suppressWarnings(rbinom(ncell(prob.v2), 1, values(sp1$suitab.raster)))
all.equal(sp1.probprev, prob.v2) # T; I think I prefer the second way

# check prevalence:
as.data.frame(sp1.probprev, xy = TRUE, na.rm = TRUE) %>%
  filter(x < -100) %>%
  summarize(nCellWest = n(),
            TruePrev = sum(layer)/n())

hist(sp1$suitab.raster)

###
#####   Species 2: target prevalence = 0.25    ####
###

##  Using multiplicative form:
par.sp2 = formatFunctions(bio1 = c(fun = 'dnorm', mean = 20, sd = 20),
                          bio5 = c(fun = 'dnorm', mean = 45, sd = 20),
                          bio12 = c(fun = 'dnorm', mean = 150, sd = 173))

sp2 = generateSpFromFun(raster.stack = curr.bio,
                        parameters = par.sp2,
                        species.type = "multiplicative",
                        plot = TRUE)
plotResponse(sp2)
names(sp2$suitab.raster) <- "TrueContSuit.sp2"

##  prevalence based on prob thresh:
prob.sp2 <- raster(sp2$suitab.raster)
values(prob.sp2) <- suppressWarnings(rbinom(ncell(prob.sp2), 1, values(sp2$suitab.raster)))
plot(prob.sp2)

# check prevalence:
as.data.frame(prob.sp2, xy = TRUE, na.rm = TRUE) %>%
  filter(x < -100) %>%
  summarize(nCellWest = n(),
            TruePrev = sum(layer)/n())


###
#####   Species 3: target prevalence = 0.5    ####
###

##  Using multiplicative form:
par.sp3 = formatFunctions(bio1 = c(fun = 'dnorm', mean = 20, sd = 46),
                          bio5 = c(fun = 'dnorm', mean = 45, sd = 46),
                          bio12 = c(fun = 'dnorm', mean = 150, sd = 460))

sp3 = generateSpFromFun(raster.stack = curr.bio,
                        parameters = par.sp3,
                        species.type = "multiplicative",
                        plot = TRUE)
plotResponse(sp3)
names(sp3$suitab.raster) <- "TrueContSuit.sp3"

##  prevalence based on prob thresh:
prob.sp3 <- raster(sp3$suitab.raster)
values(prob.sp3) <- suppressWarnings(rbinom(ncell(prob.sp3), 1, values(sp3$suitab.raster)))
plot(prob.sp3)

# check prevalence:
as.data.frame(prob.sp3, xy = TRUE, na.rm = TRUE) %>%
  filter(x < -100) %>%
  summarize(nCellWest = n(),
            TruePrev = sum(layer)/n())


###
######    Export virtual species     #####
###


####   Save virtual species in this format:   
write_rds(sp1, 'sp1_16May2018.rds')
write_rds(sp2, 'sp2_16May2018.rds')
write_rds(sp3, 'sp3_16May2018.rds')


##  Save the three true continuous rasters as a stack:
ThreeSp = stack(sp1$suitab.raster,
                sp2$suitab.raster,
                sp3$suitab.raster)
plot(ThreeSp)
write_rds(ThreeSp, 'ThreeSp_TrueContSuit_16May2018.rds')
writeRaster(ThreeSp, 'ThreeSp_TrueContSuit_16May2018.grd')
