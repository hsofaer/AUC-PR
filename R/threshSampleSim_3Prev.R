###   aucPR simulations: function to threshold continuous suitabilities, sample observations, and extract covariates

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  21 Nov 2017: based on simulateMany.R (raster version) and threshSampleSim.R (slicker df version)
##  13 Nov 2018: added USGS disclaimer

###  Function inputs:
# TrueContinuousRaster: rasterStack/Brick with 3 layers of true continuous suitability, each representing diff sp
# nsim: number of simulation runs
# nsample: number of points to sample for estimation 
# CovRaster: raster* of covariates, for estimation and prediction (validation in space)
# sample.mask: raster with NAs where I don't want to sample; here using to sample only W of -100

###  Function output:
## dataframe with list columns; nrow = nsim*nspecies (one row for each run of the simulation for each species)
## columns are:
# spID: ID of virtual species
# simnum: integer; simulation number 1:nsim
# sp_sim: concatenation of spID and simnum that gives unique row
# curr.PA: RasterLayer; 0s, 1s of current true PA after probabilistic threshold
# sim.data: df; estimation & prediction data with continuous truth, PA, InEst, covariates

threshSampleSim_3Prev <- function(TrueContRaster, nsim, nsample, CovRaster, sample.mask) {
  
  
  
  # initialize df with list columns
  outdf <- data_frame(spID = paste0("sp", rep(1:nlayers(TrueContRaster), each = nsim)),
                      simnum = rep(1:nsim, times = nlayers(TrueContRaster)),
                      sp_sim = paste0(spID, "_sim", simnum),
                      curr.PA = list(1),
                      sim.data = list(1))
  
  for (i in seq_along(1:length(unique(outdf$sp_sim)))) {
    
    ###
    #####   Threshold continous suitability    ####
    ###
    
    # empty raster:
    curr.PA <- raster(TrueContRaster)
    
    ##  Probabilistic threshold:
    values(curr.PA) <- suppressWarnings(rbinom(ncell(curr.PA), 1, values(TrueContRaster[[grep(outdf$spID[i], names(TrueContRaster))]])))
    # otherwise gives a warning that NA was generated; want NAs there (outside US borders where suitability is NA)
    names(curr.PA) <- "TruePA"
    
    ###
    #####   Sample points for model estimation    ####
    ###
    
    ##  dismo::randomPoints does wo replacement so can't get two in same cell. Also accounts for lat var in area
    # get cell numbers and use to create new 1/0/NA raster
    sample.pts = randomPoints(sample.mask, nsample, cellnumbers = TRUE)
    sample.r <- TrueContRaster[[1]]
    # make zero background
    values(sample.r) <- ifelse(is.na(values(sample.r)), NA, 0)
    sample.r[sample.pts] <- 1
    names(sample.r) <- "InEst"
    
    ###
    ######   Convert to df for estimation and prediction:    #####
    ###
    
    ##  Stack true cont suitability and PA, along with estimation covariates:
    stack.sim = stack(TrueContRaster[[grep(outdf$spID[i], names(TrueContRaster))]], 
                      curr.PA,
                      sample.r,
                      CovRaster)

    sim.data = as.data.frame(stack.sim, xy = TRUE, na.rm = TRUE) 
    sim.data$ID = outdf$sp_sim[i]
    # replace names: TrueContSuit.sp1 etc with just TrueContSuit
    names(sim.data) <- gsub("TrueContSuit.*", "TrueContSuit", names(sim.data))
    
    ###
    #####   Populate that row of the df    #####
    ###
    
    outdf$curr.PA[[i]] <- curr.PA
    outdf$sim.data[[i]] <- sim.data

  }
  
  return(outdf)
  
}