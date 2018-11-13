###   aucPR simulations
###   Script to simulate PA from true continuous suitability; fit models

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  21 Nov 2017
##  14 Feb 2018: version based on all bioclim vars
##  1 March 2018: updated based on new sp; removed code for true covs only
##  16 May 2018: reran after masking great lakes
##  13 Nov 2018: added USGS disclaimer

library(raster) # 2.6-7
library(tidyverse) # 1.2.1
theme_set(theme_bw())
library(dismo) # 1.1-4
library(glmnet) # 2.0-13
library(randomForest) # 4.6-12
# library(glmulti) # 1.0.7
library(gbm) # 2.1.3
library(PRROC) # 1.3
library(SDMTools) # 1.1-221
library(gridExtra) # 2.3
# R 3.4.3

setwd("H:/RegularizedClimChange/aucprSim")


###
#####   Import virtual species - 3 prevalences   #####
###

##  True continuous suitability:
ThreeSp = brick('ThreeSp_TrueContSuit_16May2018.grd')
plot(ThreeSp)

# Climate rasterbrick with all 19:
bio19 = read_rds("curr.bio.all19.rds")
compareRaster(ThreeSp, bio19, values = FALSE)

###
#####    Source functions      #####
###

source('threshSampleSim_3Prev.R')

###  Functions for model estimation and prediction, thin wrappers for use with list-column dataframe
# Each one puts model object and estimation time in a list; I then unlist before predicting
##  Each estimation function needs as input:
# df: dataframe with covariates and response; note that only a subset of rows (InEst == 1) should be used for estimation
# respvar: quoted name of response variable
# covvec: character vector of covariates for estimation
##  Each prediction function needs as input: model object, newdf, covvec
source('H:/RegularizedClimChange/SDMsimulations/fit_glmnet.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_glmnet.R')
source('H:/RegularizedClimChange/SDMsimulations/fit_rf.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_rf.R')
source('H:/RegularizedClimChange/SDMsimulations/fit_glm.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_glm.R')
source('H:/RegularizedClimChange/SDMsimulations/fit_gbm.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_gbm.R')
##  The function that fits all the models:
source('fitMultiSim_3Prev.R')

###
#####    Threshold, sample, and extract covariates    ####
###

##  Create a sampling mask to limit sampling to W of -100:
col100 <- colFromX(ThreeSp[[1]], -100)
SampleWest <- ThreeSp[[1]]
SampleWest[, col100:ncol(SampleWest)] <- NA
# replace all non-NA with 1s
values(SampleWest) <- ifelse(is.na(values(SampleWest)), NA, 1)
plot(SampleWest)
nCellWest = cellStats(SampleWest, sum)

##  Call function to threshold and sample
set.seed(76543)
sim3prev = threshSampleSim_3Prev(TrueContRaster = ThreeSp, 
                                 nsim = 3, nsample = 1000, 
                                 CovRaster = bio19,
                                 sample.mask = SampleWest)
head(sim3prev)

plot(sim3prev$curr.PA[[3]])
plot(sim3prev$curr.PA[[4]])
plot(sim3prev$curr.PA[[7]])
head(sim3prev$sim.data[[6]])

# check prevalence:
sim3prev$sim.data[[2]] %>%
      filter(x < -100,
             TruePA == 1) %>%
      nrow() / nCellWest
    
sim3prev$sim.data[[6]] %>%
      filter(x< -100,
             TruePA == 1) %>%
      nrow() / nCellWest
sim3prev$sim.data[[7]] %>%
  filter(x < -100,
         TruePA == 1) %>%
  nrow() / nCellWest

# number of sampled presences:
sum(sim3prev$sim.data[[3]]$TruePA == 1 & sim3prev$sim.data[[3]]$InEst == 1)
sum(sim3prev$sim.data[[5]]$TruePA == 1 & sim3prev$sim.data[[5]]$InEst == 1)
sum(sim3prev$sim.data[[7]]$TruePA == 1 & sim3prev$sim.data[[7]]$InEst == 1)

###
#####   Estimate and predict for each simulation    ####
###

set.seed(011100)
vec.all <- names(bio19)
fit.all <- fitMultiSim_3Prev(SimSp = sim3prev, covvec = vec.all, respvar = "TruePA", glmtype = "none")


####  Check number of GBM trees - want > 1000  - could tweak parameters b/c getting far more
ntrees_GBM <- function(GBMmodel) {
  GBMmodel$gbm.call$best.trees
}

fit.all %>%
  mutate(ntreesGBM = map_dbl(fit.gbm, ntrees_GBM)) %>%
  ungroup() %>%
  summarize(minTree = min(ntreesGBM),
            meanTree = mean(ntreesGBM))

##  Export:
saveRDS(fit.all, "fit.all.rds")

###
#####    Convert to long to calculate performance statistics      #####
###

#  Drop model objects:
pred.long.all <- fit.all %>%
  # keep just estimation and prediction data, and predictions
  dplyr::select(spID, simnum, sp_sim, sim.data, contains("Pred")) %>%
  unnest() %>% 
  dplyr::select(-starts_with("bio")) 
glimpse(pred.long.all) # of course not actually long yet, despite name

all.equal(pred.long.all$sp_sim, pred.long.all$ID) # T

##  Make long:
pred.long.all = pred.long.all %>%
  dplyr::select(-ID) %>%
  gather(algorithm, PredCont, contains('PredCurr')) %>%
  # drop PredCurr. part of algorithm column
  mutate(algorithm = gsub("PredCurr\\.(.*)", "\\1", algorithm))
head(pred.long.all)

pred.long.all %>% as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()


###
#####  Threshold predictions based on optimal maxS+S threshold for each model:   ####
###

optThresh.all = pred.long.all %>%
  filter(InEst == 1) %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  # optimal threshold based on estimation data: can give range of values, so taking mean to get 1 number
  summarize(optmaxSS.Est = mean(SDMTools::optim.thresh(TruePA, PredCont)$`max.sensitivity+specificity`))
head(optThresh.all)

optThresh.all %>%
  ggplot(., aes(algorithm, optmaxSS.Est)) +
  geom_point(aes(color = spID)) +
  ggtitle("Optimal threshold: each point is one simulation run") # pretty close to prevalence except higher for rf

##  Join:
pred.long.all = pred.long.all %>%
  left_join(., optThresh.all) %>%
  # threshold predictions:
  mutate(PredPA = ifelse(PredCont > optmaxSS.Est, 1, 0)) 
head(pred.long.all)

pred.long.all %>%
  filter(simnum == 1,
         InEst == 0) %>%
  sample_n(50000) %>%
  ggplot(., aes(TrueContSuit, PredCont)) +
  geom_point(aes(color = algorithm))


###
#####   Export model predictions    ####
###

saveRDS(pred.long.all, "pred.all_16May2018.rds")

