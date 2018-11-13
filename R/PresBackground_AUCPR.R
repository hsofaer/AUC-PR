###  AUC-PR evaluation with presence-background data    

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  27 July 2018
##  8 Aug 2018: added survey perf
##  13 Nov 2018: added USGS disclaimer

library(tidyverse)
library(raster)
library(dismo)
library(sf)
library(tabularaster)
library(glmnet) 
library(randomForest) 
library(gbm)
library(SDMTools)
library(PRROC)
library(gridExtra)
library(modEvA)
theme_set(theme_bw())

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Import PA layer for rare sp    ####
###

##  Previous fits, has binary raster:
fit.all = readRDS("fit.all.rds")
head(fit.all)

##  True PA for rare sp: isolate the raster
RarePA <- fit.all %>%
  filter(spID == "sp1") %>%
  dplyr::select(curr.PA)
RarePA <- stack(RarePA$curr.PA)
plot(RarePA)
names(RarePA) <- paste(c("sim1", "sim2", "sim3"), "TruePA", sep = "_")


###
#####   Red brome locations - to sample presence    ####
###

BRRU <- read.csv("BRRU_20180730_Global.csv")
head(BRRU)
table(BRRU$DataSet)

# just the 3 datasets that contribute to background:
BRRU <- BRRU %>%
  filter(DataSet %in% c("bison", "EDDMapS", "gbif")) %>%
  dplyr::select(DataSet, SciName = ITIS_AcceptedName,
                longitude, latitude, ObsYear)
range(BRRU$ObsYear) # 1980-2018

# to sf object:
BRRU <- st_as_sf(BRRU, coords = c("longitude", "latitude"),
                 crs = 4326)
# limit to US:
BRRU <- st_crop(BRRU, RarePA) # 3158 points
# cell numbers of presences:
BRRUcells <- cellnumbers(RarePA[[1]], BRRU)$cell_
# brick with rare sp PA values at those cells:
rBRRU <- brick(RarePA)
values(rBRRU) <- NA
rBRRU[BRRUcells] <- RarePA[BRRUcells]
freq(rBRRU) # differ slightly with TruePA; 129, 138, 128 pres
# make the 0s NAs too:
rBRRU[rBRRU == 0] <- NA
names(rBRRU) <- paste(c("sim1", "sim2", "sim3"), "pres", sep = "_")
plot(rBRRU, colNA = "black")

# extent from which to sample background:
brru.max = as.data.frame(rBRRU, xy = TRUE) %>% # not doing na.rm = TRUE b/c that takes out ones where one layer is NA
  filter(sim1_pres == 1 | sim2_pres == 1 | sim3_pres == 1) %>%
  summarize(brru.xmax = max(x),
            brru.ymax = max(y))

###
#####   Import and prep background locations    #####
###
##  These are from BISON, EDDMapS, and GBIF
##  Records for exotic grass sp, coarsened to 1 per 30m pixel

BG.30m <- read.csv("J:/Projects/NPS/Scripts/USGS_FORT-master/SpeciesOccurrenceData/bckgrndPntsGraminoid.csv")
head(BG.30m)
table(BG.30m$DataSet)

###   Convert background points to raster layer:
# value of 0 where point falls in cell:
rBG <- rasterize(BG.30m[, c("longitude", "latitude")], 
                 RarePA, field = 0)

# drop 1s east of -100, since I don't want those in estimation:
rBG[, colFromX(rBG, -100):ncol(rBG)] <- NA
# also east and north of max brru coordinate:
rBG[, colFromX(rBG, brru.max$brru.xmax):ncol(rBG)] <- NA
rBG[1:rowFromY(rBG, brru.max$brru.ymax), ] <- NA

# mask to get the NAs back at borders (drop MX):
rBG <- mask(rBG, RarePA[[1]])
plot(rBG, colNA = "purple")
freq(rBG)
names(rBG) <- "BG"

# ##  Sample 10K random points from these:
# # just brute force it three times, since randomPoints doesn't seem to be vectorized
# set.seed(4321)
# BG.10k.1 <- randomPoints(rBG, 10000, cellnumbers = TRUE)
# BG.10k.2 <- randomPoints(rBG, 10000, cellnumbers = TRUE)
# BG.10k.3 <- randomPoints(rBG, 10000, cellnumbers = TRUE)
# 
# # brute force assign:
# sampleBG.1 <- raster(rBG)
# values(sampleBG.1) <- NA
# sampleBG.1[BG.10k.1] <- 0
# 
# sampleBG.2 <- raster(rBG)
# values(sampleBG.2) <- NA
# sampleBG.2[BG.10k.2] <- 0
# 
# sampleBG.3 <- raster(rBG)
# values(sampleBG.3) <- NA
# sampleBG.3[BG.10k.3] <- 0
# 
# sampleBG = brick(sampleBG.1, sampleBG.2, sampleBG.3)
# plot(sampleBG, colNA = "purple")
# freq(sampleBG)
# names(sampleBG) <- paste(c("sim1", "sim2", "sim3"), "BG", sep = "_")

###
#####   3 sets of data, with covariates    ####
###

# load covariates:
bio19 <- read_rds("curr.bio.all19.rds")


##  tibble for basis of estimation and prediction:
RarePresBG <- fit.all %>%
  filter(spID == "sp1") %>%
  dplyr::select(spID, simnum, sp_sim, TruePA = curr.PA)

RarePresBG$PBest <- vector('list', 3)
RarePresBG$PBdata <- vector('list', 3)
RarePresBG$MESS <- vector('list', 3)

for (i in 1:3) {
  
  # stack up all the layers
  temp.stack = stack(bio19, 
                     RarePA[[i]],
                     rBRRU[[i]],
                     rBG)
  
  # convert to df, and drop outside of boundaries
  temp.df = as.data.frame(temp.stack, xy = TRUE) %>%
    filter(!is.na(bio1)) 
  
  # clean up the names to drop sim number:
  names(temp.df) <- gsub("sim[[:digit:]]_", "", names(temp.df))
  
  Est.df <- temp.df %>%
    gather(RespType, RespVal, pres, BG) %>%
    filter(!is.na(RespVal)) 
  
  # so I can use the existing functions
  Est.df$InEst <- 1
  
  
  ###  MESS surface:
  PB.bio <- Est.df %>%
    dplyr::select(contains('bio'))
  
  temp.mess <- mess(bio19, PB.bio, full = FALSE)
  
  # add MESS value to the df with predictions:
  temp.df$MESS <- raster::extract(temp.mess, temp.df[, c("x", "y")])
  
  # add all to df:
  RarePresBG$PBest[[i]] <- Est.df
  RarePresBG$PBdata[[i]] <- temp.df
  RarePresBG$MESS[[i]] <- temp.mess
}

head(RarePresBG$PBdata[[1]])
head(RarePresBG$PBest[[3]])
with(RarePresBG$PBest[[3]], table(RespType, RespVal))
with(RarePresBG$PBest[[3]], table(TruePA, RespVal))

RarePresBG$PBdata[[3]] %>%
  filter(!is.na(pres) & !is.na(BG)) %>%
  nrow()
# 127 overlapping points - that's basically all of them!

##  plot MESS surfaces:
plot(RarePresBG$MESS[[1]])
sum(RarePresBG$PBdata[[1]] < 0, na.rm = TRUE)

###
#####   Fit models    ####
###
# use the thin wrappers, but just do it here rather than in function
##  Each prediction function needs as input: model object, newdf, covvec
source('H:/RegularizedClimChange/SDMsimulations/fit_glmnet.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_glmnet.R')
source('H:/RegularizedClimChange/SDMsimulations/fit_rf.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_rf.R')
source('H:/RegularizedClimChange/SDMsimulations/fit_gbm.R')
source('H:/RegularizedClimChange/SDMsimulations/predict_gbm.R')

covvec <- names(bio19)
respvar <- "RespVal"

PBmodel <- RarePresBG %>%
  mutate(
    ##  lasso estimation and prediction: alpha = 1
    est.lasso = map(PBest, fit_glmnet, 
                    respvar = respvar, 
                    covvec = covvec,
                    alpha = 1),
    fit.lasso = map(est.lasso, 'model.glmnet'),
    Pred.lasso = map2(fit.lasso, PBdata, predict_glmnet, covvec = covvec),
    ##  ridge estimation and prediction: alpha = 0
    est.ridge = map(PBest, fit_glmnet,
                    respvar = respvar, 
                    covvec = covvec,
                    alpha = 0),
    fit.ridge = map(est.ridge, 'model.glmnet'),
    Pred.ridge = map2(fit.ridge, PBdata, predict_glmnet, covvec = covvec),
    # random forest estimation and prediction:
    est.rf = map(PBest, fit_rf,
                 respvar = respvar, 
                 covvec = covvec),
    fit.rf = map(est.rf, 'model.rf'),
    Pred.rf = map2(fit.rf, PBdata, predict_rf, covvec = covvec),
    # gbm estimation and prediction:
    est.gbm = map(PBest, fit_gbm,
                  respvar = respvar, 
                  covvec = covvec),
    fit.gbm = map(est.gbm, 'model.gbm'),
    Pred.gbm = map2(fit.gbm, PBdata, predict_gbm, covvec = covvec)
  ) %>%
  dplyr::select(-est.lasso, -est.ridge, -est.rf, -est.gbm) # this dropped PBest too


####  Check number of GBM trees - want > 1000  - could tweak parameters b/c getting far more
ntrees_GBM <- function(GBMmodel) {
  GBMmodel$gbm.call$best.trees
}

PBmodel %>%
  mutate(ntreesGBM = map_dbl(fit.gbm, ntrees_GBM)) %>%
  ungroup() %>%
  summarize(minTree = min(ntreesGBM),
            meanTree = mean(ntreesGBM))

##  Export:
saveRDS(PBmodel, "PBmodel_swBG.rds")
# PBmodel = readRDS("PBmodel_swBG.rds")


###
#####    Plot showing presence background locations      ######
###

states <- map_data("state")

PBab = data_frame(x = rep(-124.5, 2),
                  y = rep(49, 2),
                  label = c('a', 'b'),
                  RespType = factor(c("pres", "BG"),
                                    levels = c("pres", "BG"), 
                                    ordered = TRUE))

p.PBdata = PBmodel$PBest[[1]] %>%
  mutate(RespType = factor(RespType,
                          levels = c("pres", "BG"), 
                          ordered = TRUE)) %>%
  ggplot(.) +
  geom_point(aes(x, y, color = factor(RespVal)),
             size = .75) +
  scale_color_manual(values = c("#440154FF", "darkgoldenrod")) +
  geom_map(data = states, map = states, aes(map_id = region),
           fill = NA, color = "black") +
  facet_wrap(~ RespType) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank()) +
  geom_text(data = PBab, aes(x, y, label = label),
            size = 6) +
  xlim(-124.5, -100) +
  ylim(26, 49)
ggsave("./AUCPR_figures/Supp_PBdata.png", 
       width = 8, height = 5)

###
#####    Convert to long to calculate performance statistics      #####
###

#  Drop model objects:
PredPB <- PBmodel %>%
  # keep just estimation and prediction data, and predictions
  dplyr::select(spID, simnum, sp_sim, PBdata, contains("Pred")) %>%
  unnest() %>% 
  dplyr::select(-starts_with("bio")) 
glimpse(PredPB) 

##  Make long: 
PredPB = PredPB %>%
  mutate(InEst = ifelse(!is.na(pres) | !is.na(BG), 1, 0)) %>%
  gather(algorithm, PredCont, contains('Pred')) %>%
  # drop PredCurr. part of algorithm column
  mutate(algorithm = gsub("Pred\\.(.*)", "\\1", algorithm))
head(PredPB)

PredPB %>% as.data.frame() %>%
  dplyr::select(-pres, -BG) %>%
  filter(!complete.cases(.)) %>%
  nrow()

###
#####  Threshold predictions based on optimal maxS+S threshold for each model:   ####
###

optThresh.PB = PredPB %>%
  filter(InEst == 1) %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  # optimal threshold based on estimation data: can give range of values, so taking mean to get 1 number
  summarize(optmaxSS.Est = mean(SDMTools::optim.thresh(TruePA, PredCont)$`max.sensitivity+specificity`))
optThresh.PB # yikes; 0 for RF models


##  Join:
PredPB = PredPB %>%
  left_join(., optThresh.PB) %>%
  # threshold predictions:
  mutate(PredPA = ifelse(PredCont > optmaxSS.Est, 1, 0)) 
head(PredPB)

PredPB %>%
  group_by(sp_sim, algorithm, TruePA) %>%
  count(PredPA)

# why not all 1s for RF:
PredPB %>%
  filter(algorithm == "rf",
         PredPA == 0) %>%
  head()

###
#####   Export model predictions    ####
###

saveRDS(PredPB, "PredPB_swBG.rds")

PredPB = readRDS("PredPB_swBG.rds")

###
#######   Performance metrics for estimation extent and for West   #####
###
# Drop points used in model estimation

# function to calculate precision based on organization of confusion matrix with 0s before 1s:
precision.cm <- function(cm) {
  cm[2, 2]/(cm[2, 2] + cm[2, 1])
}

PerfWest.PB = PredPB %>%
  # West only:
  filter(InEst == 0,
         x < -100) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
            # create a confusion matrix in a list column
            ConfMat = list(confusion.matrix(TruePA, PredPA)),
            PredPrev = sum(PredPA == 1)/n(),
            TruePrev = sum(TruePA == 1)/n(),
            ncell = n()) %>%
  mutate(precision = map_dbl(ConfMat, precision.cm),
         space = "West",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich")) 
glimpse(PerfWest.PB)

PerfEst.PB = PredPB %>%
  filter(InEst == 0,
         # limit to estimation extent:
         x < brru.max$brru.xmax,
         y < brru.max$brru.ymax) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  mutate(PredCont = ifelse(PredCont < 0, 0, PredCont)) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
            # create a confusion matrix in a list column
            ConfMat = list(confusion.matrix(TruePA, PredPA)),
            PredPrev = sum(PredPA == 1)/n(),
            TruePrev = sum(TruePA == 1)/n(),
            ncell = n(),
            # calibration:
            MillerSlope = MillerCalib(obs = TruePA, pred = PredCont, plot = FALSE)$slope) %>%
  mutate(precision = map_dbl(ConfMat, precision.cm),
         space = "Estimation",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich"),
         MillerAbsDiff1 = abs(MillerSlope - 1)) 

simPerf.PB = bind_rows(PerfWest.PB, PerfEst.PB)

simPerf.PB %>% 
  dplyr::select(-PR.obj, -ConfMat) %>%
  as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()

unique(simPerf.PB$space)

simPerf.PB %>%
  ggplot(., aes(TruePrev, PredPrev, color = algorithm)) +
  geom_point() +
  facet_wrap(~ space) +
  geom_abline(intercept = 0, slope = 1)

###
#####   Plot effect of extent increase   #####
###

##  Change in AUC and aucPR with change in extent:
p.aucPR.extent = simPerf.PB %>%
  ungroup() %>%
  dplyr::select(space, spID, simnum, algorithm, AUC, aucPR) %>%
  gather(metric, value, AUC, aucPR) %>%
  unite(space_metric, space, metric) %>%
  spread(space_metric, value) %>%
  ggplot(., aes(Estimation_AUC, Estimation_aucPR,
                color = algorithm)) +
  geom_segment(aes(xend = West_AUC, yend = West_aucPR),
               arrow = arrow(length = unit(.03, "npc")),
               lwd = .75) +
  xlab("AUC-ROC: change from estimation extent to West") +
  ylab("AUC-PR:\nchange from estimation extent to West") +
  scale_y_continuous(breaks = seq(.25, .7, by = .05),
                     limits = c(.25, .7)) +
  scale_x_continuous(breaks = seq(.75, 1, by = .05),
                     limits = c(.75, 1)) +
  theme_classic() +
  theme(legend.position = c(.25, .8),
        legend.title = element_text(size = 18),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10))) 
ggsave("./AUCPR_figures/Supp_PB_extent.png", 
       width = 7, height = 6)

###
#####   Search performance     #####
###

##  Just the metrics I'll want to join to:
aucEst.PB = PerfEst.PB %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, space,
                AUC, aucPR, precision, MillerAbsDiff1) 
aucEst.PB

####   Sample based on continuous prediction   
###   Estimation extent: sample 250 sites
set.seed(34567)

range(PredPB$PredCont)

SurveyRare.Est.PB = PredPB %>%
  filter(InEst == 0,
         # just rare sp
         spID == "sp1", 
         # limit to estimation extent:
         x < brru.max$brru.xmax,
         y < brru.max$brru.ymax) %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
         space = "Estimation",
         method = "weightedSample") %>% 
  sample_n(size = 250, replace = FALSE, weight = PredCont) %>%
  group_by(spID, simnum, sp_sim, algorithm, space, method) %>%
  summarize(n.found = sum(TruePA),
            n.sample = n(),
            prop.found = n.found/n.sample) %>%
  inner_join(., aucEst.PB)  
SurveyRare.Est.PB


######   Survey based on top_n:
TopN.Est.PB = PredPB %>%
  filter(InEst == 0,
         # just rare sp
         spID == "sp1", 
         # limit to estimation extent:
         x < brru.max$brru.xmax,
         y < brru.max$brru.ymax) %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
         space = "Estimation",
         method = "topN") %>% 
  top_n(250, PredCont) %>%
  group_by(spID, simnum, sp_sim, algorithm, space, method) %>%
  summarize(n.found = sum(TruePA),
            n.sample = n(),
            prop.found = n.found/n.sample) %>%
  inner_join(., aucEst.PB)  
TopN.Est.PB

surveyPB <- bind_rows(SurveyRare.Est.PB, TopN.Est.PB)

###
#####   Plot survey perf vs metrics    #####
###

# calculate correlations to add to plot: not adding to plot anymore, so doesn't need to be this complicated
cor.PB = surveyPB %>%
  rename(AUCPR = aucPR) %>%
  gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
  # order for plot:
  ungroup() %>%
  mutate(metric = factor(metric,
                         levels = c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
                         labels = c("AUC-PR", "AUC-ROC",  "Precision", "Lack of calibration"),
                         ordered = TRUE),
         method = factor(method,
                         levels = c("topN", "weightedSample"),
                         labels = c("Highest-ranked sites", "Weighted sample of sites"),
                         ordered = TRUE)) %>%
  group_by(metric, method) %>%
  summarize(cor = round(cor(value, prop.found), 2)) %>%
  ungroup() %>%
  # brute force add 'value' and 'prop.found' to place them where I want them on the plots
  mutate(value = c(.58, .58, .88, .88, .375, .375, .5, .5),
         prop.found = rep(c(.6, .45), times = 4))


surveyPB %>%
  rename(AUCPR = aucPR) %>%
  gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
  # order for plot:
  ungroup() %>%
  mutate(metric = factor(metric,
                         levels = c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
                         labels = c("AUC-PR", "AUC-ROC",  "Precision", "Lack of calibration"),
                         ordered = TRUE),
         method = factor(method,
                         levels = c("topN", "weightedSample"),
                         labels = c("Highest-ranked sites", "Weighted sample of sites"),
                         ordered = TRUE)) %>%
  ggplot(., aes(value, prop.found)) +
  geom_point(size = 3) +
  facet_grid(method ~ metric, scale = "free", switch = "both") +
  ylab("Proportion of sampled locations with presence") +
  geom_text(data = cor.PB, aes(label = paste0("r = ", cor)),
            size = 6) +
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        strip.background = element_blank())  
ggsave("./AUCPR_figures/Supp_SurveyPerf_PB_8Aug2018.png", width = 10, height = 6)

# ###
# #####    Versions restricted to where there's no extrapolation, defined by MESS    #####
# ###
# 
# PredPB %>%
#   filter(InEst == 0,
#          x < -100,
#          MESS < 0) %>%
#   nrow()
# 
# ###   Performance when excluding negative MESS cells:
# PerfWest.PB.posMESS = PredPB %>%
#   # West only:
#   filter(InEst == 0,
#          x < -100,
#          MESS >= 0) %>%
#   group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
#   summarize(AUC = SDMTools::auc(TruePA, PredCont), 
#             PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
#                                    curve = TRUE, 
#                                    max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
#             # create a confusion matrix in a list column
#             ConfMat = list(confusion.matrix(TruePA, PredPA)),
#             PredPrev = sum(PredPA == 1)/n(),
#             TruePrev = sum(TruePA == 1)/n(),
#             ncell = n()) %>%
#   mutate(precision = map_dbl(ConfMat, precision.cm),
#          space = "West",
#          aucPR = map_dbl(PR.obj, "auc.davis.goadrich")) 
# glimpse(PerfWest.PB.posMESS)
# 
# PerfEst.PB.posMESS = PredPB %>%
#   filter(InEst == 0,
#          # limit to estimation extent:
#          x < brru.max$brru.xmax,
#          y < brru.max$brru.ymax,
#          MESS >= 0) %>%
#   group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
#   summarize(AUC = SDMTools::auc(TruePA, PredCont), 
#             PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
#                                    curve = TRUE, 
#                                    max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
#             # create a confusion matrix in a list column
#             ConfMat = list(confusion.matrix(TruePA, PredPA)),
#             PredPrev = sum(PredPA == 1)/n(),
#             TruePrev = sum(TruePA == 1)/n(),
#             ncell = n()) %>%
#   mutate(precision = map_dbl(ConfMat, precision.cm),
#          space = "Estimation",
#          aucPR = map_dbl(PR.obj, "auc.davis.goadrich")) 
# 
# simPerf.PB.posMESS = bind_rows(PerfWest.PB.posMESS, PerfEst.PB.posMESS)
# 
# simPerf.PB.posMESS %>% 
#   dplyr::select(-PR.obj, -ConfMat) %>%
#   as.data.frame() %>%
#   filter(!complete.cases(.)) %>%
#   nrow()
# 
# ##  Change in AUC and aucPR with change in extent:
# simPerf.PB.posMESS %>%
#   ungroup() %>%
#   dplyr::select(space, spID, simnum, algorithm, AUC, aucPR) %>%
#   gather(metric, value, AUC, aucPR) %>%
#   unite(space_metric, space, metric) %>%
#   spread(space_metric, value) %>%
#   ggplot(., aes(Estimation_AUC, Estimation_aucPR,
#                 color = algorithm)) +
#   geom_segment(aes(xend = West_AUC, yend = West_aucPR),
#                arrow = arrow(length = unit(.03, "npc")),
#                lwd = .75) +
#   xlab("AUC-ROC: change from est to west") +
#   ylab("AUC-PR: change from est to west") +
#   # scale_y_continuous(breaks = seq(.5, .8, by = .05),
#   #                    limits = c(.5, .77)) +
#   # scale_x_continuous(breaks = seq(.65, 1, by = .05), 
#   #                    limits = c(.69, 1)) +
#   theme_classic() +
#   theme(legend.position = c(.7, .8),
#         legend.title = element_text(size = 18),
#         axis.text = element_text(color = "black", size = 14),
#         axis.title = element_text(size = 18),
#         legend.text = element_text(size = 14),
#         axis.title.x = element_text(margin = ggplot2::margin(t = 10))) +
#   annotate(geom = "text", label = 'a', x = .69, y = .77, size = 8)
# ##  same
# 
# 
# # #####   Search performance     
# # 
# # ##  Just the metrics I'll want to join to:
# # aucWest.PB.posMESS = PerfWest.PB.posMESS %>%
# #   dplyr::select(spID, simnum, sp_sim, algorithm, space,
# #                 AUC, aucPR, precision) 
# # 
# # 
# # ####   Sample based on continuous prediction   
# # ###   West: sample 250 sites
# # set.seed(2468)
# # 
# # SurveyRare.west.PB.posMESS = PredPB %>%
# #   filter(InEst == 0,
# #          # just rare sp
# #          spID == "sp1", 
# #          # just in the west
# #          x < -100,
# #          MESS >= 0) %>%
# #   group_by(spID, simnum, sp_sim, algorithm) %>%
# #   mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
# #          space = "West") %>% 
# #   sample_n(size = 250, replace = FALSE, weight = PredCont) %>%
# #   group_by(spID, simnum, sp_sim, algorithm, space) %>%
# #   summarize(n.found = sum(TruePA),
# #             n.sample = n(),
# #             prop.found = n.found/n.sample) %>%
# #   inner_join(., aucWest.PB.posMESS)  
# # glimpse(SurveyRare.west.PB.posMESS)
# # 
# # SurveyRare.west.PB.posMESS %>% as.data.frame() %>%
# #   filter(!complete.cases(.)) %>%
# #   nrow()
# # 
# # ###
# # #####   Simple plot for ms revision: just AUC-PR, AUC-ROC, precision    #####
# # ###
# # 
# # # calculate correlations to add to plot: not adding to plot anymore, so doesn't need to be this complicated
# # SurveyRare.west.PB.posMESS %>%
# #   rename(AUCPR = aucPR) %>%
# #   gather(metric, value, AUC, AUCPR, precision) %>%
# #   group_by(metric) %>%
# #   summarise(cor = round(cor(prop.found, value), 2)) 
# # 
# # 
# # SurveyRare.west.PB.posMESS %>%
# #   rename(AUCPR = aucPR) %>%
# #   gather(metric, value, AUC, AUCPR, precision) %>%
# #   # order for plot:
# #   ungroup() %>%
# #   mutate(metric = factor(metric,
# #                          levels = c("AUCPR", "AUC", "precision"),
# #                          labels = c("AUC-PR", "AUC-ROC",  "Precision"),
# #                          ordered = TRUE)) %>%
# #   ggplot(., aes(value, prop.found,
# #                 color = algorithm)) +
# #   geom_point(size = 3) +
# #   facet_wrap( ~ metric, scale = "free_x", ncol = 3, strip.position = "bottom") +
# #   ylab("Proportion of sampled\nlocations with presence") +
# #   theme(legend.title = element_blank(),
# #         axis.text = element_text(color = "black", size = 10),
# #         axis.title.y = element_text(size = 16),
# #         legend.text = element_text(size = 12),
# #         strip.text = element_text(size = 12),
# #         strip.placement = "outside",
# #         axis.title.x = element_blank(),
# #         strip.background = element_blank()) 
# # # pretty much looks the same - excluding areas of extrapolation doesn't change the findings