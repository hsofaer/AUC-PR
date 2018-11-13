###   aucPR simulations: performance metric computation and comparison

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  17 Feb 2017: updating figures
##  1 March 2018: updated to new sp prev
##  3 May 2018: confusion matrix figure across threholds
##  17 May 2018: updated to version that masks Great Lakes; incorporated JAH figure tweaks
##  13 Nov 2018: added USGS disclaimer

library(tidyverse) # 1.2.1
theme_set(theme_bw())
library(PRROC) # 1.3
library(SDMTools) # 1.1-221
library(gridExtra) # 2.3

setwd("H:/RegularizedClimChange/aucprSim")


###
#####   Import model predictions     #####
###

pred.all = read_rds("pred.all_16May2018.rds")
head(pred.all)

pred.all %>%
  as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()

###
#######   Performance metrics for West and for L48   #####
###
# Drop points used in model estimation

# function to calculate precision based on organization of confusion matrix with 0s before 1s:
precision.cm <- function(cm) {
  cm[2, 2]/(cm[2, 2] + cm[2, 1])
}


####   For models based on all covariates:
PerfWest.all = pred.all %>%
  # West only:
  filter(InEst == 0,
         x < -100) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
            ROC.obj = list(roc.curve(PredCont, weights.class0 = TruePA, 
                                     curve = TRUE, 
                                     max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)),
            rankcor = cor(PredCont, TrueContSuit, method = "spearman"),
            # create a confusion matrix in a list column
            ConfMat = list(confusion.matrix(TruePA, PredPA)),
            PredPrev = sum(PredPA == 1)/n(),
            TruePrev = sum(TruePA == 1)/n(),
            ncell = n()) %>%
  mutate(kappa = map_dbl(ConfMat, Kappa),
         sensitivity = map_dbl(ConfMat, sensitivity),
         specificity = map_dbl(ConfMat, specificity),
         TSS = sensitivity + specificity - 1,
         PCC = 100 * map_dbl(ConfMat, prop.correct),
         precision = map_dbl(ConfMat, precision.cm),
         space = "West",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich"), # gives Davis and Goadrich 2006 method; see PRROC manual
         minAUCPR = 1 + ((1-TruePrev)*log(1-TruePrev))/TruePrev,
         aucNPR = (aucPR - minAUCPR)/(1 - minAUCPR)) 
head(PerfWest.all)

# check precision
PerfWest.all$ConfMat[[1]]
PerfWest.all$precision[1]

PerfL48.all = pred.all %>%
  filter(InEst == 0) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)), 
            ROC.obj = list(roc.curve(PredCont, weights.class0 = TruePA, 
                                     curve = TRUE, 
                                     max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)),
            rankcor = cor(PredCont, TrueContSuit, method = "spearman"),
            # create a confusion matrix in a list column
            ConfMat = list(confusion.matrix(TruePA, PredPA)),
            PredPrev = sum(PredPA == 1)/n(),
            TruePrev = sum(TruePA == 1)/n(),
            ncell = n()) %>%
  mutate(kappa = map_dbl(ConfMat, Kappa),
         sensitivity = map_dbl(ConfMat, sensitivity),
         specificity = map_dbl(ConfMat, specificity),
         TSS = sensitivity + specificity - 1,
         PCC = 100 * map_dbl(ConfMat, prop.correct),
         precision = map_dbl(ConfMat, precision.cm),
         space = "L48",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich"),
         minAUCPR = 1 + ((1-TruePrev)*log(1-TruePrev))/TruePrev,
         aucNPR = (aucPR - minAUCPR)/(1 - minAUCPR)) 

simPerf.all = bind_rows(PerfWest.all, PerfL48.all)
head(simPerf.all)

simPerf.all %>% 
  dplyr::select(-PR.obj, -ROC.obj, -ConfMat) %>%
  as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()

##  Check prevalence:
simPerf.all %>%
  ggplot(., aes(TruePrev, PredPrev)) +
  geom_point(aes(color = space,
                 shape = spID)) +
  geom_abline(slope = 1, intercept = 0)

simPerf.all %>%
  ungroup() %>%
  dplyr::select(sp_sim, space, TruePrev) %>%
  distinct()

###   Export simulation performance metrics:
saveRDS(simPerf.all, "simPerf.all.rds")


###
#####   Plot effect of extent increase   #####
###

simPerf.all = readRDS("simPerf.all.rds")


##  Change in AUC and aucPR with change in extent:
p.aucPR.extent = simPerf.all %>%
  ungroup() %>%
  # just the algorithms I'm using:
  filter(algorithm %in% c("lasso", "ridge", "rf", "gbm")) %>%
  dplyr::select(space, spID, simnum, algorithm, AUC, aucPR) %>%
  gather(metric, value, AUC, aucPR) %>%
  unite(space_metric, space, metric) %>%
  spread(space_metric, value) %>%
  mutate(spID = factor(spID, 
                       levels = c("sp1", "sp2", "sp3"), 
                       labels = c("Low", "Medium", "High"))) %>%
  ggplot(., aes(West_AUC, West_aucPR)) +
  geom_segment(aes(xend = L48_AUC, yend = L48_aucPR,
                   color = factor(spID)),
               arrow = arrow(length = unit(.03, "npc")),
               lwd = .75) +
  scale_color_manual(values = c("yellow3", "darkgoldenrod3", "brown"),
                     name = "Prevalence") +
  xlab("AUC-ROC: change from West to L48") +
  ylab("AUC-PR: change from West to L48") +
  scale_y_continuous(breaks = seq(.5, .8, by = .05),
                     limits = c(.5, .77)) +
  scale_x_continuous(breaks = seq(.65, 1, by = .05), 
                     limits = c(.69, 1)) +
  theme_classic() +
  theme(legend.position = c(.7, .8),
        legend.title = element_text(size = 18),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10))) +
  annotate(geom = "text", label = 'a', x = .69, y = .77, size = 8)
# printed below with next figure

##  which are the steepest models:
simPerf.all %>%
  ungroup() %>%
  # just the algorithms I'm using:
  filter(algorithm %in% c("lasso", "ridge", "rf", "gbm")) %>%
  dplyr::select(space, spID, simnum, algorithm, AUC, aucPR) %>%
  gather(metric, value, AUC, aucPR) %>%
  unite(space_metric, space, metric) %>%
  spread(space_metric, value) %>%
  mutate(spID = factor(spID, 
                       levels = c("sp1", "sp2", "sp3"), 
                       labels = c("Species 1", "Species 2", "Species 3"))) %>%
  ggplot(., aes(West_AUC, West_aucPR)) +
  geom_segment(aes(xend = L48_AUC, yend = L48_aucPR,
                   linetype = factor(spID),
                   color = algorithm),
               arrow = arrow(length = unit(.03, "npc")),
               lwd = .5) +
 # scale_color_manual(values = c("yellow3", "darkgoldenrod3", "brown")) +
  xlab("AUC: change from West to L48") +
  ylab("AUCPR: change from West to L48") +
  ylim(.5, .77) +
  scale_x_continuous(breaks = seq(.75, 1, by = .05), 
                     limits = c(.67, 1)) +
  theme_classic()

simPerf.all %>%
  filter(space == "West") %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, 
                TruePrev, optmaxSS.Est) %>%
  distinct() %>%
  ggplot(., aes(TruePrev, optmaxSS.Est)) +
  geom_point(aes(color = algorithm, 
                 shape = spID))

###
#####    Addition of true negatives     #####
###

# number of western cells:
pred.all %>%
  ungroup() %>%
  dplyr::select(x, y) %>%
  distinct() %>%
  filter(x < -100) %>%
  nrow() # 57465; will add this many TN predictions

# How many cells in CONUS:
pred.all %>%
  ungroup() %>%
  dplyr::select(x, y) %>%
  distinct() %>%
  nrow()


# sp_sim and algorithm combinations for setup of fake data:
SpSimAlg = pred.all %>%
  ungroup() %>%
  dplyr::select(sp_sim, algorithm) %>%
  distinct()

###  Create a fake dataset full of a single very low value
FakeTN <- data_frame(sp_sim = rep(SpSimAlg$sp_sim, times = 57465),
                     algorithm = rep(SpSimAlg$algorithm, times = 57465)) %>%
  mutate(x = -99.999,
         y = -99.999, 
         TrueContSuit = 0.00001,
         TruePA = 0,
         InEst = 0,
         PredCont = 0.00001,
         PredPA = 0,
         type = "fakeTN")
head(FakeTN)
table(FakeTN$sp_sim, FakeTN$algorithm)

##  Restrict to West and combine real and fakeTN obs
SimExtend.all = pred.all %>%
  ungroup() %>%
  filter(x < -100,
         InEst == 0) %>%
  dplyr::select(sp_sim, algorithm,
                x, y, TrueContSuit, TruePA, InEst, PredCont, PredPA) %>%
  mutate(type = "simpred") %>%
  bind_rows(., FakeTN)
head(SimExtend.all)
table(SimExtend.all$sp_sim, SimExtend.all$type)

SimExtend.all %>%
  as.data.frame %>%
  filter(!complete.cases(.)) %>%
  nrow() # 0

##  Performance with fake TN added:
perf.withFakeTN.all = SimExtend.all %>%
  group_by(sp_sim, algorithm) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = FALSE, min.compute = FALSE, rand.compute = FALSE)),
            TruePrev = sum(TruePA == 1)/n()) %>%
  mutate(aucPR = map_dbl(PR.obj, "auc.davis.goadrich"),
         perfType = "withFakeTN") %>%
  dplyr::select(-PR.obj)
head(perf.withFakeTN.all)

perf.FakeExtend = PerfWest.all %>%
  ungroup() %>%
  dplyr::select(sp_sim, algorithm, 
                AUC, TruePrev, aucPR) %>%
  mutate(perfType = "simPerf.West") %>%
  bind_rows(., perf.withFakeTN.all) 
head(perf.FakeExtend)

# quick check vs prevalence:
perf.FakeExtend %>%
  gather(aucType, AUCvalue, AUC, aucPR) %>%
  ggplot(., aes(TruePrev, AUCvalue)) +
  geom_point(aes(color = aucType))


##  plot change in AUC and aucPR:
p.TN.aucPR = perf.FakeExtend %>%
  # just the algorithms I'm using:
  filter(algorithm %in% c("lasso", "ridge", "rf", "gbm")) %>%
  dplyr::select(-TruePrev) %>%
  gather(metric, value, AUC, aucPR) %>%
  unite(space_metric, perfType, metric) %>%
  spread(space_metric, value) %>%
  # recreate spID column:
  separate(sp_sim, c("spID", "simnum"), sep = "_", remove = FALSE) %>%
  mutate(spID = factor(spID, 
                       levels = c("sp1", "sp2", "sp3"), 
                       labels = c("Low", "Medium", "High"))) %>%
  ggplot(., aes(simPerf.West_AUC, simPerf.West_aucPR)) +
  geom_segment(aes(xend = withFakeTN_AUC, yend = withFakeTN_aucPR,
                   color = factor(spID)),
               arrow = arrow(length = unit(.03, "npc")),
               lwd = .75) +
  scale_color_manual(values = c("yellow3", "darkgoldenrod3", "brown"),
                     name = "Prevalence") +
  xlab("AUC-ROC: adding true negatives") +
  ylab("AUC-PR: adding true negatives") +
  scale_y_continuous(breaks = seq(.5, .8, by = .05),
                     limits = c(.5, .77)) +
  scale_x_continuous(breaks = seq(.65, 1, by = .05), 
                     limits = c(.69, 1)) +
  theme_classic() +
  theme(legend.position = c(.2, .2),
        legend.title = element_text(size = 18),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10))) +
  annotate(geom = "text", label = 'b', x = .69, y = .77, size = 8)

pdf("./AUCPR_figures/ExtentTN_AUCvsAUCPR_17May2018.pdf", width = 14, height = 6)
grid.arrange(p.aucPR.extent, p.TN.aucPR, nrow = 1)
dev.off()

# for word version:
png("ExtentTN_AUCvsAUCPR_16May2018.png", 
    width = 14, height = 6, units = "in", res = 500)
grid.arrange(p.aucPR.extent, p.TN.aucPR, nrow = 1)
dev.off()

# ###
# #####   Plot of confusion matrix with different threshold    #####
# ###
# 
# ##  Extract cells of confusion matrix (based on max(Sens+Spec) threshold)
# simPerf.all$ConfMat[[1]]  # note 0s come before ones, so organization is off from standard
# simPerf.all$ConfMat[[1]][2, 2] # for TP
# 
# compareConfMat = simPerf.all %>%
#   filter(algorithm != "glmnet") %>%
#   dplyr::select(spID, simnum, sp_sim, algorithm, space,
#                 ConfMat, sensitivity, specificity,
#                 aucPR, AUC) %>%
#   # NOTE THE ODD ORGANIZATION OF THE CONFUSION MATRIX (b/c 0 before 1):
#   mutate(TP = map_dbl(ConfMat, ~.[2, 2]),
#          FP = map_dbl(ConfMat, ~.[2, 1]),
#          FN = map_dbl(ConfMat, ~.[1, 2]),
#          TN = map_dbl(ConfMat, ~.[1, 1])) %>%
#   dplyr::select(-ConfMat)
# head(as.data.frame(compareConfMat))
# 
# all.equal(compareConfMat$sensitivity, compareConfMat$TP/(compareConfMat$TP + compareConfMat$FN)) # T
# all.equal(compareConfMat$specificity, compareConfMat$TN/(compareConfMat$FP + compareConfMat$TN)) # T
# 
# ##  Calculate change in confusion matrix from west to L48
# SpaceChange = compareConfMat %>%
#   dplyr::select(-sensitivity, -specificity, -aucPR, -AUC) %>%
#   gather(cell, value, TP, FP, FN, TN) %>%
#   spread(space, value) %>%
#   # calculate change from west to L48:
#   mutate(Change = L48 - West,
#          # order as in confusion matrix:
#          cell = factor(cell, levels = c("TP", "FP", "FN", "TN"), ordered = TRUE))
# head(SpaceChange)
# 
# SpaceChange %>%
#   filter(Change != 0) %>%
#   ggplot(., aes(Change)) +
#   geom_histogram(aes(fill = algorithm), bins = 50) +
#   facet_grid(spID ~ cell, scales = "free_y")
# 
# # decide on limits
# SpaceChange %>%
#   filter(simnum == 1,
#          algorithm != "glmnet") %>%
#   group_by(cell) %>%
#   summarize(maxN = max(Change))
# 
# SpaceChange %>%
#   # create max coordinates - using sqrt to scale for area
#   mutate(max.x = ifelse(cell %in% c("FP", "TN"), sqrt(Change), -sqrt(Change)),
#          max.y = ifelse(cell %in% c("TP", "FP"), sqrt(Change), -sqrt(Change))) %>%
#   filter(simnum == 1) %>%
#   ggplot(., aes(ymin = 0, xmin = 0)) +
#   geom_rect(aes(xmax = max.x, ymax = max.y),
#             fill = NA, color = "black") +
#   facet_grid(spID ~ algorithm) +
#   theme_void()
#   
# ##  Add an additional threshold (e.g. .25 * opt) to get more in other boxes?
# pred.all = read_rds("pred.all_1March2018.rds")
# head(pred.all)
# 
# ##  scaling factor by which to multiply threshold
# thresh.scale = .25
# 
# pred.lowthresh = pred.all %>%
#   filter(algorithm != "glmnet",
#          InEst == 0) %>%
#   mutate(PredPA.newthresh = ifelse(PredCont > optmaxSS.Est*thresh.scale, 1, 0))
# head(pred.lowthresh)
# 
# ##  redo confusion matrices based on new threshold:
# newThresh.L48 = pred.lowthresh %>%
#   group_by(spID, simnum, sp_sim, algorithm) %>%
#   summarize(ConfMat = list(confusion.matrix(TruePA, PredPA.newthresh))) %>%
#   mutate(TP = map_dbl(ConfMat, ~.[2, 2]),
#          FP = map_dbl(ConfMat, ~.[2, 1]),
#          FN = map_dbl(ConfMat, ~.[1, 2]),
#          TN = map_dbl(ConfMat, ~.[1, 1]),
#          space = "L48") %>%
#   dplyr::select(-ConfMat)
# head(newThresh.L48)
# 
# newThresh.West = pred.lowthresh %>%
#   filter(x < -100) %>%
#   group_by(spID, simnum, sp_sim, algorithm) %>%
#   summarize(ConfMat = list(confusion.matrix(TruePA, PredPA.newthresh))) %>%
#   mutate(TP = map_dbl(ConfMat, ~.[2, 2]),
#          FP = map_dbl(ConfMat, ~.[2, 1]),
#          FN = map_dbl(ConfMat, ~.[1, 2]),
#          TN = map_dbl(ConfMat, ~.[1, 1]),
#          space = "West") %>%
#   dplyr::select(-ConfMat)
# 
# newThresh = bind_rows(newThresh.L48, newThresh.West) %>%
#   gather(cell, value, TP, FP, FN, TN) %>%
#   spread(space, value) %>%
#   # calculate change from west to L48:
#   mutate(Change = L48 - West,
#          # order as in confusion matrix:
#          cell = factor(cell, levels = c("TP", "FP", "FN", "TN"), ordered = TRUE),
#          Type = "LowThresh")
# head(newThresh)
# 
# head(SpaceChange)
# 
# # combine and plot:
# TwoThresh = SpaceChange %>%
#   ungroup() %>%
#   mutate(Type = "OptThresh") %>%
#   bind_rows(., newThresh) %>%
#   # create max coordinates - using sqrt to scale for area
#   mutate(max.x = ifelse(cell %in% c("FP", "TN"), sqrt(Change), -sqrt(Change)),
#          max.y = ifelse(cell %in% c("TP", "FP"), sqrt(Change), -sqrt(Change)),
#          # fix species labels:
#          spID = factor(spID, levels = c("sp1", "sp2", "sp3"),
#                        labels = c("Low", "Medium", "High")),
#          Type = factor(Type, levels = c("LowThresh", "OptThresh"),
#                        labels = c("Low", "Optimal"))) 
# head(TwoThresh)
# 
# p.twoThresh = TwoThresh %>%
#   filter(simnum == 1) %>%
#   ggplot(., aes(ymin = 0, xmin = 0)) +
#   geom_rect(aes(xmax = max.x, ymax = max.y,
#                 linetype = Type),
#             fill = NA, size = .7, color = "black") +
#     scale_linetype_manual(values = c("twodash", "solid"),
#                           name = "Threshold") +
#   facet_grid(spID ~ algorithm, switch = "both") +
#   theme_void() +
#   xlab("Algorithm") +
#   ylab("Simulated species prevalence") +
#   theme(strip.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.margin = unit(c(0, 0, .25, .25), "cm"))
# ggsave("ChangeConfMat_1May2018.png", width = 6, height = 5)


###
####    Confusion matrix visualization across gradient of thresholds    #####
###

##  Create sequence of thresholds and compare to continuous prediction
ThreshValue = c(.025, .05, .075, seq(.1, .9, by = .1))
names(ThreshValue) = paste0("PA_", ThreshValue)

predPA.new = as.tibble(t(sapply(pred.all$PredCont, function(x) 1*(x > ThreshValue)))) # multiply by 1 to change T/F to 1/0
head(predPA.new)

predPA.multiThresh = bind_cols(pred.all, predPA.new) %>%
  filter(algorithm != "glmnet",
         InEst == 0) %>%
  # gather threshold columns to long:
  gather(PA_thresh, value, contains("PA_")) 
head(as.data.frame(predPA.multiThresh))

# quick spot check:
predPA.multiThresh %>%
  ungroup() %>%
  select(sp_sim, PredCont, PA_thresh, value) %>%
  sample_n(30) %>%
  as.data.frame()

##  Confusion matrices at the two extents:
multiThresh.L48 = predPA.multiThresh %>%
  group_by(spID, simnum, sp_sim, algorithm, PA_thresh) %>%
  summarize(ConfMat = list(confusion.matrix(TruePA, value))) %>%
  mutate(TP = map_dbl(ConfMat, ~.[2, 2]),
         FP = map_dbl(ConfMat, ~.[2, 1]),
         FN = map_dbl(ConfMat, ~.[1, 2]),
         TN = map_dbl(ConfMat, ~.[1, 1]),
         space = "L48") %>%
  dplyr::select(-ConfMat)
head(multiThresh.L48)

multiThresh.West = predPA.multiThresh %>%
  filter(x < -100) %>%
  group_by(spID, simnum, sp_sim, algorithm, PA_thresh) %>%
  summarize(ConfMat = list(confusion.matrix(TruePA, value))) %>%
  mutate(TP = map_dbl(ConfMat, ~.[2, 2]),
         FP = map_dbl(ConfMat, ~.[2, 1]),
         FN = map_dbl(ConfMat, ~.[1, 2]),
         TN = map_dbl(ConfMat, ~.[1, 1]),
         space = "West") %>%
  dplyr::select(-ConfMat)

##  Combine and prep for plot:
multiThresh = bind_rows(multiThresh.L48, multiThresh.West) %>%
  ungroup() %>%
  gather(cell, value, TP, FP, FN, TN) %>%
  spread(space, value) %>%
  # calculate change from west to L48:
  mutate(Change = L48 - West,
         # order as in confusion matrix:
         cell = factor(cell, levels = c("TP", "FP", "FN", "TN"), ordered = TRUE),
         # create max coordinates - using sqrt to scale for area
         max.x = ifelse(cell %in% c("FP", "TN"), sqrt(Change), -sqrt(Change)),
         max.y = ifelse(cell %in% c("TP", "FP"), sqrt(Change), -sqrt(Change)),
         # fix species labels:
         spID = factor(spID, levels = c("sp1", "sp2", "sp3"),
                       labels = c("Low", "Medium", "High")),
         PA_thresh = factor(PA_thresh, levels = paste0("PA_", ThreshValue),
                            labels = ThreshValue)) 
head(multiThresh)

##  Note opt thresh in title manually: didn't do this for figures I'm not using
pred.all %>%
  filter(simnum == 1,
         algorithm != "glmnet") %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  select(optmaxSS.Est) %>%
  distinct() %>%
  arrange(algorithm) %>%
  as.data.frame()

p.gbm = multiThresh %>%
  filter(simnum == 1,
         algorithm == "gbm") %>%
  ggplot(., aes(ymin = 0, xmin = 0)) +
  geom_rect(aes(xmax = max.x, ymax = max.y),
            fill = NA, size = .7, color = "black") +
  facet_grid(spID ~ PA_thresh, switch = "both") +
  theme_void() +
  xlab("Threshold") +
  ylab("Simulated species prevalence") +
  ggtitle("GBM:") +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0, 0, .25, .25), "cm"),
        plot.title = element_text(size = 24))

p.lasso = multiThresh %>%
  filter(simnum == 1,
         algorithm == "lasso") %>%
  ggplot(., aes(ymin = 0, xmin = 0)) +
  geom_rect(aes(xmax = max.x, ymax = max.y),
            fill = NA, size = .7, color = "black") +
  facet_grid(spID ~ PA_thresh, switch = "both") +
  theme_void() +
  xlab("Threshold") +
  ylab("Simulated species prevalence") +
  ggtitle("lasso:") +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0, 0, .25, .25), "cm"),
        plot.title = element_text(size = 24))

p.rf = multiThresh %>%
  filter(simnum == 1,
         algorithm == "rf") %>%
  ggplot(., aes(ymin = 0, xmin = 0)) +
  geom_rect(aes(xmax = max.x, ymax = max.y),
            fill = NA, size = .7, color = "black") +
  facet_grid(spID ~ PA_thresh, switch = "both") +
  theme_void() +
  xlab("Threshold") +
  ylab("Simulated species prevalence") +
  ggtitle("rf:") +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0, 0, .25, .25), "cm"),
        plot.title = element_text(size = 24))

p.ridge = multiThresh %>%
  filter(simnum == 1,
         algorithm == "ridge") %>%
  ggplot(., aes(ymin = 0, xmin = 0)) +
  geom_rect(aes(xmax = max.x, ymax = max.y,
                fill = cell),
            size = .7) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(spID ~ PA_thresh, switch = "both") +
  theme_void() +
  xlab("Threshold") +
  ylab("Simulated species prevalence") +
  # ggtitle("ridge:") +
  theme(strip.text.y = element_text(size = 12, angle = 270,
                                    margin = margin(t = 0, r = 4, b = 0, l = 0)),
        strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  margin = margin(t = 6, r = 0, b = 0, l = 0)),
        plot.margin = unit(c(0, 0, .25, .25), "cm"),
        axis.title.y = element_text(size = 14, angle = 90,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank())
ggsave("multiThresh_ridge_8Nov2018.png", width = 12, height = 5)

# pdf("multiThresh_1May2018.pdf", width = 20, height = 16)
# grid.arrange(p.gbm, p.lasso, p.rf, p.ridge, nrow = 2) 
# dev.off()
