####    compare AUCPR to model performance for guiding surveys

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  6 Feb 2018
##  27 Feb 2018: added r to plots; exported L48 plot
##  1 March 2018: updated to new prev; combined West and L48 into 1 fig
##  17 May 2018: updating to version masking Great Lakes; tweaked plot & incorporated JAH suggestions
##  26 July 2018: saved as new file: RareSp_SurveyPerf.R instead of survey_perf.R; updating to just focus on rare sp within West
#                 also now just focus on AUC, aucPR, min-adj aucPR (aucNPR), precision....actually dropping min-adj version b/c adjustment doesn't matter for rare sp
##  6 Aug 2018: adding top_n search
##  13 Nov 2018: added USGS disclaimer

library(tidyverse)
theme_set(theme_bw())
library(gridExtra)
library(modEvA)

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Import model predictions and performance metrics    ####
###

pred.all = read_rds("pred.all_16May2018.rds")
head(pred.all)

simPerf = read_rds("simPerf.all.rds")
head(simPerf)

simPerf = simPerf %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, space,
                AUC, aucPR, precision) 

###
#####    Add calibration statistic     #####
###


calib.west = pred.all %>%
  filter(InEst == 0,
         spID == "sp1",
         x < -100) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont)) %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  summarize(MillerSlope = MillerCalib(obs = TruePA, pred = PredCont, plot = FALSE)$slope) %>%
  mutate(MillerAbsDiff1 = abs(MillerSlope - 1),
         space = "West")
  
simPerf = simPerf %>%
  inner_join(., calib.west)


###
#####   Sample based on continuous prediction    #####
###

head(pred.all)

###   West: sample 250 sites
set.seed(111315)

SurveyRare.west = pred.all %>%
  filter(InEst == 0,
         # just rare sp
         spID == "sp1", 
         # just in the west
         x < -100,
         algorithm != "glmnet") %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
         space = "West",
         method = "weightedSample") %>% 
  sample_n(size = 250, replace = FALSE, weight = PredCont) %>%
  group_by(spID, simnum, sp_sim, algorithm, space, method) %>%
  summarize(n.found = sum(TruePA),
            n.sample = n(),
            prop.found = n.found/n.sample) %>%
  inner_join(., simPerf)  
glimpse(SurveyRare.west)

SurveyRare.west %>% as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()


###
#####   Simple plot for ms revision: just AUC-PR, AUC-ROC, precision    #####
###

# corPerf = SurveyRare.west %>%
#   rename(AUCPR = aucPR) %>%
#   gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
#   group_by(metric) %>%
#   summarise(cor = round(cor(prop.found, value), 2)) 
# corPerf
# 
# 
# p.surveyPerf = SurveyRare.west %>%
#   rename(AUCPR = aucPR) %>%
#   gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
#   # order for plot:
#   ungroup() %>%
#   mutate(metric = factor(metric,
#                          levels = c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
#                          labels = c("AUC-PR", "AUC-ROC",  "Precision", "Lack of calibration"),
#                          ordered = TRUE)) %>%
#   ggplot(., aes(value, prop.found)) +
#   geom_point(size = 3) +
#   facet_wrap( ~ metric, scale = "free_x", ncol = 1, strip.position = "bottom") +
#   ylab("Proportion of sampled\nlocations with presence") +
#   theme(legend.title = element_blank(),
#         axis.text = element_text(color = "black", size = 10),
#         axis.title.y = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         strip.text = element_text(size = 12),
#         strip.placement = "outside",
#         axis.title.x = element_blank(),
#         strip.background = element_blank()) 
# ggsave("surveyPerf_26July2018.png", height = 3, width = 9)
# ggsave("./AUCPR_figures/surveyPerf_26July2018.pdf", height = 3, width = 9)

###
#####    Sample based on top_n      #####
###

###   West: sample top 250 sites
set.seed(56676)

TopN.west = pred.all %>%
  filter(InEst == 0,
         # just rare sp
         spID == "sp1", 
         # just in the west
         x < -100,
         algorithm != "glmnet") %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
         space = "West",
         method = "topN") %>% 
  top_n(250, PredCont) %>%
  group_by(spID, simnum, sp_sim, algorithm, space, method) %>%
  summarize(n.found = sum(TruePA),
            n.sample = n(),
            prop.found = n.found/n.sample) %>%
  inner_join(., simPerf)  
glimpse(TopN.west)

TopN.west %>% as.data.frame() %>%
  filter(!complete.cases(.)) %>%
  nrow()

# TopN.west %>%
#   group_by(spID, simnum, sp_sim, algorithm) %>%
#   rename(AUCPR = aucPR) %>%
#   gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
#   group_by(metric) %>%
#   summarise(cor = round(cor(prop.found, value), 2)) 
# 
# 
# TopN.west %>%
#   rename(AUCPR = aucPR) %>%
#   gather(metric, value, AUC, AUCPR, precision, MillerAbsDiff1) %>%
#   # order for plot:
#   ungroup() %>%
#   mutate(metric = factor(metric,
#                          levels = c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
#                          labels = c("AUC-PR", "AUC-ROC",  "Precision", "Lack of calibration"),
#                          ordered = TRUE)) %>%
#   ggplot(., aes(value, prop.found)) +
#   geom_point(size = 3) +
#   facet_wrap( ~ metric, scale = "free_x", ncol = 4, strip.position = "bottom") +
#   ylab("Proportion of sampled\nlocations with presence") +
#   theme(legend.title = element_blank(),
#         axis.text = element_text(color = "black", size = 10),
#         axis.title.y = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         strip.text = element_text(size = 12),
#         strip.placement = "outside",
#         axis.title.x = element_blank(),
#         strip.background = element_blank()) 


###
#####   Combine and plot     #####
###
head(SurveyRare.west)
head(TopN.west)

TwoTypes = bind_rows(SurveyRare.west, TopN.west)

##  correlation to add to plot:
cor.df = TwoTypes %>%
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
  mutate(value = rep(c(.63, .962, .6, .58), each = 2),
         prop.found = rep(c(.47, .27), times = 4))

# use geom_blank to set axis limits for each panel:
dummy.df <- expand.grid(metric = factor(c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
                                        levels = c("AUCPR", "AUC", "precision", "MillerAbsDiff1"),
                                        labels = c("AUC-PR", "AUC-ROC",  "Precision", "Lack of calibration"),
                                        ordered = TRUE),
                        method = factor(unique(TwoTypes$method),
                                        levels = c("topN", "weightedSample"),
                                        labels = c("Highest-ranked sites", "Weighted sample of sites"),
                                        ordered = TRUE)) %>%
  bind_cols(., data.frame(min_value = rep(c(.54, .942, .15, 0), times = 2),
                          max_value = rep(c(.65, .97, .8, .8), times = 2),
                          min_prop.found = rep(c(.25, .45), each = 4),
                          max_prop.found = rep(c(.55, .95), each = 4))) %>%
  gather(variable, limit, -metric, -method) %>%
  separate(variable, c("type", "variable"), sep = "_") %>%
  spread(variable, limit)

p.Survey = TwoTypes %>%
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
  geom_text(data = cor.df, aes(label = paste0("r = ", cor)),
            size = 6) +
  geom_blank(data = dummy.df) +
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        strip.background = element_blank()) 
ggsave("./AUCPR_figures/SurveyPerf_TopWeighted_16Aug2018.png", width = 10, height = 6)
ggsave("./AUCPR_figures/SurveyPerf_TopWeighted_16Aug2018.pdf", width = 10, height = 6)
