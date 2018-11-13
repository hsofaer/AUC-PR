##   Mini-plots for figure providing overview of workflow

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  17 may 2018
##  13 Nov 2018: added USGS disclaimer

library(tidyverse) 
theme_set(theme_bw())
library(viridis)

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Import model predictions     #####
###

pred.all = read_rds("pred.all_16May2018.rds")
head(pred.all)

###
#####   Samples for model estimation:   ####
###

states <- map_data("state")
head(states)

p.EstPts = pred.all %>%
  filter(spID == "sp2",
         simnum == 1,
         algorithm == "ridge",
         InEst == 1) %>%
  mutate(TruePA = factor(TruePA, 
                         levels = c(0, 1),
                         labels = c("absent", "present"))) %>%
  ggplot(.) +
  geom_point(aes(x, y, color = TruePA)) +
  scale_color_manual(values = c("#440154FF", "#FDE725FF")) +
  geom_map(data = states, map = states, aes(map_id = region),
           fill = NA, color = "black") +
  theme_void() +
  xlim(-125, -67) +
  ylim(24, 49.5) +
  theme(legend.position = "none")
ggsave("SamplePts.png", width = 6, height = 4)


###
#####   Predictions at West and at L48   #####
###

p.pred.L48 = pred.all %>%
  filter(spID == "sp2",
         simnum == 1,
         algorithm == "ridge") %>%
  ggplot(.) +
  geom_raster(aes(x, y, fill = PredCont)) +
  scale_fill_viridis(na.value = NA) +
  theme_void() +
  theme(legend.position = "none") +
  geom_map(data = states, map = states, aes(map_id = region),
           fill = NA, color = "black") 
ggsave("Pred_L48.png", width = 6, height = 4)

p.pred.West = pred.all %>%
  filter(spID == "sp2",
         simnum == 1,
         algorithm == "ridge",
         x < -100) %>%
  ggplot(.) +
  geom_raster(aes(x, y, fill = PredCont)) +
  scale_fill_viridis(na.value = NA) +
  theme_void() +
  theme(legend.position = "none") +
  xlim(-125, -67) +
  ylim(24, 49.5) +
  geom_map(data = states, map = states, aes(map_id = region),
           fill = NA, color = "black") 
ggsave("Pred_West.png", width = 6, height = 4)


###
####   Plot of new survey locations    ####
###

set.seed(111315)

p.surveyLoc = pred.all %>%
  filter(InEst == 0,
         x < -100,
         # subset for plot:
         algorithm == "ridge",
         simnum == 1,
         spID == "sp2") %>%
  group_by(spID, simnum, sp_sim, algorithm) %>%
  mutate(PredCont = ifelse(PredCont < 0, 1e-15, PredCont), # a few tiny neg values
         space = "West") %>% 
  sample_n(size = 250, replace = FALSE, weight = PredCont) %>%
  ggplot(.) +
  geom_point(aes(x, y, fill = factor(TruePA)), shape = 21) +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  geom_map(data = states, map = states, aes(map_id = region),
           fill = NA, color = "black") +
  theme_void() +
  xlim(-125, -67) +
  ylim(24, 49.5) +
  theme(legend.position = "none")
ggsave("SurveyLoc.png", width = 6, height = 4)
