###   MESS plots for AUC-PR paper

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  1 Aug 2018
##  13 Nov 2018: added USGS disclaimer

library(tidyverse)
library(raster)
library(dismo)
library(PRROC)
library(SDMTools)
library(gridExtra)

setwd("H:/RegularizedClimChange/aucprSim")


### 
#####    Import data and bioclim vars    #####
###

# created in SimulateAndFit_3Prev_script.R
pred.all = readRDS("pred.all_16May2018.rds")
fit.all = readRDS("fit.all.rds")
head(pred.all)
head(fit.all)
head(fit.all$sim.data[[1]])

bio19 = read_rds("curr.bio.all19.rds")


###
#####    Compute mess surface     #####
###

MESS.df <- fit.all %>%
  dplyr::select(sp_sim, sim.data) %>%
  distinct() %>%
  mutate(est.data = map(sim.data, filter, InEst == 1),
         bio.vars = map(est.data, magrittr::extract, paste0('bio', seq(1:19))),
         bio19.list = list(bio19),
         # calculate MESS
         MESS = map2(bio19.list, bio.vars, mess, full = FALSE),
         # add MESS variable to full dataset:
         sim.data = map2(MESS, sim.data, ~ mutate(.y, MESS = raster::extract(.x, .y[, c("x", "y")]))))
head(MESS.df)
head(MESS.df$est.data[[1]])
head(MESS.df$bio.vars[[1]])
plot(MESS.df$MESS[[1]])
head(MESS.df$sim.data[[1]])

# check:
as.data.frame(MESS.df$MESS[[1]], xy = TRUE) %>%
  filter(is.finite(mess)) %>%
  ggplot(.) +
  geom_raster(aes(x, y, fill = mess)) +
  scale_fill_gradient2() 

ggplot() +
  geom_point(data = MESS.df$sim.data[[1]],
             aes(x, y, color = MESS)) +
  scale_color_gradient2()
# great.

plot(stack(MESS.df$MESS))
# all very similar

###
######   Join MESS values to predictions     #####
###

head(pred.all)

MESS.xy <- MESS.df %>%
  dplyr::select(sp_sim, sim.data) %>%
  unnest() %>%
  dplyr::select(sp_sim, x, y, MESS) 
head(MESS.xy)

pred.all <- pred.all %>%
  inner_join(., MESS.xy)

###
#####    Calculate performance w/o extrapolation    #####
###

PerfWest.mess = pred.all %>%
  # West only:
  filter(InEst == 0,
         x < -100,
         # exclude negative MESS values:
         MESS >= 0) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE))) %>%
  mutate(space = "West",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich")) 
head(PerfWest.mess)

PerfL48.mess = pred.all %>%
  filter(InEst == 0,
         # exclude negative MESS values:
         MESS >= 0) %>%
  group_by(spID, simnum, sp_sim, algorithm, optmaxSS.Est) %>%
  summarize(AUC = SDMTools::auc(TruePA, PredCont), 
            PR.obj = list(pr.curve(PredCont, weights.class0 = TruePA, 
                                   curve = TRUE, 
                                   max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE))) %>%
  mutate(space = "L48",
         aucPR = map_dbl(PR.obj, "auc.davis.goadrich")) 
head(PerfL48.mess)

Perf.mess <- bind_rows(PerfWest.mess, PerfL48.mess)

###
######    Plot for ms supplement      ####
###

p.MESS = pred.all %>%
  filter(sp_sim == "sp1_sim1") %>%
  dplyr::select(sp_sim, x, y, MESS) %>%
  distinct() %>%
  mutate(alpha = factor(ifelse(MESS < 0, 0, 1))) %>%
  ggplot(.) +
  geom_raster(aes(x, y, fill = MESS,
                  alpha = alpha)) +
  scale_fill_gradient2(mid = "grey80") +
  scale_alpha_manual(values = c(.3, 1),
                     guide = FALSE) +
  theme_void() +
  geom_vline(xintercept = -100, color = "white", size = 1) +
  annotate('text', x = -124.5, y = 49.5, label = 'a', size = 8) +
  theme(legend.position = c(.93, .45))


# plot auc and aucPR change with changing extent, excluding extraolation
p.GeogExt.mess = Perf.mess %>%
  ungroup() %>%
  # just the algorithms I'm using:
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
  scale_y_continuous(breaks = seq(.55, .8, by = .05),
                     limits = c(.55, .77)) +
  scale_x_continuous(breaks = seq(.7, 1, by = .05),
                     limits = c(.7, 1)) +
  theme_classic() +
  theme(legend.position = c(.7, .8),
        legend.title = element_text(size = 18),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10)),
        axis.title.y = element_text(margin = ggplot2::margin(r = 10))) +
  annotate(geom = "text", label = 'b', x = .7, y = .77, size = 8)

png("./AUCPR_figures/Supp_MESS_sensitivity.png", 
    width = 6, height = 8, units = 'in', res = 500)
grid.arrange(p.MESS, p.GeogExt.mess, ncol = 1)
dev.off()

###  Look quickly at one other:
pred.all %>%
  filter(sp_sim == "sp3_sim1") %>%
  dplyr::select(sp_sim, x, y, MESS) %>%
  distinct() %>%
  mutate(alpha = factor(ifelse(MESS < 0, 0, 1))) %>%
  ggplot(.) +
  geom_raster(aes(x, y, fill = MESS,
                  alpha = alpha)) +
  scale_fill_gradient2(mid = "grey80") +
  scale_alpha_manual(values = c(.3, 1),
                     guide = FALSE) +
  theme_void() +
  geom_vline(xintercept = -100, color = "white", size = 1) +
  annotate('text', x = -124.5, y = 49.5, label = 'a', size = 8) +
  theme(legend.position = c(.93, .45))
# really similar