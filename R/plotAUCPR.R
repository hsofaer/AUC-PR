##  Plot AUCPR:
# comparison with AUC
# aucPR vs minimum-adjusted AUCPR

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  March 2018
##  16 May 2018: integrated Jennifer's updates to plots
##  17 May 2018: updated to version with Great Lakes clipped
##  8 Aug 2018: add comparison between West and L48
##  13 Nov 2018: added USGS disclaimer

library(tidyverse)
theme_set(theme_bw())
library(gridExtra)
library(PRROC)
library(viridis)

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Import performance metrics    ####
###

simPerf = read_rds("simPerf.all.rds")
head(simPerf)


###  and predictions:
pred.all = read_rds("pred.all_16May2018.rds")
head(pred.all)

###
#####    AUC data      ####
###

plot(simPerf$ROC.obj[[1]])
head(simPerf$ROC.obj[[1]]$curve)

renameROC <- function(ROCobj) {
  curve <- ROCobj$curve
  colnames(curve) <- c("FPR", "Sensitivity", "threshold")
  return(as.tibble(curve))
}

ROCdata = simPerf %>%
  ungroup() %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, space, ROC.obj) %>%
  mutate(curve = map(ROC.obj, renameROC)) %>%
  dplyr::select(-ROC.obj) %>%
  unnest()
head(ROCdata)

ROCdata %>% filter(!complete.cases(.)) %>%
  nrow()


###
#####    AUCPR data     ####
###

plot(simPerf$PR.obj[[1]])
str(simPerf$PR.obj[[1]]$curve)

renamePR <- function(PRobj) {
  curve <- PRobj$curve
  colnames(curve) <- c("Recall", "Precision", "threshold")
  return(as.tibble(curve))
}

PRdata = simPerf %>%
  ungroup() %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, space, PR.obj) %>%
  mutate(curve = map(PR.obj, renamePR)) %>%
  dplyr::select(-PR.obj) %>%
  unnest()
head(PRdata)

PRdata %>% filter(!complete.cases(.)) %>%
  nrow()


###
####   Decide which algorithm to show:    ####
###

alg.a = ROCdata %>%
  filter(algorithm != "glmnet",
         space == "West",
         simnum == 1) %>%
  ggplot(., aes(FPR, Sensitivity)) +
  geom_point(aes(color = spID)) +
  xlab("False positive rate (1 - specificity)") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10, color = "black")) +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~ algorithm)


# compare different algorithms:
alg.b = PRdata %>%
  filter(algorithm != "glmnet", 
         space == "West",
         simnum == 1) %>%
  ggplot(., aes(Recall, Precision)) +
  geom_point(aes(color = spID)) +
  xlab("Recall (sensitivity)") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10, color = "black")) +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~ algorithm)

png("AlgorithmByAUCandAUCPR.png", width = 10, height = 6, units = "in", res = 300)
grid.arrange(alg.a, alg.b, ncol = 2)
dev.off()


###
######    Plot AUC and AUCPR     #####
###

p.ROC = ROCdata %>%
  filter(algorithm == "ridge",
         space == "West",
         simnum == 1) %>%
  ggplot(., aes(FPR, Sensitivity)) +
  geom_point(aes(color = threshold)) +
  xlab("False positive rate (1 - specificity)") +
  scale_color_viridis(option = "inferno", 
                      name = "Threshold") +
  theme(axis.title = element_text(size = 16),
        legend.position = c(.75, .25),
        legend.justification = "center",
        axis.text = element_text(size = 10, color = "black")) +
  annotate("text", x = .15, y = .97, label = "Low") +
  annotate("text", x = .27, y = .85, label = "Medium") +
  annotate("text", x = .55, y = .75, label = "High prevalence") +
  xlim(0, 1) +
  ylim(0, 1) +
  ggtitle("a")

p.PR = PRdata %>%
  filter(algorithm == "ridge", 
         space == "West",
         simnum == 1) %>%
  ggplot(., aes(Recall, Precision)) +
  geom_point(aes(color = threshold)) +
  xlab("Recall (sensitivity)") +
  scale_color_viridis(option = "inferno",
                      name = "Threshold") +
  theme(axis.title = element_text(size = 16),
        legend.position = c(.25,  .25),
        legend.justification = "center",
        axis.text = element_text(size = 10, color = "black")) +
  annotate("text", x = .9, y = .125, label = "Low") +
  annotate("text", x = .92, y = .44, label = "Medium") +
  annotate("text", x = .9, y = .71, label = "High prevalence") +
  xlim(0, 1) +
  ylim(0, 1) +
  ggtitle("b")

png("Fig1_AUCvsAUCPR.png", width = 10, height = 5, units = "in", res = 500)
grid.arrange(p.ROC, p.PR, ncol = 2)
dev.off()

pdf("./AUCPR_figures/Fig1_AUCvsAUCPR.pdf", width = 10, height = 5)
grid.arrange(p.ROC, p.PR, ncol = 2)
dev.off()

###
#####    AUCpr vs minimum-adjusted AUCPR    #####
###

simPerf %>%
  ungroup() %>%
  filter(simnum == 1) %>%
  dplyr::select(spID, space, TruePrev) %>%
  mutate(TruePrev = round(TruePrev, 2)) %>%
  distinct()

p.MinAdj = simPerf %>%
  filter(algorithm != "glmnet",
         simnum == 1) %>%
  dplyr::select(spID, simnum, sp_sim, algorithm, space, TruePrev, aucPR, aucNPR) %>%
  unite(sp_space, spID, space, sep = "_", remove = FALSE) %>%
  mutate(sp_space = factor(sp_space,
                           levels = c("sp1_West", "sp1_L48", "sp2_West", "sp2_L48", "sp3_West", "sp3_L48"),
                           labels = c("Low West: prevalence = 0.05",
                                      "Low L48: prevalence = 0.03",
                                      "Medium West: prevalence = 0.25",
                                      "Medium L48: prevalence = 0.13", 
                                      "High West: prevalence = 0.5",
                                      "High L48: prevalence = 0.32"),
                           ordered = TRUE)) %>%
  ggplot(., aes(aucPR, aucNPR)) +
  geom_point(aes(color = sp_space,
                shape = sp_space),
             size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = rep(c("yellow3", "darkgoldenrod3", "brown"), each = 2)) +
  scale_shape_manual(values = rep(c(16, 17), times = 3)) +
  xlab("AUC-PR") +
  ylab("Minimum-adjusted AUC-PR") +
  scale_x_continuous(breaks = seq(.45, .8, by = .05),
                     labels = seq(.45, .8, by = .05),
                     limits = c(.45, .8)) +
  scale_y_continuous(breaks = seq(.45, .8, by = .05),
                     labels = seq(.45, .8, by = .05),
                     limits = c(.45, .8)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(.025, .86),
        legend.justification = "left",
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        legend.key.size = unit(.35, "cm")) +
  coord_equal(ratio = 1)
ggsave("SuppFig_minAdj_17May2018.png", width = 5, height = 5)


###
#####   S5: Plots comparing West and L48       #####
###

##  boxplot of continuous predictions vs truth:
pred.W = pred.all %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp1",
         InEst == 0,
         x < -100) %>%
  mutate(space = "West")

pred.L48 = pred.all %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp1",
         InEst == 0) %>%
  mutate(space = "L48")

p.a = bind_rows(pred.W, pred.L48) %>%
  mutate(space = factor(space,
                        levels = c("West", "L48"),
                        ordered = TRUE)) %>%
  ggplot(., aes(space, PredCont, 
                fill = factor(TruePA),
                color = factor(TruePA))) +
  geom_violin(adjust = 6) +
  ylim(0, 1) +
  ylab("Continuous predicted suitability") +
  scale_fill_manual(values = c("purple4", "darkgoldenrod"),
                    labels = c("True absence", "True presence"),
                    name = NULL) +
  scale_color_manual(values = c("purple4", "darkgoldenrod"),
                    labels = c("True absence", "True presence"),
                    name = NULL) +
  theme_classic() +
  ggtitle("a") +
  theme(legend.position = c(.2, .9),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18)) 

  
p.b = ROCdata %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp1") %>%
  ggplot(., aes(FPR, Sensitivity, group = space)) +
  geom_line(aes(color = space, size = space, alpha = space)) +
  scale_color_manual(values = c("blue4", "cornflowerblue")) +
  scale_alpha_manual(values = c(1, .4)) +
  scale_size_manual(values = c(1, 3)) +
  xlab("False positive rate (1 - specificity)") +
  theme(axis.title = element_text(size = 16),
        legend.position = c(.75, .25),
        legend.justification = "center",
        axis.text = element_text(size = 10, color = "black")) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  ggtitle("c") +
  theme(legend.position = c(.75, .4),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 18)) 

p.c = PRdata %>%
  filter(algorithm == "ridge", 
         simnum == 1,
         spID == "sp1") %>%
  ggplot(., aes(Recall, Precision, group = space)) +
  geom_line(aes(color = space, size = space, alpha = space)) +
  scale_color_manual(values = c("blue4", "cornflowerblue")) +
  scale_alpha_manual(values = c(1, .4)) +
  scale_size_manual(values = c(1, 3)) +
  xlab("Recall (sensitivity)") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  ggtitle("e") +
  theme(legend.position = c(.25, .25),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 18)) 

# pdf("./AUCPR_figures/PRvsROC_WestL48_8Aug2018.pdf", width = 5, height = 11)
# grid.arrange(p.a, p.b, p.c, ncol = 1)
# dev.off()
# 
# png("./AUCPR_figures/PRvsROC_WestL48_8Aug2018.png", width = 5, height = 11,
#     units = "in", res = 300)
# grid.arrange(p.a, p.b, p.c, ncol = 1)
# dev.off()

###   Do the same for the common species:

pred.W.med = pred.all %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp2",
         InEst == 0,
         x < -100) %>%
  mutate(space = "West")

pred.L48.med = pred.all %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp2",
         InEst == 0) %>%
  mutate(space = "L48")

p.violin.med = bind_rows(pred.W.med, pred.L48.med) %>%
  mutate(space = factor(space,
                        levels = c("West", "L48"),
                        ordered = TRUE)) %>%
  ggplot(., aes(space, PredCont, 
                fill = factor(TruePA),
                color = factor(TruePA))) +
  geom_violin(adjust = 2) +
  ylim(0, 1) +
  ylab("Continuous predicted suitability") +
  scale_fill_manual(values = c("purple4", "darkgoldenrod"),
                    labels = c("True absence", "True presence"),
                    name = NULL) +
  scale_color_manual(values = c("purple4", "darkgoldenrod"),
                     labels = c("True absence", "True presence"),
                     name = NULL) +
  theme_classic() +
  ggtitle("b") +
  theme(legend.position = 'none',
        legend.text = element_text(size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18)) 


p.roc.med = ROCdata %>%
  filter(algorithm == "ridge",
         simnum == 1,
         spID == "sp2") %>%
  ggplot(., aes(FPR, Sensitivity, group = space)) +
  geom_line(aes(color = space, size = space, alpha = space)) +
  scale_color_manual(values = c("blue4", "cornflowerblue")) +
  scale_alpha_manual(values = c(1, .4)) +
  scale_size_manual(values = c(1, 3)) +
  xlab("False positive rate (1 - specificity)") +
  theme(axis.title = element_text(size = 16),
        legend.position = c(.75, .25),
        legend.justification = "center",
        axis.text = element_text(size = 10, color = "black")) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  ggtitle("d") +
  theme(legend.position = c(.75, .4),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 18)) 

p.pr.med = PRdata %>%
  filter(algorithm == "ridge", 
         simnum == 1,
         spID == "sp2") %>%
  ggplot(., aes(Recall, Precision, group = space)) +
  geom_line(aes(color = space, size = space, alpha = space)) +
  scale_color_manual(values = c("blue4", "cornflowerblue")) +
  scale_alpha_manual(values = c(1, .4)) +
  scale_size_manual(values = c(1, 3)) +
  xlab("Recall (sensitivity)") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  ggtitle("f") +
  theme(legend.position = c(.25, .25),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 18)) 


pdf("./AUCPR_figures/PRvsROC_2sp_WestL48_18Aug2018.pdf", width = 10, height = 11)
grid.arrange(p.a, p.violin.med, 
             p.b, p.roc.med,
             p.c, p.pr.med,
             ncol = 2)
dev.off()

png("./AUCPR_figures/PRvsROC_WestL48_18Aug2018.png", width = 10, height = 11,
    units = "in", res = 300)
grid.arrange(p.a, p.violin.med, 
             p.b, p.roc.med,
             p.c, p.pr.med,
             ncol = 2)
dev.off()


# change in AUC-ROC and AUC-PR associated with these plots:
simPerf %>%
  filter(algorithm == "ridge", 
         simnum == 1,
         spID == "sp1") %>%
  dplyr::select(AUC, aucPR, space)

simPerf %>%
  filter(algorithm == "ridge", 
         simnum == 1,
         spID == "sp2") %>%
  dplyr::select(AUC, aucPR, space)
