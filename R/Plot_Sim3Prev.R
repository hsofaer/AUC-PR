###   Plot virtual species (3 prev) for AUCPR manuscript

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  29 Nov 2017: preliminary versions
##  13 Feb 2018: updated to manuscript quality
##  1 March 2018: updated after fixing prevalence; updated Fig S1 colors to match 
##  16 May 2018: updated labels; updated to version with Great Lakes clipped
##  13 Nov 2018: add USGS disclaimer; updated from Fig 2 to Fig 3 to match ms

library(raster)
library(tidyverse)
theme_set(theme_bw())
library(grid)
library(gridExtra)
library(viridis)
library(gtable)

setwd("H:/RegularizedClimChange/aucprSim")

###
#####   Fig 3: Three virtual species with different prevalence      ####
###

ThreeSp = brick('ThreeSp_TrueContSuit_16May2018.grd')
plot(ThreeSp)

# ##  apply .5 threshold:
# Thresh.5 = calc(ThreeSp, fun = function(x) ifelse(x > .5, 1, 0))
# names(Thresh.5) <- c("Thresh.5.sp1", "Thresh.5.sp2", "Thresh.5.sp3")

##  apply probabilistic threshold:
set.seed(9876)
ThreshProb = stack(calc(ThreeSp[[1]], fun = function(x) rbinom(1, 1, x)),
                   calc(ThreeSp[[2]], fun = function(x) rbinom(1, 1, x)),
                   calc(ThreeSp[[3]], fun = function(x) rbinom(1, 1, x)))
# warnings not a problem - get NAs outside of boundaries where continuous is NA
plot(ThreshProb)
names(ThreshProb) <- c(paste0("ProbThresh", c(".sp1", ".sp2", ".sp3")))



label.contsp = as_labeller(c("TrueContSuit.sp1" = "Low prevalence",
                             "TrueContSuit.sp2" = "Medium prevalence",
                             "TrueContSuit.sp3" = "High prevalence"))

p.a = as.data.frame(ThreeSp, xy = TRUE, na.rm = TRUE) %>%
  gather(variable, value, -x, -y) %>%
  ggplot(., aes(x, y)) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable, labeller = label.contsp) +
  theme_void() +
  scale_fill_viridis(na.value = NA,
                     name = "Continuous\nprobability\n") +
  theme(strip.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  ylim(24.5, 49.4) +
  annotate("segment", x = -100, xend = -100, y = 24.5, yend = 49.4, colour = "white")

# ##  Thresholded at .5: NOT going to include this
# # viridis hex values from:
# # https://www.thedataschool.co.uk/gwilym-lockwood/viridis-colours-tableau/
# p.b = as.data.frame(Thresh.5, xy = TRUE, na.rm = TRUE) %>%
#   gather(variable, value, -x, -y) %>%
#   ggplot(., aes(x, y)) +
#   geom_tile(aes(fill = factor(value))) +
#   facet_wrap(~ variable) +
#   theme_void() +
#   scale_fill_manual(values = c("#440154FF", "#FDE725FF"), 
#                     na.value = NA,
#                     name = "Continuous\nprobability\n>0.5\n") +
#   theme(strip.text = element_blank(),
#         legend.title = element_text(size = 18),
#         legend.text = element_text(size = 14)) +
#   ylim(24.5, 49.4)


p.b = as.data.frame(ThreshProb, xy = TRUE, na.rm = TRUE) %>%
  gather(variable, value, -x, -y) %>%
  mutate(value = factor(value, 
                        levels = c(0, 1),
                        labels = c(" absent", " present"))) %>%
  ggplot(., aes(x, y)) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  theme_void() +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF"), 
                    na.value = NA,
                    name = "Binary\nstatus") +
  theme(strip.text = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.justification = "left") +
  ylim(24.5, 49.4) +
  annotate("segment", x = -100, xend = -100, y = 24.5, yend = 49.4, colour = "white")

##  alignment with gtable (b/c legends are different sizes):
# following: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
g.a <- ggplotGrob(p.a)
g.b <- ggplotGrob(p.b)
g <- rbind(g.a, g.b, size = "first")
g$widths <- unit.pmax(g.a$widths, g.b$widths)

png("Fig3_Map3sp_ContProb.png", width = 18, height = 8, units = 'in', res = 500)
grid.draw(g)
dev.off()

# # For putting into word:
# png("Fig2_Map3sp_ContProb.png", width = 18, height = 8, units = 'in', res = 500)
# grid.arrange(p.a, p.b, ncol = 1)
# dev.off()

# For publication:
pdf("./AUCPR_figures/Fig3_Map3sp_ContProb.pdf", width = 18, height = 8)
grid.draw(g)
dev.off()

###
######   Fig. S1: Plot response curves that generated the species     ####
###

##  Get range of the variables:
curr.bio = read_rds("curr.bio.rds")
varRanges = cellStats(curr.bio, range) 

##  Create a dataframe with all the values, and rescale them
VirtSpResp = data_frame(bio1seq = seq(-4.5, 25, length.out = 150),
                        bio1.sp1 = dnorm(bio1seq, mean = 20, sd = 8),
                        bio1.sp1.rsc = (bio1.sp1 - min(bio1.sp1)) / (max(bio1.sp1) - min(bio1.sp1)),
                        bio1.sp2 = dnorm(bio1seq, mean = 20, sd = 20),
                        bio1.sp2.rsc = (bio1.sp2 - min(bio1.sp2)) / (max(bio1.sp2) - min(bio1.sp2)),
                        bio1.sp3 = dnorm(bio1seq, mean = 20, sd = 46),
                        bio1.sp3.rsc = (bio1.sp3 - min(bio1.sp3)) / (max(bio1.sp3) - min(bio1.sp3)),
                        bio5seq = seq(10.7, 45.4, length.out = 150),
                        bio5.sp1 = dnorm(bio5seq, mean = 45, sd = 8),
                        bio5.sp1.rsc = (bio5.sp1 - min(bio5.sp1)) / (max(bio5.sp1) - min(bio5.sp1)),
                        bio5.sp2 = dnorm(bio5seq, mean = 45, sd = 20),
                        bio5.sp2.rsc = (bio5.sp2 - min(bio5.sp2)) / (max(bio5.sp2) - min(bio5.sp2)),
                        bio5.sp3 = dnorm(bio5seq, mean = 45, sd = 46),
                        bio5.sp3.rsc = (bio5.sp3 - min(bio5.sp3)) / (max(bio5.sp3) - min(bio5.sp3)),
                        bio12seq = seq(47, 3328, length.out = 150),
                        bio12.sp1 = dnorm(bio12seq, mean = 150, sd = 77),
                        bio12.sp1.rsc = (bio12.sp1 - min(bio12.sp1)) / (max(bio12.sp1) - min(bio12.sp1)),
                        bio12.sp2 = dnorm(bio12seq, mean = 150, sd = 173),
                        bio12.sp2.rsc = (bio12.sp2 - min(bio12.sp2)) / (max(bio12.sp2) - min(bio12.sp2)),
                        bio12.sp3 = dnorm(bio12seq, mean = 150, sd = 460),
                        bio12.sp3.rsc = (bio12.sp3 - min(bio12.sp3)) / (max(bio12.sp3) - min(bio12.sp3)))
summary(VirtSpResp)

###   Plot:
png("FigS1_VirtSp_RespCurves.png", width = 10, height = 4, units = "in", res = 600)
par(las = 1, mfrow = c(1, 3), mgp = c(2.25, .75, 0), mar = c(4, 4, .5, .5), bty = "l")
# bio1:
with(VirtSpResp, plot(bio1seq, bio1.sp1.rsc, type = "l", 
                      xlab = "Mean annual temperature (C)", ylab =  "Suitability",
                      col = "yellow3", cex.lab = 1.5))
lines(VirtSpResp$bio1seq, VirtSpResp$bio1.sp2.rsc, col = "darkgoldenrod3", type = "l")
lines(VirtSpResp$bio1seq, VirtSpResp$bio1.sp3.rsc, col = "brown", type = "l")

legend(10, .4, 
       legend = c("Low", "Medium", "High"), 
       lty = 1, 
       col = c("yellow3", "darkgoldenrod3", "brown"),
       bty = 'n',
       cex = 1.25,
       title = "Prevalence")
mtext('a', side = 3, line = -1, las = 1, adj = -.15)

# bio5:
with(VirtSpResp, plot(bio5seq, bio5.sp1.rsc, type = "l", 
                      xlab = "Max. temperature warmest month (C)", ylab =  "",
                      col = "yellow3", cex.lab = 1.5))
lines(VirtSpResp$bio5seq, VirtSpResp$bio5.sp2.rsc, col = "darkgoldenrod3", type = "l")
lines(VirtSpResp$bio5seq, VirtSpResp$bio5.sp3.rsc, col = "brown", type = "l")
mtext('b', side = 3, line = -1, las = 1, adj = -.15)

# bio12:
with(VirtSpResp, plot(bio12seq, bio12.sp1.rsc, type = "l", 
                      xlab = "Annual precipitation (mm)", ylab =  "",
                      col = "yellow3", cex.lab = 1.5))
lines(VirtSpResp$bio12seq, VirtSpResp$bio12.sp2.rsc, col = "darkgoldenrod3", type = "l")
lines(VirtSpResp$bio12seq, VirtSpResp$bio12.sp3.rsc, col = "brown", type = "l")
mtext('c', side = 3, line = -1, las = 1, adj = -.15)

dev.off()