###  SDM simulation project

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Each estimation function needs as input:
# df: dataframe with covariates and response; note that only a subset of rows (InEst == 1) should be used for estimation
# respvar: quoted name of response variable
# covvec: character vector of covariates for estimation

####   GBM:
fit_gbm <- function(df, respvar, covvec) {
  time.gbm <- system.time(model.gbm <- gbm.step(data = df[df$InEst == 1, ], 
                                                gbm.x = covvec,
                                                gbm.y = respvar,
                                                family = "bernoulli", 
                                                tree.complexity = 3,
                                                learning.rate = 0.001, 
                                                bag.fraction = 0.5))
  list('model.gbm' = model.gbm,
       'time.gbm' = time.gbm)
}
