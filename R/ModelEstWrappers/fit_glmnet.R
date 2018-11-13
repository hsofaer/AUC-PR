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

####   GLMNET:
# glmnet wants a matrix for x and y should be a factor

###   glmnet estimation: default alpha = .5
fit_glmnet <- function(df, respvar, covvec, alpha = alpha) {
  
  time.glmnet <- system.time(model.glmnet <- cv.glmnet(x = as.matrix(df[df$InEst == 1, covvec]),
                                                       y = as.factor(df[df$InEst == 1, respvar]),
                                                       family = 'binomial',
                                                       alpha = alpha,
                                                       nfolds = 10,
                                                       standardize = TRUE))
  list('model.glmnet' = model.glmnet, 
       'time.glmnet' = time.glmnet)
}
