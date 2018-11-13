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

####   Random Forest:
fit_rf <- function(df, respvar, covvec) {
  
  time.rf <- suppressWarnings(system.time(model.rf <- randomForest(x = df[df$InEst == 1, covvec],
                                                                   y = df[df$InEst == 1, respvar],
                                                                   ntree = 1000)))
  # suppressing warning because I'm using regression (rather than classification) on binary data
  
  list('model.rf' = model.rf,
       'time.rf' = time.rf)
}
