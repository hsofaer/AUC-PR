###  SDM simulation project

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