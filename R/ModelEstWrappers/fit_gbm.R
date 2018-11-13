###  SDM simulation project

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