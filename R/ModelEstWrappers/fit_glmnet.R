###  SDM simulation project

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