###  SDM simulation project

##  Each prediction function needs as input: model object, newdf, covvec
# model: model object output from estimation function (and unlisted to separate it from estimation time)
# df: dataframe with covariates and response; note that only a subset of rows (InEst == 1) should be used for estimation
# covvec: character vector of covariates for estimation (and prediction)

predict_rf <- function(model, df, covvec) {
  predict(model,
          newdata = df[, covvec])
}
