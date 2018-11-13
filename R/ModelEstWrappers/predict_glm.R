###  SDM simulation project

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Each prediction function needs as input: model object, newdf, covvec
# model: model object output from estimation function (and unlisted to separate it from estimation time)
# df: dataframe with covariates and response; note that only a subset of rows (InEst == 1) should be used for estimation
# covvec: character vector of covariates for estimation (and prediction)
# glmtype: either "multi" or "regular"

predict_glm <- function(model, df, covvec, glmtype) {
  
  if (glmtype == "multi") {
    pred.glm <- t(predict(model,
                          select = "all", # multimodel inference based on the top 10 I saved above
                          type = "response", # otherwise it's on logit scale; see predict.glm defaults
                          newdata = df[, covvec])$averages) # averages = the multimodel predictions
    
  } else if (glmtype == "regular") {
    
    pred.glm <- predict(model, newdata = df[, covvec], type = "response")

  } else stop('glmtype not "multi" or "regular"')
  
  return(pred.glm)
}
