###   SDM simulations - fit multiple types of models in list-df structure

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

##  Helen Sofaer
##  1 Nov 2017
##  21 Nov 2017: for aucPR: edited fitMultiSim.R to drop predictions to future; update column name with estimation data
##  13 Nov 2018: added USGS disclaimer

##  Function inputs:
# SimSp: dataframe with list columns. Each row is one simulation, columns are sim IDs (sp and simnum), sim.data (list)
# covvec: character vector of the columns to include (bio1 through bio19)
# respvar: PA column to model; either "PA.bio1bio12" or "PA.bio6bio12"
# glmtype: "none", "regular" or "multi"; the latter uses glmulti (see fit_glm.R)
# glmnet.alpha: alpha used for elastic net (lasso and ridge have fixed values); default = .5

##  Function output:
# SimSp with additional (list) columns, containing model objects, timings, and predictions

fitMultiSim_3Prev <- function(SimSp, covvec, respvar, glmtype = "regular", glmnet.alpha = .5, ...) {
  
  if (glmtype == "none") {
    
    sim.model <- SimSp %>%
      mutate(
        ##  glmnet estimation and prediction:
        est.glmnet = map(sim.data, 
                         fit_glmnet, 
                         respvar = respvar, 
                         covvec = covvec,
                         alpha = glmnet.alpha),
        # unpack the list within est.glmnet
        fit.glmnet = map(est.glmnet, 'model.glmnet'),
        time.glmnet = map(est.glmnet, 'time.glmnet'),
        # predict:
        PredCurr.glmnet = map2(fit.glmnet, sim.data, predict_glmnet, covvec = covvec),
        ##  lasso estimation and prediction: alpha = 1
        est.lasso = map(sim.data, fit_glmnet, 
                        respvar = respvar, 
                        covvec = covvec,
                        alpha = 1),
        fit.lasso = map(est.lasso, 'model.glmnet'),
        time.lasso = map(est.lasso, 'time.glmnet'),
        PredCurr.lasso = map2(fit.lasso, sim.data, predict_glmnet, covvec = covvec),
        ##  ridge estimation and prediction: alpha = 0
        est.ridge = map(sim.data, fit_glmnet,
                        respvar = respvar, 
                        covvec = covvec,
                        alpha = 0),
        fit.ridge = map(est.ridge, 'model.glmnet'),
        time.ridge = map(est.ridge, 'time.glmnet'),
        PredCurr.ridge = map2(fit.ridge, sim.data, predict_glmnet, covvec = covvec),
        # random forest estimation and prediction:
        est.rf = map(sim.data, fit_rf,
                     respvar = respvar, 
                     covvec = covvec),
        fit.rf = map(est.rf, 'model.rf'),
        time.rf = map(est.rf, 'time.rf'),
        PredCurr.rf = map2(fit.rf, sim.data, predict_rf, covvec = covvec),
        # gbm estimation and prediction:
        est.gbm = map(sim.data, fit_gbm,
                      respvar = respvar, 
                      covvec = covvec),
        fit.gbm = map(est.gbm, 'model.gbm'),
        time.gbm = map(est.gbm, 'time.gbm'),
        PredCurr.gbm = map2(fit.gbm, sim.data, predict_gbm, covvec = covvec)
      ) %>%
      dplyr::select(-est.glmnet, -est.lasso, -est.ridge, -est.rf, -est.gbm)
    
  }  else {
    
    sim.model <- SimSp %>%
      mutate(
        ##  glmnet estimation and prediction:
        est.glmnet = map(sim.data, 
                         fit_glmnet, 
                         respvar = respvar, 
                         covvec = covvec,
                         alpha = glmnet.alpha),
        # unpack the list within est.glmnet
        fit.glmnet = map(est.glmnet, 'model.glmnet'),
        time.glmnet = map(est.glmnet, 'time.glmnet'),
        # predict:
        PredCurr.glmnet = map2(fit.glmnet, sim.data, predict_glmnet, covvec = covvec),
        ##  lasso estimation and prediction: alpha = 1
        est.lasso = map(sim.data, fit_glmnet, 
                        respvar = respvar, 
                        covvec = covvec,
                        alpha = 1),
        fit.lasso = map(est.lasso, 'model.glmnet'),
        time.lasso = map(est.lasso, 'time.glmnet'),
        PredCurr.lasso = map2(fit.lasso, sim.data, predict_glmnet, covvec = covvec),
        ##  ridge estimation and prediction: alpha = 0
        est.ridge = map(sim.data, fit_glmnet,
                        respvar = respvar, 
                        covvec = covvec,
                        alpha = 0),
        fit.ridge = map(est.ridge, 'model.glmnet'),
        time.ridge = map(est.ridge, 'time.glmnet'),
        PredCurr.ridge = map2(fit.ridge, sim.data, predict_glmnet, covvec = covvec),
        # random forest estimation and prediction:
        est.rf = map(sim.data, fit_rf,
                     respvar = respvar, 
                     covvec = covvec),
        fit.rf = map(est.rf, 'model.rf'),
        time.rf = map(est.rf, 'time.rf'),
        PredCurr.rf = map2(fit.rf, sim.data, predict_rf, covvec = covvec),
        # glm estimation and multimodel prediction via glmulti:
        est.glm = map(sim.data, fit_glm,
                      respvar = respvar, 
                      covvec = covvec,
                      glmtype = glmtype),
        fit.glm = map(est.glm, 'model.glm'),
        time.glm = map(est.glm, 'time.glm'),
        PredCurr.glm = map2(fit.glm, sim.data, predict_glm, covvec = covvec, glmtype = glmtype),
        # gbm estimation and prediction:
        est.gbm = map(sim.data, fit_gbm,
                      respvar = respvar, 
                      covvec = covvec),
        fit.gbm = map(est.gbm, 'model.gbm'),
        time.gbm = map(est.gbm, 'time.gbm'),
        PredCurr.gbm = map2(fit.gbm, sim.data, predict_gbm, covvec = covvec)
      ) %>%
      dplyr::select(-est.glmnet, -est.lasso, -est.ridge, -est.rf, -est.glm, -est.gbm)
    
  }
  
  return(sim.model)
  
}