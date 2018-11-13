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
# glmtype: either "multi" or "regular"

####   GLM:
# using glmulti for model selection

fit_glm <- function(df, respvar, covvec, glmtype = "regular") {
  
  ###  multi:
  if (glmtype == "multi") {
    # passing character vector of x variables directly
    time.glm <- system.time(model.glm <- glmulti(respvar, covvec,
                                                 data = df[df$InEst == 1, ], 
                                                 method = "g", # genetic algorithm; fastest way for glm 
                                                 level = 1, # only main effects; no interactions
                                                 confsetsize = 10, # keep top 10 models
                                                 fitfunction = "glm", family = binomial))
  } else if (glmtype == "regular") {
    
    ###  regular glm:
    glm.formula <- as.formula(paste(respvar,
                                    paste0(covvec, collapse = " + "),
                                    sep = " ~ ")) 
    
    time.glm <- system.time(model.glm <- glm(glm.formula,
                                             data = df[df$InEst == 1, ], 
                                             family = binomial))

    
  } else stop('glmtype not "multi" or "regular"')
  
  list('model.glm' = model.glm,
       'time.glm' = time.glm)
}
