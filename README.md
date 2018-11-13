# AUC-PR
Code simulates distribution of virtual species, estimates models, and evaluates use of AUC-PR for performance assessment

USGS disclaimer:
This software has been approved for release by the U.S. Geological Survey (USGS). 
Although the software has been subjected to rigorous review, 
the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
of the software and related material nor shall the fact of release constitute any such warranty. 
Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
shall be held liable for any damages resulting from its authorized or unauthorized use.

Overview of code files:
* VirtualSpCreate_3Prev: get climate data & simulate species (continuous suitability)
* Plot_Sim3Prev: create Fig 3 and Fig S1
* threshSampleSim_3Prev: function to threshold continuous suitability, sample observations, and extract covariates
* fitMultiSim_3Prev: function to fit multiple model algorithms using df with list columns. Relies on model fitting functions in ‘ModelEstWrappers’ folder
* SimulateAndFit_3Prev_script: calls threshSampleSim_3Prev & fitMultiSim_3Prev. Fits models and generates predictions
* Perf_3PrevSim: performance metric computation and comparison, including extent change
* RareSp_SurveyPerf: compare perf metrics to survey performance
* plotAUCPR: create Fig 1, S3, S5
* MESS_ PlotandSensitivity: sensitivity analysis excluding areas of extrapolation in prediction
* WorkflowFigure: for Fig 2
* PresBackground_AUCPR: sensitivity analysis with presence background data
