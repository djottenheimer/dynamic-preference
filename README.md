# dynamic-preference
Data and code from the 2020 manuscript on satiety-sensitive preference encoding in ventral pallidum.

Raw data from ephys sessions (behavioral and neural spike timestamps) are in the RAWCh and RAWCh_extra files.
RAWCh contains the main sessions from the manuscript.
RAWCh_extra contains 4 extra sessions used in the supplement.
Sessions < 10 are Specific Cues.
Sessions > 10 are Uncertain Outcome.

Raw data from optogenetic sessions are in the ICSSFiles and OptoFiles folders.

MainAnalysis creates PSTHs (behavior and ephys) and fits the GLMs. It produces R, saved as R_choice.
Note: R_choice is too big for GitHub, but you can easily generate it by running MainAnalysis.

DataForModeling prepares data for model fitting and produces NeuronInfo.
fitModels uses NeuronInfo to fit the models, producing os, saved as ModelFits.
fitModelsCV uses NeuronInfo to fit models for cross-validation, saved as ModelFits_CV.

EphysFigures plots data based on R.
ModelingFigure plots data based on ModelFits.
ModelingFigure_extra plots data based on ModelFits_extra.
OptoFigures plots data from the ChR2 sessions.

Note: initially, the GLM predictor was preference instead of time, so all associated variables are named "preference" or "outbypref"