# Longitudinal BMI and cancer risk
Source code for the paper "Longitudinal body mass index and cancer risk: a cohort study of 2.6 million Catalan adults"

The code is intended as informative since it cannot be run without access to patients level data. In order to get access to it, you can consider request a new research project to SIDIAP (www.sidiap.org).

The code is organized in data preparation and data analysis part.

### 1) Data preparation

#### BMI trajectories imputation
dataPrepTrayectories is used to prepare the whole population before applying multilevel time raster multiple imputation in order to impute the trajectories of BMI.
In this file you need to create definitions of cancer, type 2 diabetes mellitus and hypertension, and the files cancer_def, cvd, hta i dm will be helpful.
Multiple imputation is performed in the script MI and post-processed in join imp file.

#### Time to event population preparation



### 2) Data analysis
