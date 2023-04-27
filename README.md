# Longitudinal BMI and cancer risk
Source code for the paper "Longitudinal body mass index and cancer risk: a cohort study of 2.6 million Catalan adults"

The code is intended as informative since it cannot be run without access to patients level data. In order to get access to it, you can consider request a new research project to SIDIAP (www.sidiap.org).

The code run in R and it uses the tidy format. 

The code is organized in data preparation and data analysis part.

### 1) Data preparation

#### BMI trajectories imputation
dataPrepTrayectories is used to prepare the whole population before applying multilevel time raster multiple imputation in order to impute the trajectories of BMI.
In this file you need to create definitions of cancer, type 2 diabetes mellitus and hypertension, and the files cancer_def, cvd, hta i dm will be helpful.
Multiple imputation is performed in the script MI and post-processed in join imp file.

#### Time to event population preparation
To create the population for the study, you can use the script dataPrepAnalysis. Before doing that, you will need to define the pre and post-menopausal status, using the scripts menopause_definition and pre-postmeno, and the BMI at index, using the corresponding script.


### 2) Data analysis

#### Main analysis
The main analysis is performed in the descriptive analysis, analysis_linear and analysis_spine scripts. Names are representative of the content of the scripts.

#### Sensitivity analysis
In this folder you will find all the sensitivity analyses performed. 

For any doubt you can write and email to apistillo@idiapjgol.info or through this github.
