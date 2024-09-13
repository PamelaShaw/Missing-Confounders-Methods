# Reproducing numerical experiments

The code in this folder can be used to reproduce the numerical experiments in ["Assessing treatment effects in observational data with missing confounders: A comparative study of practical doubly-robust and traditional missing data methods"](). 

## Overview

The key R files are:
* `00_install_packages.R`: installs the required packages for the experiments. Run this prior to attempting to run any experiments.
* `00_utils.R`: generally useful functions
* `01_generate_data.R`: functions for generating a dataset
* `02_methods_design.R`: design-based estimators (oracle, confounded model, complete-case, IPW, Raking)
* `02_methods_mi.R`: multiple imputation-based estimators (currently only MICE)
* `02_methods_tmle.R`: TMLE-based estimators
* `02_methods_xgb.R`: MI-XGB estimators
* `03_estimate.R`: run a single iteration of the simulation [generate data, run estimator(s), get output]
* `03_true_values.R`: creates the true values of each estimand. Must be run prior to running `05_summarize.R`
* `04_main.R`: set up a given simulation (based on Y, missing, X models and sample size) and run it a set number of times (e.g., 2500). 
* `05_summarize.R`: obtain summary statistics over the simulation replicates

The full set of experiments was originally run on several Windows Virtual machines, each with 24 cores and 512GB of RAM; using R version >= 4.0. The provided `run_simple_outcome_simple_MAR.bat` file allows you to run the code from the Windows command prompt, and we will use this structure below. If you use a different cluster computing environment, you can use the `.bat` file as a guide. If you prefer to run R interactively (e.g., in RStudio), you will need to run `04_main.R` and modify the arguments according to the simulation you wish to reproduce.

## Reproducing results for the synthetic base case simulation: 12% outcome proportion, 40% missing-data proportion, simple outcome model and simple MAR missingness

The provided file `run_simple_outcome_simple_MAR.bat` can be used to reproduce the results for this scenario. Run this code from the `sims/` subdirectory.

Specifically, the `.bat` file sets `yscenario=1.1` and `mscenario=1.1`, then runs the simulation for 2500 replications for each of the following estimators: 
* `cc_oracle`: the oracle estimator, which uses the true data-generating outcome model and can access the confounder data for all observations (i.e., there is no missing data)
* `cc_population`: the best-case census estimator, which can access the confounder data for all observations (i.e., there is no missing data)
* `cc_noW`: the confounded model, which drops the confounders that are prone to missingness
* `cc`: the complete-case estimator, which drops observations with missing values
* `ipw`: the inverse probability weighted estimator
* `gr`: generalized raking
* `mice`: multiple imputation via chained equations (MICE) (using predictive mean matching)
* `rf`: MICE using random forests
* `xgb`: MICE using gradient boosted trees
* `ipcw-tmle_m`: inverse probability of coarsening-weighted (IPCW) targeted maximum likelihood estimation (TMLE) using a Super Learner only for the missing-data model (M); generalized linear models are used for the outcome regression and treatment propensity models
* `ipcw-tmle_mto`: IPCW-TMLE using a Super Learner for the missing-data model, outcome regression model, and treatment propensity model (MTO)
* `ipcw-a-tmle_m`: augmented-data IPCW-TMLE-M, passing a confounded estimate of the outcome regression as an additional data column
* `ipcw-a-tmle_mto`: augmented-data IPCW-TMLE-MTO, passing a confounded estimate of the outcome regression as an additional data column

Once the simulation completes, there will be a folder `sims/Results/m1.1_y1.1_x1_s1/`, which contains two `.csv` files for each estimator: `coefs*.csv` contains the regression coefficients for all variables in the census working logistic regression model, while `output*.csv` contains the point estimates (and other information) for the conditional odds ratio and the marginal risk difference, relative risk, and odds ratio.

Next, run `03_true_values.R`. This will create the true values of each estimand. 

Finally, run `05_summarize.R`. This will read in the results and produce LaTeX tables with the summaries, matching those for this scenario in the main manuscript and supplementary material.

## More detailed description of R files

### Generally-useful functions

`00_utils.R` contains generally useful functions that may be helpful across other `R` scripts. This file gets `source`d into `04_main.R`.

### Data generation

`01_generate_data.R` generates the data. This file gets `source`d into `04_main.R`.

Key inputs here are:
* the outcome regression model (`YScenario`), which can range from a simple glm to more complicated functional form
* the missing-data model (`missScenario`), which can range from simple MAR scenarios (e.g., a glm that doesn't depend on the outcome; a glm that does depend on the outcome) to more complicated scenarios
* the covariate model (`XScenario`), which determines correlation and treatment assignment

The main function is `gen_data`, which generates a dataset of a given sample size `N` based on a given `YScenario`, `missScenario`, `XScenario`, and correlation.

### Methods

`02_methods_*.R` houses the individual methods. These files gets `source`d into `04_main.R`.

Methods are grouped into three groups: design-based methods (e.g., complete-case, IPW, raking; `02_methods_design.R`), multiple imputation-based methods (e.g., MICE; `02_methods_mi.R`), TMLE-based methods (`02_methods_tmle.R`), and MI-XGB (`02_methods_xgb.R`). Each group has its own `02_methods` file. The main function for each method is of the form `run_*`, which runs the method based on a given dataset.

### Estimate parameter(s) of interest

`03_estimate.R` runs the simulation a single time, from data generation through estimation and returning results. Can also be used to return only the dataset, without running any estimators. This file gets `source`d into `04_main.R`.

### Running the simulation

`04_main.R` takes in command-line arguments and runs the simulations. It `source`s all of the previous files and saves output to a directory tree differentiated by Y and missing scenario, followed by estimator. The file name includes Y and missing scenario, sample size, and estimator.

To run multiple estimators at once, pass in a semicolon-separated character string with the multiple estimators (e.g., `"cc_oracle;cc_noW"` will run both the oracle complete-case estimator [using confounders that we don't have access to for other estimators] and the confounded complete-case estimator dropping the columns with missing data).

If you are running `04_main.R` interactively (e.g., in RStudio), you will need to modify the appropriate elements of the `args` list after running through line 82. For example, if you want to run Y scenario 1.1 rather than Y scenario 1, you will need to set `args$yscenario <- 1.1`.

## R `sessionInfo()`

Here is the output of `sessionInfo()` on a representative virtual machine used for the analysis:

```{r}
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] splines   grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] marginaleffects_0.20.0 mixgb_1.0.2            tmle_2.0.1             readr_2.1.5           
 [5] dplyr_1.1.4            future.apply_1.11.2    future_1.33.2          optparse_1.7.5        
 [9] dbarts_0.9-28          ranger_0.16.0          xgboost_1.7.7.1        glmnet_4.1-8          
[13] SuperLearner_2.0-29    gam_1.22-3             foreach_1.5.2          nnls_1.5              
[17] mice_3.16.0            sandwich_3.1-0         survey_4.4-2           survival_3.6-4        
[21] Matrix_1.7-0           pROC_1.18.5            mvtnorm_1.2-4         

loaded via a namespace (and not attached):
 [1] gtable_0.3.5       shape_1.4.6.1      ggplot2_3.5.1      lattice_0.22-6     tzdb_0.4.0        
 [6] vctrs_0.6.5        tools_4.4.0        generics_0.1.3     parallel_4.4.0     getopt_1.20.4     
[11] tibble_3.2.1       fansi_1.0.6        pan_1.9            pkgconfig_2.0.3    jomo_2.7-6        
[16] data.table_1.15.4  RcppParallel_5.1.7 lifecycle_1.0.4    compiler_4.4.0     munsell_0.5.1     
[21] mitools_2.4        codetools_0.2-20   pillar_1.9.0       nloptr_2.0.3       tidyr_1.3.1       
[26] MASS_7.3-60.2      iterators_1.0.14   rpart_4.1.23       boot_1.3-30        mitml_0.4-5       
[31] nlme_3.1-164       parallelly_1.37.1  tidyselect_1.2.1   digest_0.6.35      purrr_1.0.2       
[36] listenv_0.9.1      colorspace_2.1-0   cli_3.6.2          magrittr_2.0.3     Rfast_2.1.0       
[41] utf8_1.2.4         broom_1.0.5        scales_1.3.0       backports_1.4.1    RcppZiggurat_0.1.6
[46] globals_0.16.3     nnet_7.3-19        lme4_1.1-35.3      zoo_1.8-12         hms_1.1.3         
[51] rlang_1.1.3        Rcpp_1.0.12        glue_1.7.0         DBI_1.2.2          rstudioapi_0.16.0 
[56] minqa_1.2.6        jsonlite_1.8.8     R6_2.5.1           plyr_1.8.9        
```