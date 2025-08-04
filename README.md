# Missing-Confounders-Methods
This is a repository of code for the paper "Assessing treatment effects in observational data with missing confounders: A comparative study of practical doubly-robust and traditional missing data methods" currently under development as part of the FDA Sentinel Initiative (https://www.fda.gov/safety/fdas-sentinel-initiative). The paper can be found at https://arxiv.org/abs/2412.15012.

## reports

This directory contains tables of simulation results; some tables are repeated from the main supplement, while others describe simulation scenarios that were not featured in the manuscript.

## sims

This directory contains the code that we used to run simulations, and a script that reproduces one of the synthetic-data results in the manuscript. Please see the README for more details.

## vignettes

This directory contains two vignettes demonstrating how to estimate conditional and marginal parameters discussed in our paper.

1. The file `Marginal EstimationVignette_20240403.pdf` demonstrates how to get marginal estimates and confidence intervals for the IPW, generalized raking, and multiple imputation procedures.
2. The files `tmle_vignette.pdf` and `tmle_vignette.Rmd` demonstrate how to obtain marginal and conditional estimates and confidence intervals using TMLE.

