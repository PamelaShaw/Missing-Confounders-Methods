# Methods - design-based

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 20230725       bdw       harmonize methods into a single file
# 20230802       bdw       break out into design-based, MI, TMLE files
# 20230807       bdw       harmonize with latest version from Pam
# 20230816       bdw       add marginal parameter point estimates
# 20230912       bdw       change return to return all coefs
# 20231108       cak       using marginaleffects package for marginal estimates and SEs
# 20240306       bdw       use new version (>= 4.4) of survey

# complete-case analysis -------------------------------------------------------
# @param data the dataset
# @param formula a character string specifying the formula for the model
# @param procedure the CC procedure ("oracle" = fully-observed everything [unobserved], "noW" = confounded model, leaving out covariates subject to missingness, "CC" = only considering participants with )
# @param coefficient_of_interest the coefficient of interest, either a character string (e.g., "X") or its position in a formula (e.g., 2, behind the intercept)
# @param fam outcome regression family
# @return a list with the parameter estimate and standard error (in a data.frame) and the fitted model
run_cc <- function(data, formula = "Y ~ X", procedure = "oracle", coefficient_of_interest = "X", fam = "binomial") {
  if (fam != "survival") {
    fit <- glm(formula = as.formula(formula), family = fam, data = data)
  }
  est <- coef(fit)[coefficient_of_interest] # don't hard code this, get name of variable that we want
  SE <- summary(fit)$coef[coefficient_of_interest, 2]
  conditional_est <- data.frame("procedure" = procedure, "estimand" = "canonical_parameter", 
                                "est" = est, "SE" = SE)
  
  # marginal_params <- get_marginal_estimates(fit = fit, data = data, estimator = procedure)
  marginal_params <- marginal_effects_new(fit=fit, weights=NULL, coefficient_of_interest = coefficient_of_interest)
  
  marginal_est <- data.frame("procedure" = procedure, marginal_params)
  
  coefs <- data.frame("procedure" = procedure, summary(fit)$coefficients)
  names(coefs) <- c("procedure", "Estimate", "SE", "test_statistic", "p_value")
  return(list("results" = rbind.data.frame(conditional_est, marginal_est),
              "coefs" = cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL)))
  )
}

# IPW analysis -----------------------------------------------------------------
# @param data the dataset (assumed to have a column "is.complete" that identifies complete (1) vs missing (0) observations)
# @param formula a character string specifying the formula for the model
# @param coefficient_of_interest the coefficient of interest, either a character string (e.g., "X") or its position in a formula (e.g., 2, behind the intercept)
# @param missing_indicator the name of the missingness indicator
# @param outcome regression family
# @return a list with the parameter estimate(s), procedure, estimand, and standard error(s) (in a data.frame) and the fitted model
run_ipw <- function(data, outcome_formula = "Y ~ X", miss_formula = "is.complete ~ X + Y", coefficient_of_interest = "X", missing_indicator = "is.complete", fam = "binomial") {
  # estimate propensity score, compute IP weights
  miss_fit <- glm(formula = as.formula(miss_formula), family = "binomial", data = data)
  propensity_score <- predict(miss_fit, type = "response")
  ip_weights <- 1 / propensity_score

  # get IPW estimator and robust SE
  data_cc <- data
  data_cc$ipw_wt <- ip_weights
  data_cc <- data_cc[data_cc[[missing_indicator]] == 1, ]
  if (fam == "binomial") {
    ipw_fit <- glm(formula = as.formula(outcome_formula), family = "quasibinomial", data = data_cc, weights = data_cc$ipw_wt)
  }
  est <- coef(ipw_fit)[coefficient_of_interest]
  SE <- sqrt(diag(sandwich::sandwich(ipw_fit)))[coefficient_of_interest]
  # marginal_params <- get_marginal_estimates(fit = ipw_fit, data = data_cc, estimator = "ipw")
  marginal_params <- marginal_effects_new(fit=ipw_fit, weights=weights(ipw_fit), coefficient_of_interest = coefficient_of_interest)
  
  conditional_est <- data.frame("procedure" = "IPW", "estimand" = "canonical_parameter", 
                                     "est" = est, "SE" = SE)
  marginal_est <- data.frame("procedure" = "IPW", marginal_params)
  coefs <- data.frame("procedure" = "IPW", summary(ipw_fit)$coefficients)
  names(coefs) <- c("procedure", "Estimate", "SE", "test_statistic", "p_value")
  return(list("results" = rbind.data.frame(conditional_est, marginal_est),
              "coefs" = cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL))))
}

# @param data the dataset (assumed to have a column "is.complete" that identifies complete (1) vs missing (0) observations)
# @param formula a character string specifying the formula for the model
# @param stable should we use stabilized weights?
# @param coefficient_of_interest the coefficient of interest, either a character string (e.g., "X") or its position in a formula (e.g., 2, behind the intercept)
# @return a list with the parameter estimate(s), procedure, estimand, and standard error(s) (in a data.frame) and the fitted model
run_iptw <- function(data, outcome_formula = "Y ~ X", miss_formula = "is.complete ~ X + Y", stable = TRUE, coefficient_of_interest = "X") {
  # estimate propensity score, compute IP weights
  fit <- glm(formula = as.formula(miss_formula), family = "binomial", data = data)
  propensity_score <- predict(fit, type = "response")
  ipw_wt <- 1 / propensity_score
  ipw_wtStable <- mean(data$X)*data$X / propensity_score + mean(1 - data$X) * (1 - data$X) / propensity_score
  
  #Based on input var stable, choose which version of the weights
  if(stable == TRUE){
    these_weights <- ipw_wtStable
    this_procedure <- "IPTW_Stable"
  } else {
    these_weights <- ipw_wt
    this_procedure <- "IPTW"
  }
  
  # get IPW estimator and robust SE
  data_cc <- data
  data_cc$weights <- these_weights
  data_cc <- data_cc[data_cc$is.complete == 1, ]
  ipw_fit <- glm(formula = as.formula(outcome_formula), family = "quasibinomial", data = data_cc, weights = data_cc$these_weights)

  est <- summary(ipw_fit)$coef[coefficient_of_interest, 1]
  SE <- sqrt(diag(sandwich::sandwich(ipw_fit)))[coefficient_of_interest]
  
  conditional_est <- data.frame("procedure" = this_procedure, 
                   "estimand" = "canonical_parameter",
                   "est" = est, "SE" = SE)
  marginal_params <- get_marginal_estimates(fit = ipw_fit, data = data_cc, estimator = "iptw")
  marginal_est <- data.frame("procedure" = this_procedure, marginal_params)

  return(list("results" = rbind.data.frame(conditional_est, marginal_est), 
  "coefs" = summary(ipw_fit)$coefficients, "ipw_wt" = ipw_wt, "ipw_wtStable"=ipw_wtStable))
}

# RRZ calibration
# @param data the dataset
# @param outcome_formula a character string specifying the formula for the model, should be phase1 data only
# @param cal_formula the calibration formula
# @param coefficient_of_interest the name of the coefficient of interest
# @param missing_indicator the name of the missingness indicator
# @param fam the outcome regression family
# @return a list, with a data.frame called "results" (contains procedure, estimand, estimate, standard error) and with an object called "fit" that can contain anything else you'd like to save
run_RRZ_lr <- function(data = NULL, outcome_formula = "Y~X+Zs+Zw+Ws_obs", cal_formula="~X+Zs+Zw", coefficient_of_interest = "X", missing_indicator = "is.complete", fam = "binomial") {
 	
	designEstWt <- survey::estWeights(data, formula = as.formula(cal_formula), subset =~ I(data[[missing_indicator]] == 1))
  if (fam == "binomial") {
    fit <- survey::svyglm(as.formula(outcome_formula), design = designEstWt, family = quasibinomial)
  }
  est <- summary(fit)$coef[coefficient_of_interest, 1]
 	SE <- sqrt(diag(sandwich(fit)))[coefficient_of_interest]
  
  conditional_est <- data.frame("procedure" = "RRZ", "estimand" = "canonical_parameter", 
                   "est" = est, "SE" = SE)
  marginal_params <- get_marginal_estimates(fit = fit, data = data[data[[missing_indicator]] == 1, ], estimator = "RRZ")
  marginal_est <- data.frame("procedure" = "RRZ", marginal_params)
  coefs <- data.frame("procedure" = "RRZ", summary(fit)$coefficients)
  names(coefs) <- c("procedure", "Estimate", "SE", "test_statistic", "p_value")
  return(list("results" = rbind.data.frame(conditional_est, marginal_est), 
  "coefs" = cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL))))
}

# raking procedures ------------------------------------------------------------
# @param data the dataset
# @param formula a character string specifying the formula for the model
# @param NimpRaking number of imputations for raking variable
# @param calOption which calibration option to use?
# @param start_from_ipw should we start from the IPW weights?
# @param coefficient_of_interest the coefficient of interest, either a character string (e.g., "X") or its position in a formula (e.g., 2, behind the intercept)
# @param missing_indicator the name of the missingness indicator
# @param fam the outcome regression family
# @return a list with the parameter estimate and standard error (in a data.frame) and the fitted model
run_raking_lr <- function(data = NULL, formula = "Y ~ X", NimpRaking = 50, calOption = 1, 
                          start_from_ipw = FALSE, rake_on_y = FALSE,
                          miss_formula = "is.complete ~ X + Y",
                          coefficient_of_interest = "X", missing_indicator = "is.complete", 
                          fam = "binomial") {
  # get initial weights
  miss_fit <- glm(formula = as.formula(miss_formula), family = "binomial", data = data)
  p_obs <- predict(miss_fit, type = "response")
  ip_weights <- 1 / p_obs
  # Count parameters in outcome model
  if (fam == "binomial") {
    initfit <- glm(formula = as.formula(formula), family = "binomial", data = data)
  }
 	coefnum <- length(coef(initfit))

  #define a dataset with formula covariates and the outcome variable and impute missing data
 	data.mi <- model.frame(as.formula(formula), data = data, na.action = NULL)
 	phase_2_indx <- (data$is.complete == 1)
 	# add on fake phase 2 data with missing values; this is for averaging
 	miss_phase2 <- data.mi[data$is.complete == 1, ]
 	miss_cols <- which(colSums(is.na(data.mi)) > 0)
 	miss_phase2[, miss_cols] <- NA
 	data.mi2 <- rbind(data.mi, miss_phase2)
 	full_phase_1 <- (1:nrow(data.mi2) %in% which(data$is.complete == 0))
 	fake_phase_2 <- ((1:nrow(data.mi2)) %in% (nrow(data.mi) + 1):nrow(data.mi2))
 	data.mi <- data.mi2
  init <- mice::mice(data.mi, maxit = 0)
  pred.matrix <- init$predictorMatrix
  data.imputed <- mice::mice(data.mi, predictorMatrix = pred.matrix, m = NimpRaking, maxit = 20, print = FALSE)
  # get expected value of influence functions by averaging over MI datasets
  infMat_all <- array(data = 0, dim = c(nrow(data), coefnum, NimpRaking))
	for (iter in 1:NimpRaking) {
	  # get imputations for actual phase 1, phase2 observations (i.e., remove real phase 2)
	  imp_init <- mice::complete(data.imputed, iter)
		impData_i <- imp_init[1:nrow(data), ]
		impData_i[phase_2_indx, ] <- imp_init[fake_phase_2, ]
    if (fam == "binomial") {
      mifit <- glm(formula = as.formula(formula), family = binomial, data = impData_i)
    }
    infMat_all[, , iter] <- inf.fun.logit(mifit)
	}
  infMat <- rowMeans(infMat_all, dims = 2)
  
  # choose raking variables and put relevant raking variables into data frame
  # default: generalized raking with all influence functions
  # rakeformula = ~ inf1 + ... + infk, where k = # coef fit by regression
  if (calOption == 1) {
    this_procedure <- "GR_Vanilla"
    rakeformula <- "~ inf1"
    for (i in 1:coefnum) {
      varname <- paste0("inf", i)
      data$inf <- infMat[, i]
      names(data)[names(data) == "inf"] <- varname 
      if (i > 1) {
        rakeformula <- paste0(rakeformula, "+", varname)
      }
    }
  }
  # option 2: calibrate using single influence function for X
  if (calOption == 2) {
    this_procedure <- "GR_X"
    rakeformula <- "~ infX"
    trtColNum <- (1:length(coef(initfit)))[names(coef(initfit)) == coefficient_of_interest]
    data$infX <- infMat[, trtColNum]
  }
	# option 3: also calibrate on y
	if (rake_on_y) {
	  this_procedure <- paste0(this_procedure, "Y")
	  rakeformula <- paste0(rakeformula," + ", "Y")
	}
  data2 <- data
  data2$ipw_wts <- ip_weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~I(data2[[missing_indicator]] == 1), 
    prob = list(NULL, ~I(1 / ipw_wts)), data = data2,
    pps = list(NULL, poisson_sampling(1 / data2$ipw_wts[data2[[missing_indicator]] == 1]))
  )
  infcal <- survey::calibrate(mydesign, formula = as.formula(rakeformula), phase = 2, calfun = "raking")
  
  if (fam == "binomial") {
    rakefit <- survey::svyglm(formula, design = infcal, family = quasibinomial)
  }
	est <- coef(rakefit)[coefficient_of_interest]
	se <- summary(rakefit)$coef[coefficient_of_interest, 2]
  conditional_est <- data.frame("procedure" = this_procedure, "estimand" = "canonical_parameter", "est" = est, "SE" = se)
  # marginal_params <- get_marginal_estimates(fit = rakefit, data = data[data[[missing_indicator]] == 1, ], estimator = this_procedure)
  marginal_params <- marginal_effects_new(fit=rakefit, weights=weights(rakefit), coefficient_of_interest = coefficient_of_interest)
  marginal_est <- data.frame("procedure" = this_procedure, marginal_params)
  coefs <- data.frame("procedure" = this_procedure, summary(rakefit)$coefficients)
  names(coefs) <- c("procedure", "Estimate", "SE", "test_statistic", "p_value")
  return(list("results" = rbind.data.frame(conditional_est, marginal_est),
              "coefs" = cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL))))
}