# Generally helpful functions across input, methods, results
# Update Jan 11 2024: added RF to procedures

# get the latest version number for a given set of functions
# @param all_files the list of files to search through
# @param file the file of interest, e.g., "00_utils"
# @return the latest version of that file, e.g., "20230802"
get_latest_version <- function(all_files = NULL, file = "00_utils") {
  these_files <- all_files[grepl(file, all_files)]
  dates_chr <- gsub(paste0(file, "_"), "", these_files)
  dates <- as.Date(dates_chr, format = "%Y%m%d")
  today <- Sys.Date()
  latest <- dates_chr[which.min(today - dates)]
  return(latest)
}

# Get oracle, working formulas for outcome regression, missing-data model, 
# and propensity score.
# ** THESE SHOULD BE EDITED WHEN NEW SCENARIOS ARE ADDED **
# @param YScenario the outcome scenario
# @param plasmode logical, whether or not this is a plasmode data analysis
# @param rare_outcome logical, whether or not this is a rare-outcome plasmode data analysis
get_outcome_formulas <- function(YScenario = "1", plasmode = FALSE, rare_outcome = FALSE) {
  if (!plasmode) {
    if (grepl("1.", YScenario, fixed = TRUE) | (grepl("1", YScenario) & !grepl(".1", YScenario, fixed = TRUE))) {
      oracle_formula <- "Y ~ X + Zs + Zw + Ws + Ww"
      population_formula <- oracle_formula
      outcome_formula <- "Y ~ X + Zs + Zw + Ws_obs + Ww_obs"
      outcome_formula_confounded <- "Y ~ X + Zs + Zw"
    } else if (grepl("2.", YScenario, fixed = TRUE) | (grepl("2", YScenario) & !grepl(".2", YScenario, fixed = TRUE))) {
      # oracle_formula <- "Y ~ X + Ww + I(Zw < -3) + I(Zw > 3) + Ww * Ws + I(Zs < 3) + Ws + Ws * I(Zs < -3) + Ws * I(Zw > 3)" 
      # Updated 8/28/2024 (Chloe)
      oracle_formula <- "Y ~ X + Zs + Zw + Ws + Ww + Us" 
      population_formula <- "Y ~ X + Zs + Zw + Ws + Ww"
      outcome_formula <- "Y ~ X + Zs + Zw + Ws_obs + Ww_obs"
      outcome_formula_confounded <- "Y ~ X + Zs + Zw"
    } else if (grepl("3.", YScenario, fixed = TRUE)) {
      # y scenario 3: functional form misspecification
      if (grepl(".1", YScenario, fixed = TRUE)) {
        oracle_formula <- "Y ~ X + I(Zs^2) + I(exp(Zw)) + Ws + Ww"
      } else if (grepl(".2", YScenario, fixed = TRUE)) {
        oracle_formula <- "Y ~ X + I(Zs^2) + I(exp(Zw)) + Ws + Ww"
      } else if (grepl(".3", YScenario, fixed = TRUE)) {
        oracle_formula <- "Y ~ X + Zs + Zw + I(Ws^2) + I(exp(Ww))"
      } else {
        oracle_formula <- "Y ~ X + Zs + Zw + I(Ws^2) + I(exp(Ww))"
      }
      population_formula <- "Y ~ X + Zs + Zw + Ws + Ww"
      outcome_formula <- "Y ~ X + Zs + Zw + Ws_obs + Ww_obs"
      outcome_formula_confounded <- "Y ~ X + Zs + Zw"
    } else if (grepl("4.", YScenario, fixed = TRUE) | (grepl("4", YScenario) & !grepl(".4", YScenario, fixed = TRUE))) {
      # y scenario 4: more complex functional form misspecification, interactions
      oracle_formula <- "Y ~ X + Ww + I(Zw < -0.5) + I(Zw > 2) + Ww * Ws + I(Zs < -1) + Ws + Ws * I(Zs < -1) + Ws * I(Zw > 2)"  
      population_formula <- "Y ~ X + Zs + Zw + Ws + Ww"
      outcome_formula <- "Y ~ X + Zs + Zw + Ws_obs + Ww_obs"
      outcome_formula_confounded <- "Y ~ X + Zs + Zw"
    } else {
      # not yet implemented
    }
    # lists out the outcome formula with factors explicit; for TMLE
    outcome_formula_factor <- outcome_formula
  } else {
    oracle_formula <- paste0("Y ~ numEpiType + sex + AgeAtIndex + AgeIndexsq + ",
                             "Charlson_cat + aud_dx_priorYr + anx_dx_priorYr + ",
                             "priorSH + MHIP_prior5yr + ",
                             "IndexTrtPHQ8_score_cat + IndexTrtPHQ_item9_score_cat + ",
                             "Charlson_cat*anx_dx_priorYr + Charlson_cat*AgeAtIndex + ",
                             "sex*AgeAtIndex + sex*priorSH + priorSH*AgeAtIndex + ",
                             "IndexTrtPHQ_item9_score_cat*sex + ",
                             "IndexTrtPHQ_item9_score_cat*priorSH")
    population_formula <- paste0("Y ~ numEpiType + sex + AgeAtIndex_cat + Charlson_cat + ",
                                 "aud_dx_priorYr + anx_dx_priorYr + ",
                                 "priorSH + MHIP_prior5yr + IndexTrtPHQ8_score_cat + ",
                                 "IndexTrtPHQ_item9_score_cat +",
                                 "sex * AgeAtIndex + sex * priorSH + ",
                                 "priorSH * AgeAtIndex + ",
                                 "sex * IndexTrtPHQ_item9_score_cat + ",
                                 "priorSH * IndexTrtPHQ_item9_score_cat")
    outcome_formula_confounded <- paste0("Y ~ numEpiType + sex + AgeAtIndex_cat + Charlson_cat + ",
                                         "aud_dx_priorYr + anx_dx_priorYr + ",
                                         "priorSH + MHIP_prior5yr + ",
                                         "sex * AgeAtIndex + sex * priorSH + ",
                                         "priorSH * AgeAtIndex")
    terms_factor_age <- paste0("AgeAtIndex_cat", 
                               c("25_34", "35_44", "45_54", "55_64", "65plus"))
    if (rare_outcome) {
      terms_factor_charlson <- "Charlson_cat" 
    } else {
      terms_factor_charlson <- paste0("Charlson_cat.", c("L", "Q", "C"))
    }
    terms_factor_phq8 <- paste0("IndexTrtPHQ8_score_cat_obs", 
                                c("6_10", "11_15", "16_20", "20plus"))
    terms_factor_item9 <- paste0("IndexTrtPHQ_item9_score_cat_obs",
                                 c("some_days", "most_days", "nearly_everyday"))
    outcome_formula <- paste0("Y ~ numEpiType + sex + Charlson_cat + ",
                              "aud_dx_priorYr + anx_dx_priorYr + ",
                              "priorSH + MHIP_prior5yr + IndexTrtPHQ8_score_cat_obs + ",
                              "IndexTrtPHQ_item9_score_cat_obs + ",
                              "sex * AgeAtIndex + sex * priorSH + ",
                              "priorSH * AgeAtIndex + ",
                              "sex * IndexTrtPHQ_item9_score_cat_obs")
    outcome_formula_factor <- paste0("Y ~ numEpiType + sex + ",
                                     paste0(terms_factor_charlson, collapse = " + "), " + ",
                                     "aud_dx_priorYr + anx_dx_priorYr + priorSH + MHIP_prior5yr + ",
                                     paste0(terms_factor_phq8, collapse = " + "), " + ",
                                     paste0(terms_factor_item9, collapse = " + "), " + ",
                                     "sex * AgeAtIndex + ", 
                                     "sex * priorSH + ",
                                     "priorSH * AgeAtIndex + ",
                                     paste0(paste0("sex * ", terms_factor_item9, collapse = " + ")))
    if (!rare_outcome) {
      outcome_formula <- paste0(outcome_formula, " + AgeAtIndex_cat",
                                " + priorSH * IndexTrtPHQ_item9_score_cat_obs")
      outcome_formula_factor <- paste0(outcome_formula_factor, " + ",
                                       paste0(terms_factor_age, collapse = " + "), " + ",
                                       paste0(paste0("priorSH * ", terms_factor_item9, collapse = " + ")))
    }
    
  }
  
  return(list("oracle" = oracle_formula, "population" = population_formula,
              "outcome" = outcome_formula, "confounded" = outcome_formula_confounded,
              "outcome_factor" = outcome_formula_factor))
}
# @param mScenario the missing-data scenario
# @param plasmode logical, whether or not this is a plasmode data analysis
get_miss_formula <- function(mScenario = "1", xScenario = "1", plasmode = FALSE, rare_outcome = FALSE) {
  if (!plasmode) {
    miss_formula <- "is.complete ~ X + Zs + Zw"  
    if (grepl("1.1", mScenario, fixed = TRUE)) {
      miss_formula <- paste0(miss_formula, " + Y")
    }
    if (grepl("2.2", mScenario, fixed = TRUE) | grepl("2.4", mScenario, fixed = TRUE)) { # don't know the missing-data mechanism has an interaction
      miss_formula <- paste0(miss_formula, " + Y")
    }  
    if (xScenario == "1.1") { # allow auxiliary variables
      miss_formula <- paste0(miss_formula, " + As + Aw")
    }
    miss_formula_factor <- miss_formula
  } else {
    miss_formula <- paste0("is.complete ~ numEpiType + sex + AgeAtIndex_cat + ",
                           "Charlson_cat + aud_dx_priorYr + anx_dx_priorYr + ",
                           "priorSH + MHIP_prior5yr + ",
                           "AgeAtIndex_cat * sex + priorSH * sex + ",
                           "AgeAtIndex_cat * priorSH")
    terms_factor_age <- paste0("AgeAtIndex_cat", 
                               c("25_34", "35_44", "45_54", "55_64", "65plus"))
    if (rare_outcome) {
      terms_factor_charlson <- "Charlson_cat" 
    } else {
      terms_factor_charlson <- paste0("Charlson_cat.", c("L", "Q", "C"))
    }
    terms_factor_phq8 <- paste0("IndexTrtPHQ8_score_cat_obs", 
                                c("6_10", "11_15", "16_20", "20plus"))
    terms_factor_item9 <- paste0("IndexTrtPHQ_item9_score_cat_obs",
                                 c("some_days", "most_days", "nearly_everyday"))
    miss_formula_factor <- paste0("is.complete ~ numEpiType + sex + ",
                                  paste0(terms_factor_age, collapse = " + "), " + ",
                                  paste0(terms_factor_charlson, collapse = " + "), " + ",
                                  "aud_dx_priorYr + anx_dx_priorYr + priorSH + MHIP_prior5yr + ",
                                  paste0(paste0("sex * ", terms_factor_age, collapse = " + ")), " + ",
                                  "sex * priorSH + ",
                                  paste0(paste0("priorSH * ", terms_factor_age, collapse = " + ")))
  }
  return(list("miss_formula" = miss_formula, "miss_formula_factor" = miss_formula_factor))  
}
#XScenario 1.1 uses same tx formula as 1.  No changes made for it.
# @param XScenario the outcome scenario
# @param plasmode logical, whether or not this is a plasmode data analysis
get_tx_formula <- function(XScenario = "1", plasmode = FALSE, rare_outcome = FALSE) {
  if (!plasmode) {
    if (grepl("1", XScenario)) {
      tx_formula <- "X ~ Zs + Zw + Ws_obs + Ww_obs + As + Aw"
    } else {
      tx_formula <- "X ~ Zs + Zw + Ws_obs + Ww_obs + As + Aw"
    }  
    tx_formula_factor <- tx_formula
  } else {
    tx_formula <- paste0("numEpiType ~ sex + AgeAtIndex_cat + ",
                         "Charlson_cat + aud_dx_priorYr + anx_dx_priorYr + ",
                         "priorSH + MHIP_prior5yr + ",
                         "IndexTrtPHQ8_score_cat_obs + IndexTrtPHQ_item9_score_cat_obs + ",
                         "priorSH * sex + priorSH * AgeAtIndex_cat + ",
                         "sex * AgeAtIndex_cat + ",
                         "IndexTrtPHQ_item9_score_cat_obs * sex + ",
                         "IndexTrtPHQ_item9_score_cat_obs * priorSH")
    terms_factor_age <- paste0("AgeAtIndex_cat", 
                               c("25_34", "35_44", "45_54", "55_64", "65plus"))
    if (rare_outcome) {
      terms_factor_charlson <- "Charlson_cat" 
    } else {
      terms_factor_charlson <- paste0("Charlson_cat.", c("L", "Q", "C"))
    }
    terms_factor_phq8 <- paste0("IndexTrtPHQ8_score_cat_obs", 
                                c("6_10", "11_15", "16_20", "20plus"))
    terms_factor_item9 <- paste0("IndexTrtPHQ_item9_score_cat_obs",
                                 c("some_days", "most_days", "nearly_everyday"))
    tx_formula_factor <- paste0("numEpiType ~ sex + ",
                                paste0(terms_factor_age, collapse = " + "), " + ",
                                paste0(terms_factor_charlson, collapse = " + "), " + ",
                                "aud_dx_priorYr + anx_dx_priorYr + priorSH + MHIP_prior5yr + ",
                                paste0(terms_factor_phq8, collapse = " + "), " + ",
                                paste0(terms_factor_item9, collapse = " + "), " + ",
                                paste0(paste0("sex * ", terms_factor_age, collapse = " + ")), " + ",
                                paste0(paste0("sex * ", terms_factor_item9, collapse = " + ")), " + ",
                                paste0(paste0("priorSH * ", terms_factor_item9, collapse = " + ")))
  }
  return(list("tx_formula" = tx_formula, "tx_formula_factor" = tx_formula_factor))
}

# calculate the expit (equivalent to stats::qlogis)
# @param x a vector
# @return application of the expit function to x
expit<-function(x){
  exp(x)/(1+exp(x))
}
get_ci <- function(est, SE, scale = "identity") {
  if (length(scale) == 1) {
    ci_ll <- switch(as.numeric(scale != "identity") + 1, est - qnorm(0.975) * SE, exp(log(est) - qnorm(0.975) * SE))
    ci_ul <- switch(as.numeric(scale != "identity") + 1, est + qnorm(0.975) * SE, exp(log(est) + qnorm(0.975) * SE))  
  } else {
    ci_ll <- ifelse(scale == "identity", est - qnorm(0.975) * SE, exp(log(est) - qnorm(0.975) * SE))
    ci_ul <- ifelse(scale == "identity", est + qnorm(0.975) * SE, exp(log(est) + qnorm(0.975) * SE))  
  }
  # needs to be able to return a vector or a matrix
  return(cbind(ci_ll, ci_ul))
}
# coverage of a 95\% confidence interval
# @param mu the true parameter value
# @param est the estimate
# @param SE the standard error
# @param scale what scale should the CI be computed on? 
# @return TRUE if mu is contained in the CI
cover <- function(mu,est,SE,scale = "identity"){
  # Edited by Chloe April 26, 2024: believe issue with line below finding ci_ll and ci_ul
  cis <- get_ci(est = est, SE = SE, scale = scale)
  ci_ll <- cis[, 1]
  ci_ul <- cis[, 2]
  ## return indicator mu in the interval (est-1.96SE, est+1.96SE)
  return((ci_ll <= mu) * (ci_ul >= mu))
  # return((mu > est-qnorm(.975)*SE)* (mu<est+qnorm(0.975)*SE))
}
power <- function(mu, est, SE, scale = "identity") {
  cis <- get_ci(est = est, SE = SE, scale = scale)
  ci_ll <- cis[, 1]
  ci_ul <- cis[, 2]
  greater_0 <- ci_ll > mu
  less_0 <- ci_ul < mu
  return(greater_0 + less_0)
}

#### Logistic regression influence function  - for Raking/survey calibration
### Tong's binary influence function for binary data
inf.fun.logit <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat) / nrow(dm)
  infl
}
#### Boe el Sandwich paper influence function calculation in Figure1 
### Fit must be a svyglm object to use the fit$naiv.cov
### Can create design object for svyglm for srs by this statement
### mydata.svy<-svydesign(id=~1, data = mydata)

inf.fun.boe<-function(fit){
  if (inherits(fit, "svyglm")) {
    A <- -solve(fit$naive.cov) / nrow(fit$data)
  } else {
    A <- -solve(summary(fit)$cov.unscaled) / nrow(fit$data)
  }
  U <- as.matrix(estfun(fit))
  U %*% solve(A) / nrow(fit$data)
}


# NEW November 8 (added by Chloe): getting marginal effects using marginaleffects package.
marginal_effects_new <- function(fit, weights, coefficient_of_interest = "X") {
  if (inherits(fit, "svyglm")) {
    marginaleffects_vcov <- TRUE
  } else {
    marginaleffects_vcov <- sandwich::sandwich
  }
  # For simplicity and following others' advice, I'm using the Hajek estimator for IPW, but an option for HT can be added here.
  rd_comparison <- avg_comparisons(fit, comparison="differenceavg", variables = coefficient_of_interest, wts=weights,
                                   vcov = marginaleffects_vcov)
  rr_comparison <- avg_comparisons(fit, comparison="lnratioavg", variables = coefficient_of_interest, wts=weights,
                                   vcov = marginaleffects_vcov)
  or_comparison <- avg_comparisons(fit, comparison="lnoravg", variables = coefficient_of_interest, wts=weights,
                                   vcov = marginaleffects_vcov)
  
  results <- rbind(c(rd_comparison[1,"estimate"], rd_comparison[1,"std.error"]),
                   c(rr_comparison[1,"estimate"], rr_comparison[1,"std.error"]),
                   c(or_comparison[1,"estimate"], or_comparison[1,"std.error"]))
  estimand <- c("RD","logRR","logOR")
  results <- cbind.data.frame(estimand,results)
  colnames(results) <- c("estimand","est","SE")
  return(results)
}

# get marginal parameter estimates from conditional approach
# @param fit the conditional model fit
# @param data the dataset
# @param estimator the name of the estimator, e.g., "oracle"
### Comments/edits added by Chloe to this function September 15, 2023.
get_marginal_estimates <- function(fit, data, estimator = "oracle") {
  data_0 <- data
  data_0$X <- 0
  data_1 <- data
  data_1$X <- 1
  if (inherits(fit, "glm")) {
    # preds_0 <- predict(fit, newdata = data_0, type = "response")
    preds_0 <- predict(fit, newdata = data_0, type = "response",na.action=na.omit) 
    # I *think* this is needed to be safe...the default output is N long for some estimators, and length of complete cases for others.
    # Since the weights object is always of the length of complete cases, I believe this could have led to incorrect estimates for IPW previously.
    # preds_1 <- predict(fit, newdata = data_1, type = "response")
    preds_1 <- predict(fit, newdata = data_1, type = "response",na.action=na.omit)
    # Added for general use, requires modification for IPW.
    denom <- length(preds_1)
  } else {
    stop("Non-glm regression functions currently not supported.")
  }
  
  if (grepl("oracle", estimator, ignore.case = TRUE) | grepl("population", estimator, ignore.case = TRUE) | grepl("mi", estimator, ignore.case = TRUE) | grepl("noW", estimator, ignore.case = TRUE)) {
    weights <- rep(1, nrow(data))
  } else if (grepl("cc", estimator, ignore.case = TRUE) & !grepl("oracle", estimator, ignore.case = TRUE)) {
    weights <- rep(1, denom)
    # Added by Chloe 9/15/2023
  } else if (grepl("rak", estimator, ignore.case = TRUE) | grepl("gr", estimator, ignore.case = TRUE)) { 
    # weights <- fit$weights
    weights <- weights(fit) # Could also use the weights from the survey object, i.e. weights(infcal), the weights used
    # to fit the regression, and divide below by n, i.e. replicating what is being done with the IPW estimator below.
    # However, weights(fit) are re-scaled to the complete sample size.
  } else { # Intended for IPW estimators (I've investigated standard, not stabilized.)
    # weights <- fit$weights
    weights <- weights(fit)
    denom <- n # Need to divide by total original sample size for IPW; also true for raking using fitting weights,
    # but re-standardized weights from weights(fit) can be divided by the sample size with complete data.
  } # Note from Chloe 9/18/2023: I have not looked into by-hand calculation for IPTW, but presumably it would be the same
  # as for IPW in last "else" statement above.
  # mu_1 <- mean(weights * preds_1)
  mu_1 <- sum(weights * preds_1)/denom
  # mu_0 <- mean(weights * preds_0)
  mu_0 <- sum(weights * preds_0)/denom
  rd <- mu_1 - mu_0
  rr <- mu_1 / mu_0
  or <- (mu_1 / (1 - mu_1)) / (mu_0 / (1 - mu_0))
  return(data.frame("estimand" = c("RD", "RR", "OR"),
                    "est" = c(rd, rr, or),
                    "SE" = c(NA, NA, NA)))
}

# fix factor variables, for TMLE -----------------------------------------------
# @param data the dataset, possibly with factors
fix_factor_variables <- function(data) {
  factor_vars <- sapply(data, is.factor)
  if (!any(factor_vars)) {
    return(data)
  } else {
    data2 <- data
    for (i in 1:sum(factor_vars)) {
      tmp <- data[, factor_vars][, i]
      tmp2 <- plyr::revalue(tmp, replace = c("18-24" = "18_24", "25-34" = "25_34",
                                             "35-44" = "35_44", "45-54" = "45_54",
                                             "55-64" = "55_64", "65+" = "65plus",
                                             "3+" = "3plus",
                                             "<5" = "less5", "6-10" = "6_10",
                                             "11-15" = "11_15", "16-20" = "16_20",
                                             "20+" = "20plus",
                                             "Some of the days" = "some_days",
                                             "Most of the days" = "most_days",
                                             "Nearly everyday" = "nearly_everyday"),
                            warn_missing = FALSE)
      data2[, factor_vars][, i] <- tmp2
    }
    return(data2) 
  }
}
# @param data the dataset, possibly with factors
# @param outcome_formula the outcome regression formula
# @return a list of datasets with factors turned into dummy variables
make_dummy_variables <- function(data, outcome_formula, tx_formula, W, V, w_sub,
                                 phase1_covars_no_AY, phase2_covars) {
  factor_vars <- sapply(data, is.factor)
  if (!any(factor_vars)) {
    W_msm <- W
    V_msm <- V
    W_sub_msm <- w_sub
  } else {
    q_model_mat <- model.matrix(as.formula(outcome_formula), 
                                model.frame(~ ., data = data, na.action = "na.pass"))[, -1]
    g_model_mat <- model.matrix(as.formula(tx_formula), 
                                model.frame(~ ., data = data, na.action = "na.pass"))[, -1]
    data2_y <- data.frame(q_model_mat)
    data2_a <- data.frame(g_model_mat)
    data2 <- cbind(data.frame(q_model_mat), data.frame(g_model_mat)[!(colnames(g_model_mat) %in% colnames(q_model_mat))])
    phase1_covars_no_AY_msm_W <- sapply(names(data2), function(x) {
      any(sapply(phase1_covars_no_AY, function(z) {
        grepl(z, x) & !any(sapply(phase2_covars, function(w) grepl(w, x)))
      }))
    })
    phase1_covars_no_AY_msm_V <- sapply(names(data2_y), function(x) {
      any(sapply(phase1_covars_no_AY, function(z) {
        grepl(z, x) & !any(sapply(phase2_covars, function(w) grepl(w, x)))
      }))
    })
    phase2_covars_msm <- sapply(names(data2_y), function(x) {
      any(sapply(phase2_covars, function(z) grepl(z, x)))
    })
    W_msm <- data2[, phase1_covars_no_AY_msm_W, drop = FALSE]
    V_msm <- data2_y[, phase1_covars_no_AY_msm_V | phase2_covars_msm, drop = FALSE]
    phase_2_msm <- data2_y[, phase2_covars_msm, drop = FALSE]
    W_sub_msm <- phase_2_msm[complete.cases(phase_2_msm), , drop = FALSE]
  }
  return(list("W" = W_msm, "V" = V_msm, "W_sub" = W_sub_msm))
}

# for nice tables --------------------------------------------------------------
get_nice_procedure <- function(proc) {
  dplyr::case_when(
    proc == "cc" ~ "cc",
    proc == "cc_noW" ~ "cc_noW",
    proc == "cc_oracle" ~ "cc_oracle",
    proc == "cc_population" ~ "cc_population",
    proc == "gr" ~ "GR_Vanilla",
    proc == "rak" ~ "GR_Vanilla",
    proc == "ipw" ~ "IPW",
    proc == "mice" ~ "MICE",
    proc == "RF" | proc == "rf" ~ "RF",
    proc == "xgb" ~ "XGB",
    proc == "tmle_m" ~ "tmle_m",
    proc == "tmle_mto" ~ "tmle_mto",
    proc == "ipcw-tmle_m" ~ "ipcw-tmle_m",
    proc == "ipcw-tmle_mto" ~ "ipcw-tmle_mto",
    proc == "ipcw-a-tmle_m" ~ "ipcw-a-tmle_m",
    proc == "ipcw-a-tmle_mto" ~ "ipcw-a-tmle_mto",
    proc == "r-ipcw-tmle_m" ~ "r-ipcw-tmle_m",
    proc == "r-ipcw-tmle_mto" ~ "r-ipcw-tmle_mto",
    proc == "r-ipcw-a-tmle_m" ~ "r-ipcw-a-tmle_m",
    proc == "r-ipcw-a-tmle_mto" ~ "r-ipcw-a-tmle_mto",
    # proc == "tmle_m" ~ "TMLE-M",
    # proc == "tmle_mto" ~ "TMLE-MTO",
    # proc == "ipcw-a-tmle_m" ~ "IPCW-a-TMLE-M",
    # proc == "ipcw-a-tmle_mto" ~ "IPCW-a-TMLE-MTO",
    # proc == "r-ipcw-tmle_m" ~ "r-IPCW-TMLE-M",
    # proc == "r-ipcw-tmle_mto" ~ "r-IPCW-TMLE-MTO",
    # proc == "r-ipcw-a-tmle_m" ~ "r-IPCW-a-TMLE-M",
    # proc == "r-ipcw-a-tmle_mto" ~ "r-IPCW-a-TMLE-MTO",
    proc != "cc" & proc != "cc_noW" & proc != "cc_oracle" & proc != "cc_population" &
      proc != "gr" & proc != "rak" & proc != "ipw" & proc != "mice" & proc != "RF" & proc != "xgb" &
      proc != "tmle_m" & proc != "tmle_mto" & proc != "ipcw-tmle_m" & proc != "ipcw-tmle_mto" & proc != "ipcw-a-tmle_m" & proc != "ipcw-a-tmle_mto" &
      proc != "r-ipcw-tmle_m" & proc != "r-ipcw-tmle_mto" & proc != "r-ipcw-a-tmle_m" & proc != "r-ipcw-a-tmle_mto" ~ "Unknown procedure"
  )
}

get_label <- function(xscenario = "1", yscenario = "1",
                      mscenario = "1", estimand = "cOR") {
  paste0("x", xscenario, "_m", mscenario, "_y", yscenario,
         "_", estimand)
}

# get a nice table caption for latex tables
# @param xscenario the X scenario
# @param yscenario the Y scenario
# @param mscenario the M scenario
# @param estimand the estimand
# @param n the sample size
# @param nreps the number of replications
# @param tmle_nreps the number of replications for TMLE
# @param rd_mult the multiplier for risk difference results
# @param cont_mce, bin_mce the maximum observed Monte-Carlo error for continuous, binary summaries, respectively
# @param oracle are we comparing to the oracle truth or the population truth?
# @param plasmode is this a plasmode simulation?
# @param true_value the true value that estimators are aiming to estimate
# @param bold use bold (TRUE) or star (FALSE) to emphasize mismatched/matched estimators?
# @param oracle_census_equal are the oracle and census equal? If so, suppress printing the emphasis text
get_caption <- function(xscenario = "1", yscenario = "1",
                        mscenario = "1", estimand = "cOR",
                        n = 1000, nreps = 2500, tmle_nreps = 1000, rd_mult = 10,
                        cont_mce = 0, bin_mce = 0, oracle = TRUE,
                        plasmode = FALSE, true_value = 0, bold = FALSE,
                        oracle_census_equal = FALSE) {
  nice_estimand <- case_when(
    estimand == "cOR" ~ "conditional odds ratio (cOR)",
    estimand == "mOR" ~ "marginal odds ratio (mOR)",
    estimand == "mRD" ~ "marginal risk difference (mRD)",
    estimand == "mRR" ~ "marginal relative risk (mRR)"
  )
  estimand_descr <- ifelse(estimand == "mRD", paste0(" Note: Bias, ESE, ASE, MAD, and RMSE scaled by a factor of ", rd_mult, " to facilitate comparisons across estimands."), "")
  nreps_descr <- paste0(nreps, " simulation replications", 
                        ifelse(tmle_nreps != nreps & estimand != "cOR", 
                               paste0(" (only ", tmle_nreps, " replications for TMLE)."), ""))
  abbreviations <- paste0(" ESE = empirical standard error, ASE = asymptotic standard error, MAD = mean absolute deviation, ",
                          "RMSE = root mean squared error, rRMSE = robust RMSE (using median bias and MAD), ",
                          "Oracle coverage = coverage of a confidence interval based on the ESE, ",
                          "Nominal coverage = coverage of a confidence interval based on the ASE.")
  if (plasmode) {
    # nice_y <- gsub("_", "\\_", gsub("tree_", "", gsub("glm_", "", as.character(yscenario))), fixed = TRUE)
    if (grepl("SH", yscenario)) {
      nice_y_txt <- " self-harm or hospitalization"
    } else {
      nice_y_txt <- " self-harm"
    }
    if (grepl("365", yscenario)) {
      nice_time <- "1-year"
    } else {
      nice_time <- "5-year"
    }
    nice_y <- paste0(nice_time, nice_y_txt)
    est_txt <- ifelse(grepl("glm", yscenario), "glms", "trees")
    if (is.na(true_value)) {
      truth_descr <- "The true value of the estimand is undefined, because the oracle cOR does not exist for tree-based data-generating models."   
    } else {
      truth_descr <- paste0(" The value of the estimand is ", round(as.numeric(true_value), 4), ".")
    }
    if (oracle_census_equal) {
      bold_descr <- NULL
    } else {
      bold_descr <- paste0(" Estimators that are mismatched with the estimand (i.e., are estimating a different parameter) are emphasized using ", 
                           switch(as.numeric(bold) + 1, "a star", "bold"), ".")  
    }
    n2 <- n
    if (n == "full") {
      n2 <- "50,337"
    }
    paste0("\\textbf{Plasmode data simulation: ", nice_y, ", regression functions are ", est_txt,
           ", ", ifelse(oracle, "oracle ", "census "), nice_estimand, "}.", 
           " Relative performance of estimators with sample size n = ", n2, " and ", nreps_descr, ".", estimand_descr,
           truth_descr, abbreviations, bold_descr)
    # paste0("\\textbf{Scenario: plasmode, Y = ", nice_y, ", regression functions are ", est_txt,
    #        ", ", nice_estimand, ifelse(oracle, " (oracle truth)", " (census truth)"), "}.", 
    #        " Relative performance of estimators with sample size n = ", n2, " and ", nreps_descr, estimand_descr,
    #        truth_descr, bold_descr)
  } else {
    nice_x_descr <- switch(as.numeric(xscenario == "1") + 1,
                           "unspecified",
                           "simple treatment model")
    if (yscenario == "1" | yscenario == "1.1" | yscenario == "1.15" | yscenario == "1.16" | yscenario == "1.17") {
      if (grepl(".", yscenario, fixed = TRUE)) {
        short_y_descr <- "simple outcome"  
      } else {
        short_y_descr <- "simple outcome (no treatment effect)"
      }
      full_y_descr <- ""
    } else if (grepl("3.", yscenario, fixed = TRUE)) {
      prefix <- " The semi-complex "
      short_y_descr <- "semi-complex outcome"
      if (grepl(".1", yscenario, fixed = TRUE) | grepl(".2", yscenario, fixed = TRUE)) {
        middle <- "outcome is a function of exponentiated and squared terms on Z."
      } else {
        middle <- "outcome is a function of exponentiated and squared terms on W."
      }
      full_y_descr <- paste0(prefix, middle)
    } else if (grepl("4", yscenario, fixed = TRUE) & !grepl(".4", yscenario, fixed = TRUE)) {
      short_y_descr <- "complex outcome"
      prefix <- "complex outcome model with dependence on binary or categorical Z and interactions between the Ws "
      if (grepl(".1", yscenario, fixed = TRUE)) {
        middle <- "with conditional OR of treatment = log(1.5)"
      } else {
        middle <- "and between the strong W and both Zs, with no treatment effect"
      }
      full_y_descr <- ""
    } else if (grepl("2.1", yscenario, fixed = TRUE)) {
      short_y_descr <- "simple outcome (unobserved covariate)"
      full_y_descr <- ""
    } else {
      # not yet implemented
    }
    
    outcome_proportion <- ifelse(grepl(".17", yscenario, fixed = TRUE), "5\\%", "12\\%")
    missing_proportion <- ifelse(mscenario %in% c("3.1", "2.3", "2.4", "2.7", "2.8"), "80\\%", "40\\%")
    if (mscenario %in% c("1", "1.1", "3", "3.1")) {
      if (grepl(".1", mscenario)) {
        nice_m_descr <- "simple MAR"
      } else {
        nice_m_descr <- "simple MAR (no dependence on Y)"
      }
    } else if (mscenario %in% c("2.2", "2.4")) {
      nice_m_descr <- "complex MAR"
    } else if (mscenario %in% c("2.1", "2.3")) {
      nice_m_descr <- "complex MAR (no dependence on Y)"  
    } else if (mscenario %in% c("2.5", "2.7")){
      nice_m_descr <- "MNAR-unobserved"
    } else {
      nice_m_descr <- "MNAR-value"
    }
    # nice_scenario_descr <- paste0(nice_x_descr, ", ", nice_m_descr, ", and ",
    #                               nice_y_descr)
    nice_scenario_descr <- paste0(short_y_descr, " and ", nice_m_descr)
    nice_cont_mce <- ifelse(cont_mce < 1e-4, "< 0.0001", sprintf("%.3f", cont_mce))
    nice_bin_mce <- ifelse(bin_mce < 1e-4, "< 0.0001", sprintf("%.3f", bin_mce))
    mce_descr <- paste0(" Maximum observed Monte-Carlo error over the ", nreps_descr, " was ", nice_cont_mce,
                        " for all summaries besides coverage and ", nice_bin_mce,
                        " for coverage.")
    truth_descr <- paste0(" The value of the estimand is ", round(as.numeric(true_value), 4), ".")
    if (oracle_census_equal) {
      bold_descr <- NULL
    } else {
      bold_descr <- paste0(" Estimators that are mismatched with the estimand (i.e., are estimating a different parameter) are emphasized using ", 
                           switch(as.numeric(bold) + 1, "a star", "bold"), ".")  
    }
    if (mscenario %in% c("2.5", "2.6", "2.7", "2.8")) {
      m_descr <- "MNAR"
    } else {
      m_descr <- "MAR"
    }
    estimand_descr <- paste0(ifelse(oracle, "oracle ", "census "), nice_estimand)
    paste0("\\textbf{Synthetic data ", m_descr, " simulation: ", 
           estimand_descr, ", ", outcome_proportion, " outcome proportion, ", missing_proportion, " missing proportion}.",
           " Comparing estimators under the \\textbf{", nice_scenario_descr, "} scenario.",
           truth_descr, " The sample size is n = ", n, ".", mce_descr,
           # " Comparing estimators of the ", ifelse(oracle, "oracle ", "census "), estimand, " estimand under a ", nice_scenario_descr,
           abbreviations, bold_descr, full_y_descr)
    # paste0("\\textbf{Scenario: Y = ", yscenario, ", M = ", mscenario, ", X = ", xscenario,
    #        ", ", nice_estimand, ifelse(oracle, " (oracle truth)", " (census truth)"), "}.", 
    #        " Relative performance of estimators under a ",
    #        nice_scenario_descr, " with sample size n = ", n, " and ", nreps_descr, estimand_descr,
    #        mce_descr, truth_descr, bold_descr)
  }
}

# put emphasis on estimators that are matched (or mismatched) with a particular estimand
# @param ests the vector of estimator names
# @param oracle if TRUE, the oracle estimand; else the census estimand
# @param matched whether we should emphasize matches (TRUE) or mismatches with the estimand
# @param oracle_equals_census whether the oracle and census estimands are the same
# @param bold if TRUE, bold the text; if FALSE, add a star next to the text
# @param cOR if TRUE, this is the cOR; otherwise, a marginal parameter
emphasize_estimator_estimand_pairs <- function(ests = NULL, oracle = TRUE, matched = TRUE,
                                               oracle_equals_census = FALSE,
                                               bold = FALSE, cOR = FALSE) {
  ests_chr <- as.character(ests)
  if (bold) {
    prefix <- "\\textbf{"
    suffix <- "}"
  } else {
    prefix <- ""
    suffix <- "${}^*$"
  }
  is_benchmark <- grepl("benchmark", ests_chr, ignore.case = TRUE)
  if (matched) {
    # emphasize when they are matched
    if (oracle_equals_census) {
      # everything gets emphasized
      ret <- paste0(prefix, ests_chr, suffix)
    } else {
      ret <- ests_chr
      is_tmle <- grepl("TMLE", ests_chr, ignore.case = TRUE) 
      if (oracle & !cOR) {
        ret[is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[is_tmle & !is_benchmark], suffix)
      } else if (cOR) {
        # everything is targeted towards the census cOR
        ret <- paste0(prefix, ests_chr, suffix)
      } else {
        ret[!is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[!is_tmle & !is_benchmark], suffix)
      }
    }
  } else {
    # then we emphasize if they are mismatched
    if (oracle_equals_census) {
      ret <- ests_chr
    } else {
      ret <- ests_chr
      is_tmle <- grepl("TMLE", ests_chr, ignore.case = TRUE)
      if (oracle & !cOR) {
        ret[!is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[!is_tmle & !is_benchmark], suffix)
      } else if (oracle & cOR) {
        # everything is mismatched for oracle cOR if oracle not equal census
        # ret[is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[is_tmle & !is_benchmark], suffix)
        ret[!is_benchmark] <- paste0(prefix, ests_chr[!is_benchmark], suffix)
      } else if (cOR) {
        # everything is matched for census cOR
        # ret[is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[is_tmle & !is_benchmark], suffix)
      } else {
        # tmle is mismatched for census non-cOR
        ret[is_tmle & !is_benchmark] <- paste0(prefix, ests_chr[is_tmle & !is_benchmark], suffix)
      }
    }
  }
  return(ret)
}

# summarize across the variables of interest for reporting ---------------------
get_summaries <- function(output, num_finished_and_reasonable) {
  summary_tib <- output %>% 
    group_by(yscenario, mscenario, xscenario, nice_procedure, estimand, seed,
             truth_oracle, truth_pop) %>% 
    summarize(mean_bias_pop_truth = mean(bias_pop_init, na.rm = TRUE),
              mean_bias_oracle_truth = mean(bias_oracle_init, na.rm = TRUE),
              median_bias_pop_truth = median(bias_pop_init, na.rm = TRUE),
              median_bias_oracle_truth = median(bias_oracle_init, na.rm = TRUE),
              ESE = sd(est, na.rm = TRUE), ASE = mean(SE, na.rm = TRUE),
              MAD = mad(est, na.rm = TRUE),
              median_SE = median(SE, na.rm = TRUE),
              nominal_coverage_pop = mean(nominal_cover_pop_init, na.rm = TRUE),
              nominal_coverage_oracle = mean(nominal_cover_oracle_init, na.rm = TRUE),
              oracle_coverage_pop = mean(oracle_cover_pop_init, na.rm = TRUE),
              oracle_coverage_oracle = mean(oracle_cover_oracle_init, na.rm = TRUE),
              mad_coverage_pop = mean(mad_cover_pop_init, na.rm = TRUE),
              mad_coverage_oracle = mean(mad_cover_oracle_init, na.rm = TRUE),
              power = mean(power_init, na.rm = TRUE),
              oracle_power = mean(oracle_power_init, na.rm = TRUE),
              mad_power = mean(mad_power_init, na.rm = TRUE),
              .groups = "drop") %>% 
    mutate(rmse_oracle_pop = sqrt(mean_bias_pop_truth ^ 2 + ESE ^ 2),
           rmse_oracle_oracle = sqrt(mean_bias_oracle_truth ^ 2 + ESE ^ 2),
           rmse_pop = sqrt(mean_bias_pop_truth ^ 2 + ASE ^ 2),
           rmse_oracle = sqrt(mean_bias_oracle_truth ^ 2 + ASE ^ 2),
           rmse_mad_med_bias_pop = sqrt(median_bias_pop_truth ^ 2 + MAD ^ 2),
           rmse_mad_med_bias_oracle = sqrt(median_bias_oracle_truth ^ 2 + MAD ^ 2),
           rmse_med_SE_med_bias_pop = sqrt(median_bias_pop_truth ^ 2 + median_SE ^ 2),
           rmse_med_SE_med_bias_oracle = sqrt(median_bias_oracle_truth ^ 2 + median_SE ^ 2)) %>% 
    left_join(num_finished_and_reasonable, by = c("yscenario", "mscenario", 
                                                  "xscenario", "nice_procedure", 
                                                  "estimand", "seed"))
  return(summary_tib)
}

# create a tibble with the relevant summary data, depending on if robust statistics (i.e., medians) are requested ----
get_summary_tib <- function(output, rd_mult = 1, reasonable_threshold = log(10),
                            threshold_scale = "log", include_medse = FALSE) {
  if (threshold_scale == "log") {
    threshold_txt <- paste0("log(", exp(reasonable_threshold), ")")
  } else {
    threshold_txt <- as.character(reasonable_threshold)
  }
  reasonable_nm <- eval(parse(text = paste0("c(`Prop. < ", threshold_txt, "` = 'prop_reasonable')")))
  summary_tib_init <- output %>%
    filter(seed == 1) %>% 
    select(-seed) %>% 
    mutate(across(contains(c("bias", "rmse", "ESE", "ASE", "MAD"), ignore.case = FALSE),
                  .fns = ~ .x * rd_mult))
  if (!include_medse) {
    summary_tib_init <- summary_tib_init %>% 
      select(-contains("med_SE_med_bias"), -contains("median_SE"))  
  }
  summary_tib <- summary_tib_init %>% 
    mutate(across(where(is.double), 
                  .fns = ~ round(.x, digits = 3))) %>% 
    rename(`Y scenario` = yscenario, `M scenario` = mscenario,
           `X scenario` = xscenario, `Estimator` = nice_procedure,
           `Estimand` = estimand, 
           `Mean bias (pop. truth)` = mean_bias_pop_truth, 
           `Mean bias (oracle truth)` = mean_bias_oracle_truth,
           `Median bias (pop. truth)` = median_bias_pop_truth,
           `Median bias (oracle truth)` = median_bias_oracle_truth, 
           `Nominal coverage (pop. truth)` = nominal_coverage_pop,
           `Nominal coverage (oracle truth)` = nominal_coverage_oracle,
           `ESE RMSE (pop. truth)` = rmse_oracle_pop,
           `ESE RMSE (oracle truth)` = rmse_oracle_oracle,
           `RMSE (pop. truth)` = rmse_pop,
           `RMSE (oracle truth)` = rmse_oracle,
           `rRMSE (pop. truth)` = rmse_mad_med_bias_pop,
           `rRMSE (oracle truth)` = rmse_mad_med_bias_oracle,
           `Oracle coverage (pop. truth)` = oracle_coverage_pop, 
           `Oracle coverage (oracle truth)` = oracle_coverage_oracle,
           `MAD coverage (pop. truth)` = mad_coverage_pop,
           `MAD coverage (oracle truth)` = mad_coverage_oracle,
           `Power` = power,
           `ESE power` = oracle_power,
           `MAD power` = mad_power,
           `Prop. completed` = prop_complete) %>%
    rename(all_of(reasonable_nm)) %>% 
    mutate(old_estimand = Estimand,
           Estimand = gsub("log", "", Estimand)) %>%  # Keeping log in name to keep track.
    mutate(across(contains(c("bias", "rmse", "ESE", "ASE", "MAD")),
                  .fns = ~ ifelse(abs(as.numeric(.x)) > reasonable_threshold, paste0("abs > ", threshold_txt), .x)))
  return(summary_tib)
}

# create a knitr kable for printing as part of report --------------------------
create_report_kable <- function(output, oracle_census_equal = TRUE, estimand = "cOR", oracle = TRUE,
                                matched = FALSE, this_scenario = list("xscenario" = 1, "yscenario" = 1, "mscenario" = 1),
                                ns = 10000, nreps = 2500, tmle_nreps = 2500, 
                                max_binary_mce = 0.001, max_continuous_mce = 0.001, 
                                true_val = 0, col_width = "2cm", robust = TRUE,
                                zs = FALSE, ws = FALSE, plasmode = FALSE) {
  z_w_txt <- paste0(ifelse(zs, "_zs", ""), ifelse(ws, "_ws", ""))
  col_width_units <- ifelse(grepl("cm", col_width), "cm", "in")
  double_col_width <- paste0(2 * as.numeric(gsub("in", "", gsub("cm", "", col_width))), col_width_units)
  if (robust) {
    stats_of_interest <- output %>%
      select(Estimator, contains("median"), contains("mad"), contains("rr"), contains("prop")) %>% 
      select(Estimator, contains("bias"), MAD, contains("rmse"), contains("coverage"), contains("power"), everything())
  } else {
    stats_of_interest <- output %>% 
      select(Estimator, contains("mean"), contains("median"), ESE, ASE, MAD,
             contains("rmse") & !contains("rrmse") & !contains("ESE"), contains("rr"), contains("coverage") & !contains("mad"), 
             contains("power") & !contains("mad") & !contains("ESE"), contains("prop")) %>% 
      select(Estimator, contains("bias"), ESE, ASE, MAD,
             contains("rmse") & !contains("ESE"), contains("rmse") & contains("ESE"), 
             contains("coverage") & contains("oracle"), contains("coverage") & !contains("oracle"), 
             contains("power") & contains("ESE"), contains("power") & !contains("ESE"), 
             everything())
  }
  if (oracle) {
    oracle_census_txt <- paste0("_oracle", z_w_txt)
    this_kable_init <- stats_of_interest %>% 
      select(!contains("pop. truth")) %>% 
      rename_with(~str_remove(., fixed(" (oracle truth)"))) %>% 
      filter(`Estimator` != "Population model") %>%
      mutate(`Estimator` = forcats::fct_recode(Estimator, "Benchmark model" = "Oracle model"),
             `Estimator` = emphasize_estimator_estimand_pairs(ests = `Estimator`,
                                                              oracle = oracle, 
                                                              matched = matched,
                                                              oracle_equals_census = oracle_census_equal,
                                                              cOR = (estimand == "cOR")))
  } else {
    oracle_census_txt <- paste0("_census", z_w_txt)
    this_kable_init <- stats_of_interest %>% 
      select(!contains("oracle truth")) %>% 
      rename_with(~str_remove(., fixed(" (pop. truth)"))) %>% 
      filter(`Estimator` != "Oracle model") %>%
      mutate(`Estimator` = forcats::fct_recode(Estimator, "Benchmark model" = "Population model"),
             `Estimator` = emphasize_estimator_estimand_pairs(ests = `Estimator`,
                                                              oracle = oracle, 
                                                              matched = matched,
                                                              oracle_equals_census = oracle_census_equal,
                                                              cOR = (estimand == "cOR")))
  }
  this_kable <- this_kable_init %>% 
    knitr::kable(format = "latex",  escape = FALSE,
                 label = paste0(get_label(xscenario = this_scenario$xscenario,
                                          yscenario = this_scenario$yscenario,
                                          mscenario = this_scenario$mscenario,
                                          estimand = estimand),
                                oracle_census_txt),
                 caption = get_caption(xscenario = this_scenario$xscenario,
                                       yscenario = this_scenario$yscenario,
                                       mscenario = this_scenario$mscenario,
                                       estimand = as.character(estimand),
                                       plasmode = plasmode,
                                       n = ns[1], nreps = nreps, tmle_nreps = tmle_nreps,
                                       cont_mce = max_continuous_mce, bin_mce = max_binary_mce,
                                       oracle = oracle,
                                       true_value = true_val,
                                       oracle_census_equal = oracle_census_equal)) %>% 
    kableExtra::column_spec(2:13, width = col_width) %>% 
    kableExtra::column_spec(1, width = double_col_width) %>% 
    kableExtra::kable_styling(full_width = FALSE, font_size = 9)
  return(this_kable)
}


# get true values for simulation -----------------------------------------------
get_preds <- function(fit, df) {
  if (isTRUE(inherits(fit, "ranger"))) {
    pred_obj <- predict(fit, data = df, type = "response")$predictions
    nms <- colnames(pred_obj)
    if (!any(grepl("1", nms))) {
      colnames(pred_obj) <- as.character(as.numeric(as.logical(nms)))
    }
    preds <- pred_obj[, "1"]
  } else {
    preds <- predict(fit, newdata = df, type = "response")
  }
  return(preds)
}
get_oracle_cOR <- function(fit, coef_name) {
  if (isTRUE(inherits(fit, "glm"))) {
    oracle_cor <- coefficients(fit)[coef_name]
  } else {
    oracle_cor <- NA
  }
  return(oracle_cor)
}
get_true_values <- function(all_xym = NULL, iter = 1, sample_size = 1e6, nreps = 10,
                            lowcor = 0.2, midcor = 0.4, highcor = 0.7, gencor = 0.2,
                            y_only = FALSE, plasmode = FALSE, 
                            census_formula = "Y ~ X + Zs + Zw + Ws + Ww",
                            plasmode_folder = paste0(getwd(), "/Data/plasmode data sets/")) {
  this_scenario <- all_xym[iter, ]
  xscenario <- this_scenario$xscenario
  yscenario <- this_scenario$yscenario
  mscenario <- this_scenario$mscenario
  outcome_formulas <- get_outcome_formulas(YScenario = yscenario, plasmode = plasmode)
  census_formula <- outcome_formulas$population
  outcome_name <- strsplit(census_formula, "~", fixed = TRUE)[[1]][1]
  
  if (isTRUE(plasmode)) {
    # load the original dataset
    load(paste0(plasmode_folder, "AD_PT_plasmode_base_dataset.Rdata"))
    if (isTRUE(grepl("glm", yscenario))) {
      model_txt <- "glm"
    } else {
      model_txt <- "tree"
    }
    miss_fit_name <- paste0("miss_", model_txt, "_fit")
    ps_fit_name <- paste0("ps_", model_txt, "_fit")
    outcome_fit_name <- paste0(gsub(paste0(model_txt, "_"), "", yscenario), "_", model_txt, "_fit")
    # load the missing-data fit; it's called miss_glm_fit
    load(paste0(plasmode_folder, miss_fit_name, ".Rdata"))
    miss_fit <- eval(parse(text = miss_fit_name))
    # load the propensity score of treatment; it's called 
    load(paste0(plasmode_folder, ps_fit_name, ".Rdata"))
    ps_fit <- eval(parse(text = ps_fit_name))
    # load the outcome regression
    load(paste0(plasmode_folder, outcome_fit_name, ".Rdata"))
    outcome_fit <- eval(parse(text = outcome_fit_name))
  } else {
    outcome_model <- tx_model <- miss_model <- NULL
  }
  
  logrr <- vector("numeric", length = nreps)
  rd <- vector("numeric", length = nreps)
  logor <- vector("numeric", length = nreps)
  beta <- vector("numeric", length = nreps)
  logrr_zs <- vector("numeric", length = nreps)
  rd_zs <- vector("numeric", length = nreps)
  logor_zs <- vector("numeric", length = nreps)
  beta_zs <- vector("numeric", length = nreps)
  logrr_ws <- vector("numeric", length = nreps)
  rd_ws <- vector("numeric", length = nreps)
  logor_ws <- vector("numeric", length = nreps)
  if (isTRUE(plasmode)) {
    beta_ws <- t(replicate(3, vector("numeric", length = nreps)))
  } else {
    beta_ws <- vector("numeric", length = nreps)  
  }
  
  pop_logrr <- vector("numeric", length = nreps)
  pop_rd <- vector("numeric", length = nreps)
  pop_logor <- vector("numeric", length = nreps)
  pop_beta <- vector("numeric", length = nreps)
  pop_logrr_zs <- vector("numeric", length = nreps)
  pop_rd_zs <- vector("numeric", length = nreps)
  pop_logor_zs <- vector("numeric", length = nreps)
  pop_logrr_ws <- vector("numeric", length = nreps)
  pop_rd_ws <- vector("numeric", length = nreps)
  pop_logor_ws <- vector("numeric", length = nreps)
  if (isTRUE(plasmode)) {
    pop_beta_zs <- t(replicate(5, vector("numeric", length = nreps)))
    pop_beta_ws <- t(replicate(3, vector("numeric", length = nreps)))
  } else {
    pop_beta_zs <- vector("numeric", length = nreps)
    pop_beta_ws <- vector("numeric", length = nreps)  
  }
  for (j in 1:nreps) {
    # generate a dataset
    if (isTRUE(!plasmode)) {
      df <- gen_data(N = sample_size, XScenario = xscenario,
                     YScenario = yscenario, missScenario = mscenario,
                     lowcor = lowcor, midcor = midcor, highcor = highcor,
                     gencor = gencor)  
    } else {
      samp <- sample(x = 1:nrow(data), size = sample_size, replace = TRUE)
      df <- data[samp, ]
      ps_preds <- get_preds(fit = ps_fit, df = df)
      df$numEpiType <- as.numeric(runif(sample_size) < ps_preds)
      miss_preds <- get_preds(fit = miss_fit, df = df)
      df$is.complete <- as.numeric(runif(sample_size) < miss_preds)
      outcome_preds <- get_preds(fit = outcome_fit, df = df)
      df$Y <- as.numeric(runif(sample_size) < outcome_preds)
      df$AgeAtIndex_cat <- factor(
        case_when(
          df$AgeAtIndex <= 24 ~ 0,
          df$AgeAtIndex > 24 & df$AgeAtIndex <= 34 ~ 1,
          df$AgeAtIndex > 34 & df$AgeAtIndex <= 44 ~ 2,
          df$AgeAtIndex > 44 & df$AgeAtIndex <= 54 ~ 3,
          df$AgeAtIndex > 54 & df$AgeAtIndex <= 64 ~ 4,
          df$AgeAtIndex > 64 ~ 5
        ), levels = 0:5, labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+")
      )
    }
    
    if (y_only) {
      print(this_scenario)
      print(mean(df$Y))
      next
    }
    # datasets setting X = 0, X = 1 for everyone
    newdata_0 <- df
    newdata_1 <- df
    
    if (isTRUE(plasmode)) {
      newdata_0$numEpiType <- 0
      newdata_0$Y <- as.numeric(runif(sample_size) < get_preds(outcome_fit, newdata_0))
      newdata_1$numEpiType <- 1
      newdata_1$Y <- as.numeric(runif(sample_size) < get_preds(outcome_fit, newdata_1))
    } else {
      newdata_0$X <- 0
      newdata_0$Y <- YFunc(data = newdata_0, YScenario = yscenario)$Y
      newdata_1$X <- 1
      newdata_1$Y <- YFunc(data = newdata_1, YScenario = yscenario)$Y  
    }
    
    # datasets setting confounders to their mean value (or another sensible value) 
    # Z is an always-observed continuous/categorical confounder
    # W is a not-always-observed continuous/categorical confounder
    newdata_Zs_0 <- df
    newdata_Zs_1 <- df
    newdata_Ws_0 <- df
    newdata_Ws_1 <- df
    if (isTRUE(plasmode)) {
      # Z is age at index
      newdata_Zs_0$AgeAtIndex <- mean(df$AgeAtIndex)
      newdata_Zs_0$AgeAtIndexsq <- mean(df$AgeAtIndex) ^ 2
      newdata_Zs_0$AgeAtIndex_cat <- cut(newdata_Zs_0$AgeAtIndex,
                                         breaks = c(24, 34, 44, 54, 64, 200, 201), right = TRUE,
                                         labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))
      newdata_Zs_1$AgeAtIndex_cat <- cut(newdata_Zs_1$AgeAtIndex,
                                         breaks = c(24, 34, 44, 54, 64, 200, 201), right = TRUE,
                                         labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))
      newdata_Zs_1$AgeAtIndex <- mean(df$AgeAtIndex) + 1
      newdata_Zs_1$AgeAtIndexsq <- (mean(df$AgeAtIndex) + 1) ^ 2
      # W is item 9 score
      newdata_Ws_0$IndexTrtPHQ_item9_score <- mean(df$IndexTrtPHQ_item9_score)
      newdata_Ws_0$IndexTrtPHQ_item9_score_cat <- cut(newdata_Ws_0$IndexTrtPHQ_item9_score,
                                                      breaks = c(0, 1, 2, 3, 4), right = FALSE,
                                                      labels = c("None", "Some of the days", "Most of the days", "Nearly everyday"))
      newdata_Ws_1$IndexTrtPHQ_item9_score <- mean(df$IndexTrtPHQ_item9_score) + 1
      newdata_Ws_1$IndexTrtPHQ_item9_score_cat <- cut(newdata_Ws_1$IndexTrtPHQ_item9_score,
                                                      breaks = c(0, 1, 2, 3, 4), right = FALSE,
                                                      labels = c("None", "Some of the days", "Most of the days", "Nearly everyday"))
    } else {
      # for fully-synthetic data, Zs = mean(Zs), Zs = mean(Zs) + 1 for everyone (unless Y scenario is 4 or 4.1, in which case Zs is a binary variable)
      if (!(yscenario %in% c(4, 4.1))) {
        newdata_Zs_0$Zs <- mean(df$Zs)
        newdata_Zs_1$Zs <- mean(df$Zs) + 1
      } else { # function of Zs is to check if less than -1
        newdata_Zs_0$Zs <- -2
        newdata_Zs_1$Zs <- 0
      }
      newdata_Zs_0$Y <- YFunc(data = newdata_Zs_0, YScenario = yscenario)$Y
      newdata_Zs_1$Y <- YFunc(data = newdata_Zs_1, YScenario = yscenario)$Y
      # datasets setting Ws = mean(Ws), Ws = mean(Ws) + 1 for everyone 
      newdata_Ws_0$Ws <- mean(df$Ws)
      newdata_Ws_1$Ws <- mean(df$Ws) + 1
      newdata_Ws_0$Y <- YFunc(data = newdata_Ws_0, YScenario = yscenario)$Y
      newdata_Ws_1$Y <- YFunc(data = newdata_Ws_1, YScenario = yscenario)$Y
    }
    
    # Oracle values, i.e. using "true" outcome model, so based on observed outcomes. 
    e_y1 <- mean(newdata_1$Y)
    e_y0 <- mean(newdata_0$Y)
    logrr[j] <- log(e_y1 / e_y0)
    rd[j] <- e_y1 - e_y0
    logor[j] <- log((e_y1 / (1 - e_y1)) / (e_y0 / (1 - e_y0)))
    if (isTRUE(plasmode)) {
      beta[j] <- get_oracle_cOR(outcome_fit, "numEpiType")
    } else {
      beta[j] <- YFunc(data = newdata_0[, 1:10], YScenario = yscenario)$beta[1]  
    }
    
    e_y1_zs <- mean(newdata_Zs_1$Y)
    e_y0_zs <- mean(newdata_Zs_0$Y)
    logrr_zs[j] <- log(e_y1_zs / e_y0_zs)
    rd_zs[j] <- e_y1_zs - e_y0_zs
    logor_zs[j] <- log((e_y1_zs / (1 - e_y1_zs)) / (e_y0_zs / (1 - e_y0_zs)))
    if (isTRUE(plasmode)) {
      beta_zs[j] <- get_oracle_cOR(outcome_fit, "AgeAtIndex")
    } else {
      beta_zs[j] <- YFunc(data = newdata_Zs_0[, 1:10], YScenario = yscenario)$beta[2]  
    }
    
    e_y1_ws <- mean(newdata_Ws_1$Y)
    e_y0_ws <- mean(newdata_Ws_0$Y)
    logrr_ws[j] <- log(e_y1_ws / e_y0_ws)
    rd_ws[j] <- e_y1_ws - e_y0_ws
    logor_ws[j] <- log((e_y1_ws / (1 - e_y1_ws)) / (e_y0_ws / (1 - e_y0_ws)))
    if (isTRUE(plasmode)) {
      coef_names <- names(coefficients(outcome_fit))
      beta_ws[, j] <- get_oracle_cOR(outcome_fit, grepl("IndexTrtPHQ_item9_score_cat", coef_names) & !grepl(":", coef_names, fixed = TRUE))
    } else {
      beta_ws[j] <- YFunc(data = newdata_Ws_0[, 1:10], YScenario = yscenario)$beta[3]
    }   
    # Population values, i.e. using fitted values to assumed outcome model. 
    # For now, we're using the same outcome formula for each data generating mechanism
    # (couldn't use output of get_outcome_formulas, since it depended on W_obs, not W).
    # Depending on outcome model used in future, might need to change.
    # I am also not specifically considering the confounded model, where the outcome formula
    # does not include confounders W (to match what was done previously), but this is inconsistent.
    fit <- glm(as.formula(census_formula), family = "binomial", data = df)
    preds_0 <- predict(fit, newdata = newdata_0, type = "response", na.action = na.omit) 
    preds_1 <- predict(fit, newdata = newdata_1, type = "response", na.action = na.omit) 
    preds_0_zs <- predict(fit, newdata = newdata_Zs_0, type = "response", na.action = na.omit) 
    preds_1_zs <- predict(fit, newdata = newdata_Zs_1, type = "response", na.action = na.omit) 
    preds_0_ws <- predict(fit, newdata = newdata_Ws_0, type = "response", na.action = na.omit) 
    preds_1_ws <- predict(fit, newdata = newdata_Ws_1, type = "response", na.action = na.omit) 
    
    pop_e_y1 <- mean(preds_1)
    pop_e_y0 <- mean(preds_0)
    pop_logrr[j] <- log(pop_e_y1 / pop_e_y0)
    pop_rd[j] <- pop_e_y1 - pop_e_y0
    pop_logor[j] <- log((pop_e_y1 / (1 - pop_e_y1)) / (pop_e_y0 / (1 - pop_e_y0)))
    if (isTRUE(!plasmode)) {
      tx_name <- "X"
      zs_name <- "Zs"
      ws_name <- "Ws"
      pop_beta_ws[j] <- coefficients(fit)[ws_name]
      pop_beta_zs[j] <- coefficients(fit)[zs_name]
    } else {
      tx_name <- "numEpiType"
      census_coef_names <- names(coefficients(fit))
      zs_name <- census_coef_names[grepl("AgeAtIndex_cat", census_coef_names) & !grepl(":", census_coef_names, fixed = TRUE)]
      ws_name <- census_coef_names[grepl("IndexTrtPHQ_item9_score_cat", census_coef_names) & !grepl(":", census_coef_names, fixed = TRUE)]
      pop_beta_ws[, j] <- coefficients(fit)[ws_name]
      pop_beta_zs[, j] <- coefficients(fit)[zs_name]
    }
    pop_beta[j] <- coefficients(fit)[tx_name]
    
    pop_e_y1_zs <- mean(preds_1_zs)
    pop_e_y0_zs <- mean(preds_0_zs)
    pop_logrr_zs[j] <- log(pop_e_y1_zs / pop_e_y0_zs)
    pop_rd_zs[j] <- pop_e_y1_zs - pop_e_y0_zs
    pop_logor_zs[j] <- log((pop_e_y1_zs / (1 - pop_e_y1_zs)) / (pop_e_y0_zs / (1 - pop_e_y0_zs)))
    
    pop_e_y1_ws <- mean(preds_1_ws)
    pop_e_y0_ws <- mean(preds_0_ws)
    pop_logrr_ws[j] <- log(pop_e_y1_ws / pop_e_y0_ws)
    pop_rd_ws[j] <- pop_e_y1_ws - pop_e_y0_ws
    pop_logor_ws[j] <- log((pop_e_y1_ws / (1 - pop_e_y1_ws)) / (pop_e_y0_ws / (1 - pop_e_y0_ws)))
  }
  mc_id <- 1:nreps
  ws_mult <- length(ws_name)
  zs_mult <- length(zs_name)
  ws_name_marginal <- ifelse(isTRUE(plasmode), "IndexTrtPHQ_item9_score_cat", "Ws")
  zs_name_marginal <- ifelse(isTRUE(plasmode), "AgeAtIndex_cat", "Zs")
  zs_name_oracle <- ifelse(isTRUE(plasmode), "AgeAtIndex", "Zs")
  oracle_values_x <- tibble::tibble(
    mc_id = rep(mc_id, 4),
    estimand = rep(c("logRR", "RD", "logOR", "canonical_parameter"), each = nreps),
    type = "oracle",
    variable = tx_name,
    value = c(logrr, rd, logor, beta)
  )
  oracle_values_z <- tibble::tibble(
    mc_id = rep(mc_id, 4),
    estimand = rep(c("logRR", "RD", "logOR", "canonical_parameter"), each = nreps),
    type = "oracle",
    variable = zs_name_oracle,
    value = c(logrr_zs, rd_zs, logor_zs, beta_zs)
  )
  oracle_values_w_marginal <- tibble::tibble(
    mc_id = rep(mc_id, 3),
    estimand = rep(c("logRR", "RD", "logOR"), each = nreps),
    type = "oracle",
    variable = ws_name_marginal,
    value = c(logrr_ws, rd_ws, logor_ws)
  )
  oracle_values_w_conditional <- tibble::tibble(
    mc_id = rep(mc_id, each = ws_mult),
    estimand = "canonical_parameter",
    type = "oracle",
    variable = rep(ws_name, nreps),
    value = as.numeric(beta_ws)
  )
  oracle_values_w <- dplyr::bind_rows(oracle_values_w_marginal, oracle_values_w_conditional)
  pop_values_x <- tibble::tibble(
    mc_id = rep(mc_id, 4),
    estimand = rep(c("logRR", "RD", "logOR", "canonical_parameter"), each = nreps),
    type = "census",
    variable = tx_name,
    value = c(pop_logrr, pop_rd, pop_logor, pop_beta)
  )
  pop_values_z_marginal <- tibble::tibble(
    mc_id = rep(mc_id, 3),
    estimand = rep(c("logRR", "RD", "logOR"), each = nreps),
    type = "census",
    variable = zs_name_marginal,
    value = c(pop_logrr_zs, pop_rd_zs, pop_logor_zs)
  )
  pop_values_z_conditional <- tibble::tibble(
    mc_id = rep(mc_id, each = zs_mult),
    estimand = "canonical_parameter",
    type = "census",
    variable = rep(zs_name, nreps),
    value = as.numeric(pop_beta_zs)
  )
  pop_values_z <- dplyr::bind_rows(pop_values_z_marginal, pop_values_z_conditional)
  pop_values_w_marginal <- tibble::tibble(
    mc_id = rep(mc_id, 3),
    estimand = rep(c("logRR", "RD", "logOR"), each = nreps),
    type = "census",
    variable = ws_name_marginal,
    value = c(pop_logrr_ws, pop_rd_ws, pop_logor_ws)
  )
  pop_values_w_conditional <- tibble::tibble(
    mc_id = rep(mc_id, each = ws_mult),
    estimand = "canonical_parameter",
    type = "census",
    variable = rep(ws_name, nreps),
    value = as.numeric(pop_beta_ws)
  )
  pop_values_w <- dplyr::bind_rows(pop_values_w_marginal, pop_values_w_conditional)
  result <- dplyr::bind_rows(oracle_values_x, oracle_values_w, oracle_values_z,
                             pop_values_x, pop_values_w, pop_values_z) %>% 
    mutate(xscenario = xscenario, yscenario = yscenario, mscenario = mscenario)
  return(result)
}
