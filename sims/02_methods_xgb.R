# Methods - XGB

# NOTES
# Chloe 11/16/2023: getting marginal estimates using marginaleffects package;
# using Eric's trick, replacing MICE dummy object datasets with XGBoost generated datasets.

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 19 Feb 2024   EAJ         Changed run_xgb to have initial.fac="sample"

# XGB analysis ------------------------------------------------------------------

run_xgb <- function(data = NULL, outcome_name = "Y", 
                    outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", 
                    n_imp = 25, coefficient_of_interest = "X",
                    missing_indicator = "is.complete", fam = "binomial",
                    initial.fac="sample") {
  data.mi <- data[, !grepl(missing_indicator, names(data))]
  init <- mice::mice(data.mi, maxit = 0, print = F)
  
  # Chloe 11/16/2023: Creating dummy MICE object; 
  # imputed datasets to be replaced with imputed datasets from XGBoost later.
  dummy.mice <- mice::mice(data.mi, m=n_imp, maxit = 1, print = F)
  # And index of missing observations to be used later.
  missing_index <- !complete.cases(data)
  missing_cols <- which(colSums(is.na(data)) > 0)

  xb.cv <- mixgb_cv(data.mi, verbose=FALSE)
  xgb_result <- mixgb(data.mi, m=n_imp, nrounds=xb.cv$best.nrounds) # Imputed datasets.
  
  res_xgb <- results_xgb(xgb_result, outcome_formula = outcome_formula, 
                         coefficient_of_interest = coefficient_of_interest, fam = fam)
  res_xgb_all <- results_xgb_all(xgb_result, outcome_formula = outcome_formula, fam = fam)
  
  names(res_xgb_all) <- c("Estimate", "SE")
  conditional_est <- data.frame(
    "procedure" = "XGB",
    "estimand" = "canonical_parameter",
    "est" = res_xgb[, 1],
    "SE" = res_xgb[, 2]
  )
  marginal_est <- data.frame("procedure" = "XGB", marginal_results_xgb(xgb_data = xgb_result, 
                                                                       outcome_formula = outcome_formula, 
                                                                       coefficient_of_interest = coefficient_of_interest, 
                                                                       fam = fam, dummy.mice=dummy.mice, 
                                                                       missing_index=missing_index,
                                                                       missing_cols = missing_cols))
  coefs <- data.frame("procedure" = "XGB", res_xgb_all, "test_statistic" = NA, "p_value" = NA)
  return(list(
    "results" = rbind.data.frame(conditional_est, marginal_est),
    "coefs" = cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL))
  ))
}

# helper functions for MI
# get results for the specified coefficient of interest
# @param mice_data the multiply-imputed data (a mids object)
# @param outcome_formula the formula of interest
# @param coefficient_of_interest the name of the coefficient of interest
# @param fam the outcome regression family
# @return the point estimate and SE from a regression model
results_xgb <- function(xgb_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", coefficient_of_interest = "X", fam = "binomial") {
  xgb_results <- results_xgb_all(xgb_data = xgb_data, outcome_formula = outcome_formula, 
                                 fam = fam)
  return(xgb_results[coefficient_of_interest, , drop = FALSE])
}
# obtain marginal estimates
# @param mice_data the multiply-imputed data (a mids object)
# @param outcome_formula the formula of interest
# @param coefficient_of_interest the name of the coefficient of interest
# @param fam the outcome regression family
# @param missing_index the row indices for missing data
# @param missing_cols the column indices for missing data
# @return the point estimate and SE of the marginal effects, after applying Rubin's rules
marginal_results_xgb <- function(xgb_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", coefficient_of_interest = "X", fam = "binomial",dummy.mice = NULL, missing_index = NULL, missing_cols = NULL) {
  
  ## NEW 11/16/2023: replacing imputed datasets in MICE dummy object with XGBoost-imputed datasets.
  xgb_mice <- dummy.mice
  missing_index <- !complete.cases(dummy.mice$data)
  missing_cols <- which(colSums(is.na(dummy.mice$data)) > 0)
  for (j in 1:length(missing_cols)) {
    for (i in 1:length(xgb_data)) {
      xgb_mice$imp[[missing_cols[j]]][, i] <- xgb_data[[i]][[missing_cols[j]]][missing_index]  
    }
  }
  
  regress_mice_xgb <- with(xgb_mice,glm(as.formula(outcome_formula), family = fam))
  # Testing new code: 
  # pool(regress_mice_xgb): matches original conditional results.
  
  marginal_params <- marginal_effects_new(fit=regress_mice_xgb, weights=NULL, coefficient_of_interest = coefficient_of_interest)
  # Will differ slightly from original for OR/RR since we're now looking log scale now.
  
  ## Original code to get results.
  # if (fam != "survival") {
  #   regress_xgb <- foreach(i=(1:length(xgb_data))) %do% {
  #     glm(as.formula(outcome_formula),family=fam, data=xgb_data[[i]])  
  #   }
  # }
  # #going to use get_marginal_estimates for MI, since XGB also uses Rubin's rules
  # 
  # imps <- length(xgb_data)
  # all_marginal_ests <- foreach(i=1:imps, .combine=rbind) %do% {
  #   data.frame("m" = i, get_marginal_estimates(fit = regress_xgb[[i]], data = xgb_data[[i]], estimator = "MI"))
  # }
  # # apply Rubin's rules
  # marginal_ests_xb <- all_marginal_ests %>%
  #   group_by(estimand) %>%
  #   summarize(SE = sqrt(
  #     mean(SE^2) + (1 + 1 / imps) * var(est)
  #   ), est = mean(est)) %>%
  #   select(estimand, est, SE)
  # return(as.data.frame(marginal_ests_xb))
}
# get all coefficients
# @param xgb_data the multiply-imputed data (a list with n_imp dataframes)
# @param outcome_formula the formula of interest
# @param fam the outcome regression family
# @return the point estimates and SEs from a regression model
results_xgb_all <- function(xgb_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", fam = "binomial") {
  if (fam != "survival") {
    regress_xgb <- foreach(i=(1:length(xgb_data))) %do% {
      glm(as.formula(outcome_formula),family=fam, data=xgb_data[[i]])  
    }
  }
  
  pooled_est <- mice::pool(regress_xgb) # Default rule is "rubin1987", try calculating by hand in future.
  xgb_results <- summary(pooled_est)
  rownames(xgb_results) <- xgb_results[, 1]
  return(xgb_results[, c("estimate", "std.error")])
}