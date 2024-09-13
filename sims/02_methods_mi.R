# Methods - MI

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 20230725       bdw       harmonize methods into a single file
# 20230802			 bdw       break out into design, mi, tmle files
# 20230807       bdw       use coefficient_of_interest
# 20230816       bdw       use a single method, defaults to 20 imputations, using outcome
# 20230912       bdw       change return to return all coefs from MI
# 20231115       bdw       specify default vs specific method for MI
# 20240111       eaj       add RF to the procedure list
# 20240321       cak       makes interactions in outcome model more congenial with imputation model

# MI analysis ------------------------------------------------------------------

# @param data the dataset
# @param outcome_name the outcome name
# @param outcome_formula the formula of interest
# @param n_imp the number of imputations
# @param maxit the number of iterations
# @param coefficient_of_interest the coefficient name of interest
# @param method the method (e.g., "pmm")
# @param missing_indicator the name of the missingness indicator
# @param fam the outcome regression family
# @return a set of results for different mice tuning parameter choices
run_mice <- function(data = NULL, outcome_name = "Y", outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", n_imp = 10, maxit = 20, coefficient_of_interest = "X",
                     method = "default", missing_indicator = "is.complete", fam = "binomial", plasmode = FALSE) {
  data.mi <- data[, !grepl(missing_indicator, names(data))]
  init <- mice::mice(data.mi, maxit = 0, print = F)
  pred.matrix <- init$predictorMatrix
  
  
  #### ADDED 3/18/2024 BY CHLOE: CODE CHUNK TO MODIFY FORMULAS FOR PLASMODE WITH INTERACTIONS IN OUTCOME MODEL #####
  # Find variables with missing values.
  # Change formulas just for those, swapping outcome and variable being imputed.
  # This code chunk should not change anything if there are no interactions in the model.
  
  # Finding names of variables subject to missingness.
  have_missing <- names(init$nmis)[which(init$nmis>0)]
  
  for (i in 1:length(have_missing)){
    formula_new <- gsub(have_missing[i],outcome_name, outcome_formula)
    formula_new <- as.formula(gsub(paste0(outcome_name," ~"),paste0(have_missing[i]," ~"),formula_new))
    index <- which(names(init$formulas)==have_missing[i])
    init$formulas[[index]] <- formula_new
  }
  
  formulas_use <- init$formulas
  
  ##################################################################################################################
  
  pred.matrix.excludeY <- pred.matrix
  pred.matrix.excludeY[, outcome_name] <- 0
  
  #Edit from Eric, 11 January: allow different procedure names based on method used
  # Edit from Chloe, March 21, 2024: now use formulas specified above rather than predictor matrix (used previously)
  # for plasmode.
  
  if (plasmode == TRUE){
    if (method == "default") {
      mice_result <- mice::mice(
        data.mi, formulas=formulas_use, m = n_imp, maxit = maxit, 
        print = FALSE
      )
      procedure_value <- "MICE"
    } else {
      mice_result <- mice::mice(
        data.mi, formulas=formulas_use, m = n_imp, method = method, 
        maxit = maxit, print = FALSE
      )
      procedure_value <- "RF"
    }
    } else{
      if (method == "default") {
        mice_result <- mice::mice(
          data.mi, predictorMatrix = pred.matrix, m = n_imp, maxit = maxit, 
          print = FALSE
        )
        procedure_value <- "MICE"
      } else {
        mice_result <- mice::mice(
          data.mi, predictorMatrix = pred.matrix, m = n_imp, method = method, 
          maxit = maxit, print = FALSE
        )
        procedure_value <- "RF"
      }
  }
  
  # Chloe 11/8/2023: moved this out of results_mi_all so that I didn't have to create it twice
  # (once for MI results, once for object for marginal estimation).
  if (fam != "survival") {
    regress_mi <- with(
      # mice_data,
      mice_result,
      glm(as.formula(outcome_formula), family = fam)
    )
  }
  
  res_mi <- results_mi(regress_mi=regress_mi, outcome_formula = outcome_formula, coefficient_of_interest = coefficient_of_interest, fam = fam)
  
  # res_mi_all <- results_mi_all(mice_result, outcome_formula = outcome_formula, fam = fam)
  res_mi_all <- results_mi_all(regress_mi=regress_mi, outcome_formula = outcome_formula, fam = fam)
  
  names(res_mi_all) <- c("Estimate", "SE")
  conditional_est <- data.frame(
    "procedure" = procedure_value,
    "estimand" = "canonical_parameter",
    "est" = res_mi[, 1],
    "SE" = res_mi[, 2]
  )
  # marginal_est <- data.frame("procedure" = procedure_value, marginal_results_mi(mice_data = mice_result, outcome_formula = outcome_formula, coefficient_of_interest = coefficient_of_interest, fam = fam))
  marginal_params <- marginal_effects_new(fit=regress_mi, weights=NULL, coefficient_of_interest = coefficient_of_interest)
  marginal_est <- data.frame("procedure" = procedure_value,marginal_params)
  coefs <- data.frame("procedure" = procedure_value, res_mi_all, "test_statistic" = NA, "p_value" = NA)
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
# results_mi <- function(mice_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", coefficient_of_interest = "X", fam = "binomial") {
results_mi <- function(regress_mi, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", coefficient_of_interest = "X", fam = "binomial") {
  mi_results <- results_mi_all(regress_mi = regress_mi, outcome_formula = outcome_formula, 
                               fam = fam)
  return(mi_results[coefficient_of_interest, , drop = FALSE])
}
# obtain marginal estimates
# @param mice_data the multiply-imputed data (a mids object)
# @param outcome_formula the formula of interest
# @param coefficient_of_interest the name of the coefficient of interest
# @param fam the outcome regression family
# @return the point estimate and SE of the marginal effects, after applying Rubin's rules
marginal_results_mi <- function(mice_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", coefficient_of_interest = "X", fam = "binomial") {
  if (fam != "survival") {
    regress_mi <- with(
      mice_data,
      glm(as.formula(outcome_formula), family = fam)
    )
    
  }
  completed_data <- mice::complete(mice_data, action = "all")
  all_marginal_ests <- do.call(rbind, lapply(as.list(1:mice_data$m), function(i) {
    data.frame("m" = i, get_marginal_estimates(fit = regress_mi$analyses[[i]], data = completed_data[[i]], estimator = "MI"))
  }))
  # apply Rubin's rules
  marginal_ests_mi <- all_marginal_ests %>%
    group_by(estimand) %>%
    summarize(SE = sqrt(
      mean(SE^2) + (1 + 1 / mice_data$m) * var(est)
    ), est = mean(est)) %>%
    select(estimand, est, SE)
  return(as.data.frame(marginal_ests_mi))
}
# get all coefficients
# @param mice_data the multiply-imputed data (a mids object)
# @param outcome_formula the formula of interest
# @param fam the outcome regression family
# @return the point estimates and SEs from a regression model

# Chloe 11/8/2023: I changed this around, so that the fit occurs within the primary function above.
# results_mi_all <- function(mice_data, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", fam = "binomial") {
results_mi_all <- function(regress_mi, outcome_formula = "Y ~ X + Zs + Zw + Ww_obs + Ws_obs", fam = "binomial") {
  
  # Moved into main MICE function.
  # if (fam != "survival") {
  #   regress_mi <- with(
  #   mice_data,
  #   glm(as.formula(outcome_formula), family = fam)
  # )
  # }
  
  pooled_est <- mice::pool(regress_mi) # Default rule is "rubin1987", try calculating by hand in future.
  mi_results <- summary(pooled_est)
  rownames(mi_results) <- mi_results[, 1]
  return(mi_results[, c("estimate", "std.error")])
}