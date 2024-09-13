# Methods - TMLE

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 20230725       bdw       harmonize methods into a single file
# 20230802			 bdw       break out into design, mi, tmle files
# 20230811       bdw       fix EIFs and move utility functions from utils to here
#                          fix predicting missingness -- SL hard-codes "Y" so have to trick it
# 20230815       bdw       add subcalTMLE option
# 20230818       bdw       remove quasibinomial family
# 20230907       bdw       use subcalTMLE and tmle package > 2.0, implement cOR estimation
# 20240502       bdw       updates from Susan: possibility to change default library for rare outcomes,
#                          augment the dataset for better missing-data probability estimation
# 20240503       bdw       estimate prob of tx for augmentation term

# TMLE analysis ----------------------------------------------------------------
# run subcalTMLE
# @param data the dataset (assumes in order and has is.complete variable)
# @param outcome_name the name of the outcome
# @param tx_name the name of the treatment assignment variable
# @param outcome_formula the outcome regression formula, only used if the outcome regression library is only glm
# @param outcome_formula_factor the outcome regression formula with explicit factor dummy variables, for TMLE-M plasmode
# @param tx_formula the propensity score formula, only used if the propensity score library is only glm
# @param tx_formula_factor the propensity score formula, only used if the propensity score library is only glm, for TMLE-M plasmode
# @param miss_formula the missing-data formula, only used if the missing data library is only glm
# @param miss_formula_factor the missing-data formula with explicit factor dummy variables, used if these are passed in for the MSM
# @param delta_name the name of the variable used to designate complete observations
# @param phase1_covars the names of the covariates measured on everyone
# @param phase2_covars the names of the covariates measured only on a subset
# @param g_learner_lib the SL library for estimating the (treatment assignment) propensity score
# @param miss_learner_lib the SL library for estimating the probability of missingness
# @param q_learner_lib the SL library for estimating the outcome regression
# @param K the number of cross-validation folds for the Super Learner
# @param procedure the type of procedure, e.g., TMLE-M
# @param fam the outcome regression family
# @param augment_w should we augment the variables for initial estimation of the missing-data probability?
# @param rare_outcome should we use a special SL library for rare outcomes?
# @return the TMLE estimators and SEs for the risk difference and odds ratio
run_subcaltmle <- function(data, g_learner_lib = c("SL.glm"), miss_learner_lib = c("SL.glm"),
                           q_learner_lib = c("SL.glm"), K = 5,
                           outcome_name = "Y", tx_name = "X",
                           outcome_formula = "Y ~ X + Zs + Zw + Ws_obs + Ww_obs",
                           outcome_formula_factor = "Y ~ X + Zs + Zw + Ws_obs + Ww_obs",
                           tx_formula = "X ~ Zs + Zw + Ws_obs + Ww_obs + As + Aw",
                           tx_formula_factor = "X ~ Zs + Zw + Ws_obs + Ww_obs + As + Aw",
                           miss_formula = "is.complete ~ X + Zs + Zw",
                           miss_formula_factor = "is.complete ~ X + Zs + Zw",
                           delta_name = "is.complete",
                           phase1_covars = c("Y", "X", "Zs", "Zw"),
                           condition_on_auxiliary = FALSE,
                           phase2_covars = c("Ws", "Ww"), procedure = "TMLE-M", fam = "binomial",
                           augment_w = FALSE, rare_outcome = FALSE) {
  any_factor <- any(sapply(data, is.factor))
  # if fitting a glm for any portion, use the passed-in formula
  # if there are factor variables, use the formula for factor variables
  Qform <- switch(as.numeric(all(q_learner_lib == "SL.glm")) + 1, NULL, 
                  gsub(tx_name, "A", 
                       switch(as.numeric(any_factor) + 1, outcome_formula, outcome_formula_factor)))
  gform_init <- switch(as.numeric(any_factor) + 1, tx_formula, tx_formula_factor)
  gform_init2 <- gsub("A(?![g|t])", "X", gform_init, perl = TRUE)
  gform <- switch(as.numeric(all(g_learner_lib == "SL.glm")) + 1, NULL, 
                  gsub(paste0(tx_name, " ~ "), "A ~ ", gform_init2))
  piform <- switch(as.numeric(all(miss_learner_lib == "SL.glm")) + 1, NULL,
                   gsub(outcome_name, "real_Y", gsub(delta_name, "Delta.W", gsub(tx_name, "A", miss_formula))))  
  
  # get phase 2 data, w_sub, and V (for MSM)
  phase_2_data <- data[, phase2_covars, drop = FALSE]
  w_sub <- phase_2_data[complete.cases(phase_2_data), , drop = FALSE]
  phase1_covars_no_AY <- phase1_covars[!(phase1_covars %in% outcome_name) & 
                                         !(phase1_covars %in% tx_name)]
  W <- data[, phase1_covars_no_AY, drop = FALSE]
  v_vars <- c(phase1_covars_no_AY, phase2_covars)
  if (!condition_on_auxiliary) {
    v_vars <- v_vars[!(grepl("A", v_vars, ignore.case = FALSE) & !grepl("Age", v_vars, ignore.case = FALSE))]
  }
  V <- data[, v_vars, drop = FALSE]
  names(W) <- gsub("A(?![g|t])", "X", names(W), perl = TRUE)
  names(V) <- gsub("A(?![g|t])", "X", names(V), perl = TRUE)
  cond_set_names <- "W"
  if (grepl(tx_name, miss_formula)) {
    cond_set_names <- c(cond_set_names, "A")
  }
  if (grepl(outcome_name, miss_formula)) {
    cond_set_names <- c(cond_set_names, "Y")
  }
  # estimate marginal parameters
  marg_subcalTMLE <- subcalTMLE(
    Y = data[[outcome_name]], A = data[[tx_name]], W = W,
    Delta.W = data[[delta_name]], W.sub = w_sub, pi.SL.library = miss_learner_lib, pi.discreteSL = FALSE,
    condSetNames = cond_set_names, Q.family = fam, V.pi = K, V.Q = K, V.g = K,
    Q.SL.library = q_learner_lib, g.SL.library = g_learner_lib, verbose = FALSE,
    Qform = Qform, gform = gform, piform = piform, augmentW = augment_w,
    rareOutcome = rare_outcome
  )
  if (isTRUE(is.null(marg_subcalTMLE$tmle))) {
    marginal_est <- data.frame(
      "procedure" = procedure,
      "estimand" = c("RD", "OR", "RR"),
      "est" = rep(NA, 3),
      "SE" = rep(NA, 3)
    )
  } else {
    marginal_est <- data.frame(
      "procedure" = procedure,
      "estimand" = c("RD", "OR", "RR"),
      "est" = c(
        marg_subcalTMLE$tmle$estimates$ATE$psi, 
        marg_subcalTMLE$tmle$estimates$OR$psi,
        marg_subcalTMLE$tmle$estimates$RR$psi
      ),
      "SE" = c(
        sqrt(as.numeric(marg_subcalTMLE$tmle$estimates$ATE$var.psi)), 
        sqrt(as.numeric(marg_subcalTMLE$tmle$estimates$OR$var.log.psi)),
        sqrt(as.numeric(marg_subcalTMLE$tmle$estimates$RR$var.log.psi))
      )
    )
  }
  # need to convert factors to dummy variables for MSM
  msm_data_list <- make_dummy_variables(data = data, outcome_formula = outcome_formula,
                                        tx_formula = tx_formula,
                                        W = W, V = V, w_sub = w_sub,
                                        phase1_covars_no_AY = phase1_covars_no_AY,
                                        phase2_covars = phase2_covars)
  W_msm <- msm_data_list$W
  V_msm <- msm_data_list$V
  W_sub_msm <- msm_data_list$W_sub
  msm_form <- paste0("A + ",  paste0(names(V_msm), collapse = " + "))
  piform_msm <- switch(as.numeric(all(miss_learner_lib == "SL.glm")) + 1, NULL,
                       gsub(delta_name, "Delta.W", gsub(tx_name, "A", miss_formula_factor)))
  cond_subcalTMLE <- subcalTMLEMSM(
    Y = data[[outcome_name]], A = data[[tx_name]], W = W_msm, V = V_msm, MSM = msm_form,
    Delta.W = data[[delta_name]], W.sub = W_sub_msm, pi.SL.library = miss_learner_lib, pi.discreteSL = FALSE,
    condSetNames = cond_set_names, Q.family = fam, V.pi = K, V.Q = K, V.g = K,
    Q.SL.library = q_learner_lib, g.SL.library = g_learner_lib, verbose = FALSE,
    Qform = Qform, gform = gform, piform = piform_msm, augmentW = augment_w,
    rareOutcome = rare_outcome
  )
  if (isTRUE(is.null(cond_subcalTMLE$tmleMSM))) {
    n_params <- ncol(V_msm) + 2
    coefs <- data.frame("procedure" = procedure, "Estimate" = rep(NA, n_params),
                        "SE" = rep(NA, n_params), "test_statistic" = NA, p_value = NA)
    rownames(coefs) <- c("(Intercept)", "numEpiType", names(V_msm))
  } else {
    summ <- summary(cond_subcalTMLE$tmleMSM)$estimates[, 1:3]
    coefs <- data.frame("procedure" = procedure, summ[, 1:2], "test_statistic" = NA,
                        "p_value" = summ[, 3])
  }
  names(coefs) <- c("procedure", "Estimate", "SE", "test_statistic", "p_value")
  final_coefs <- cbind.data.frame("Variable" = rownames(coefs), data.frame(coefs, row.names = NULL))
  final_coefs$Variable <- gsub("A", tx_name, final_coefs$Variable)
  conditional_est <- data.frame(
    "procedure" = procedure,
    "estimand" = c("cOR"),
    "est" = final_coefs$Estimate[final_coefs$Variable == tx_name],
    "SE" = final_coefs$SE[final_coefs$Variable == tx_name]
  )
  return(list(
    "results" = rbind.data.frame(conditional_est, marginal_est),
    "coefs" = final_coefs
  ))
}

# utility functions for subcalTMLE ---------------------------------------------
#--------------------------------------------------------------------------
# Susan Gruber
# sgruber@TLrevolution.com
# August 10, 2023

# functions for implementing a TMLE for subset calibration
# (based  on two-stage TMLE Rose&vanderLaan, 2011)
# that evaluates a point treatment effect or the marginal mean outcome for a
# single arm study using a 2-stage sampling design
#
# For a two-arm study the full data consists of baseline covariates X, and potential outcomes Y0, Y1
# Observed data consists of baseline covariates (W, Delta.W, Delta.W W'),
# treatment A, outcome (Delta, Delta Y), where covariates X can be divided into two sets: W measured
# on everyone, and W' measured on only a subset of individuals.  Delta.W is an indicator of being in
# the selected subset. A is a unary (single-arm study) or binary treatment indicator, Delta indicates
# whether the outcome is measured and DeltaY is the outcome value (when Delta = 1), otherwise missing.

# The functions defined below take an initial observed dataset and transform it into a suitable
# representation for applying targeted maximum likelihood estimation (TMLE) using the tmle package (>= 2.0.0).
# Marginal parameter estimates (ATE/RD, RR, OR) can be obtained by calling the tmle function.
# Conditional parameter estimates are obtained by calling the tmleMSM function. Note this requires
# specifying the desired MSM. Data on observations in the subset are passed in to these functions,
# along with sampling weights that account for the possibly biased selection into the sample.
#
# Step 1. Simple imputation of miscellaneous missing values in W + augment W with binary indicators of
# 		  which values were present/absent
# Step 2. Construct inverse probability of censoring weights to pass into the tmle or tmleMSM functions
# 	2a. Model conditional probability of being included in the subset using SuperLearner or logistic regression.
# 	2b. Construct IP weights
# Step 3. Call the tmle or tmleMSM function to evaluate the parameter, estimate variance and CI bounds
# 	using influence curve-based inference.  Optionally, the tmle function can also provide targeted
# bootstrap-based variance estimates (recommended when there are near positivity violations).
#--------------------------------------------------------------------------

# For NOW - no missing data allowed in W, W.sub
#-------------------------------
# Utilities
summary.subcal <- function(object, ...) {
  if (!is.null(object$coef)) {
    picoef <- object$coef
    if (inherits(picoef, "matrix")) {
      piterms <- colnames(picoef)
    } else {
      piterms <- names(picoef)
    }
    pimodel <- paste("Delta.W ~ 1")
    if (length(piterms) > 1) {
      pimodel <- paste("Delta.W ~ ", paste(piterms, collapse = " + "))
    }
  } else {
    pimodel <- piterms <- picoef <- NULL
  }
  return(list(
    pimodel = pimodel, piterms = piterms, picoef = picoef, pitype = object$type,
    pidiscreteSL = object$discreteSL
  ))
}

print.summary.subcalTMLE <- function(x, ...) {
  if (inherits(x, "summary.subcalTMLE")) {
    cat("\n Estimation of Pi (subset sampling mechanism)\n")
    cat("\t Procedure:", x$subcal$pitype)
    if (!(is.null(x$subcal$pidiscreteSL))) {
      if (x$subcal$pidiscreteSL) {
        cat(", discrete")
      } else {
        cat(", ensemble")
      }
    }
    if (!(is.null(x$subcal$piAUC))) {
      cat("\t Empirical AUC =", round(x$subcal$piAUC, 4), "\n")
    }
    cat("\n")
    if (!(is.na(x$subcal$picoef[1]))) {
      cat("\t Model:\n\t\t", x$subcal$pimodel, "\n")
      cat("\n\t Coefficients: \n")
      terms <- sprintf("%15s", x$subcal$piterms)
      extra <- ifelse(x$subcal$picoef >= 0, "  ", " ")
      for (i in 1:length(x$subcal$picoef)) {
        cat("\t", terms[i], extra[i], x$subcal$picoef[i], "\n")
      }
    }
    cat("\n")
    print(x$tmle)
  }
}

print.subcalTMLE <- function(x, ...) {
  cat("Subset calibration TMLE\n")
  if (inherits(x, "subcalTMLE")) {
    print(x$tmle)
  }
}

summary.subcalTMLE <- function(x, ...) {
  # summary for estimating Pi here
  if (inherits(x, "subcalTMLE")) {
    sum.subcalTMLE <- list()
    sum.subcalTMLE$subcal <- summary.subcal(x$subcal)
    sum.subcalTMLE$tmle <- summary(x$tmle)
    class(sum.subcalTMLE) <- "summary.subcalTMLE"
    return(sum.subcalTMLE)
  }
}

#-----------------------------------------
# Set the number of cross-validation
# folds as a function of n.effective
# See Phillips 2023 doi.org/10.1093/ije/dyad023
#-----------------------------------------
setV <- function(n.effective) {
  if (n.effective <= 30) {
    V <- n.effective
  } else if (n.effective <= 500) {
    V <- 20
  } else if (n.effective <= 1000) {
    V <- 10
  } else if (n.effective <= 10000) {
    V <- 5
  } else {
    V <- 2
  }
  return(V)
}


#----------------
# construct column names
#-----------------
.getColNames <- function(condSetNames, Wnames, Vnames = NULL) {
  orig.colnames <- NULL
  for (i in 1:length(condSetNames)) {
    # if (condSetNames[i] %in% c("A", "Y")) {
    if (condSetNames[i] == "A") {
      orig.colnames <- c(orig.colnames, condSetNames[i])
    } else if (condSetNames[i] == "Y") {
      orig.colnames <- c(orig.colnames, "real_Y")
    } else if (condSetNames[i] == "W") {
      orig.colnames <- c(orig.colnames, Wnames)
    } else if (condSetNames[i] == "V") {
      orig.colnames <- c(orig.colnames, Vnames)
    }
  }
  return(orig.colnames)
}

#-------------------------
# Function subcalTMLE
# purpose: subset calibration TMLE
# Arguments:
#  	Y : outcome of interest (missingness allowed)
# 	A : binary treatment indicator
# 	W : matrix or data.frame of covariates measured on entire population
# 	Delta.W : Indicator of inclusion in subset with additional information
# 	W.sub : matrix or data.frame of covariates measured in subset population
# 			  (in same order as the corresponding rows in W. DO NOT include
# 			   rows for subjects not included in the subset)
# 	Z : Only required when the goal is estimating a controlled direct effect
# 		mediated by a binary variable,Z,  between treatment and outcome
# 	Delta : binary indicator that outcome Y is observed
# 	pi : optional vector of sampling probabilities
# 	piform : optional parametric regression model for estimating pi
# 	pi.SL.library : optional SL library specification for estimating pi (ignored when piform or pi is provided)
# 	V.pi : optional number of cross-validation folds for super learning (ignored when piform or pi is provided)
# 	pi.discreteSL: flag to indicate whether to use ensemble or discrete super learning (ignored when piform or pi is provided)
# 	condSet : variables to condition on when estimating pi. Covariates in W are always included. Can optionally also
# 			   condition on A and/or Y.
# 	id : optional indicator of independent units of observation
# 	Q.family : "gaussian" or "binomial", used to evaluate the TMLE
#   augmentW : TRUE to augment W with predicted values when A = 0 and A = 1
#			   from fitting the outcome regression conditional on A and W on all observations 
# 	rareOutcome: When true calls TMLE with  V.Q = 20, Q.discreteSL = TRUE,
#   		   Q.SL.library = c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2")
# 	verbose : set to TRUE to print out informative messages
# 	...	other arugments passed to the tmle function (see tmle help files)
# Return value
# 	object of class "subcalTMLE" (a list):
# 	 subcal - a list of details on estimating the conditional sampling probabilities
# 		pi - predicted probabilities
# 		coef - coefficients in the parametric model or SL ensemble weights
# 		type - "user-supplied regression formula" or "Super Learner"
# 		"discreteSL" - FALSE when ensemble SL used, TRUE when best single algorithm selected
# 	tmle - results returned from the call to the tmle function with
# 		   observation weights equal to 1/pi. (class "tmle", see tmle package documentation)
#-------------------------
subcalTMLE <- function(Y, A, W, Delta.W, W.sub, Z = NULL, Delta = rep(1, length(Y)), 
                       pi = NULL, piform = NULL, 
                       pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet", 
                                         "tmle.SL.dbarts.k.5"), 
                       V.pi = 10, pi.discreteSL = TRUE, condSetNames = c("A", "W"), 
                       id = NULL, Q.family = "gaussian", augmentW = TRUE, 
                       rareOutcome = FALSE, verbose = FALSE, ...) {
  if (!requireNamespace("tmle", versionCheck = ">=2.0", quietly = FALSE)) {
    stop("Loading required tmle package (>=2.0) failed", call. = FALSE)
  }
  
  if (is.null(id)) {
    id <- 1:length(Y)
  }
  
  if (is.vector(W.sub)) {
    W.sub <- as.matrix(W.sub)
    colnames(W.sub) <- "W.sub"
  }
  
  # if augmenting W call TMLE to get initial estimate of Q in full sample	
  # do this through the tmle function to get cross-validated initial estimates
  # that are not overfit to the data and are appropriately bounded
  if (augmentW) {
    L <- list(...)
    n.effective <- sum(Delta)
    if (Q.family == "binomial") {
      n.effective <- min(c(table(Y[Delta == 1]) * 5, n.effective))
    }
    V.Q <- setV(n.effective)
    d.g <- data.frame(Y = Y, mget(condSetNames))
    W.g <- tmle:::estimateG(d = d.g, gform = "A ~ .", SL.library = "SL.glm", 
                            outcome = "A", message = "treatment assigment probability",
                            g1W = NULL, V = V.pi, discreteSL = TRUE,
                            id = id, obsWeights = rep(1, nrow(W)), verbose = verbose)
    W.Q <- tmle(Y, A, W, Z=Z, Delta = Delta, id = id, family = Q.family,
                g1W = W.g$g1W, pDelta1 = cbind(A0 = rep(1,n), A1 = rep(1,n)),
                V.Q = V.Q, Q.SL.library = L$Q.SL.library,
                Q.discreteSL = TRUE)$Qinit$Q
    colnames(W.Q) <- c("W.Q0", "W.Q1")			
  } else {
    W.Q <- NULL
  }
  
  # Evaluate conditional sampling probabilities
  if (is.null(pi)) {
    validCondSetNames <- c("A", "W", "Y")
    if (!(all(condSetNames %in% validCondSetNames))) {
      stop("condSetNames must be any combination of 'A', 'W', 'Y'")
    }
    if (any(condSetNames == "Y")) {
      if (any(is.na(Y))) {
        stop("Cannot condition on the outcome to evaluate sampling probabilities when some outcome values are missing")
      }
    }
    
    if (is.null(W.Q)) {
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames))
      colnames(d.pi)[-1] <- .getColNames(condSetNames, Wnames = colnames(W))
    } else {
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames), W.Q)
      colnames(d.pi)[-c(1, tail(1:ncol(d.pi), n = 2L))] <- .getColNames(condSetNames, Wnames = colnames(W))
    }
    res.subcal <- tmle:::estimateG(
      d = d.pi, g1W = NULL, gform = piform, SL.library = pi.SL.library,
      id = id, V = V.pi, message = "sampling weights", outcome = "Delta.W", discreteSL = pi.discreteSL,
      obsWeights = rep(1, nrow(W)), verbose = verbose
    )
    names(res.subcal)[1] <- "pi"
  } else {
    res.subcal <- list()
    res.subcal$pi <- pi
    res.subcal$type <- "User supplied values"
    res.subcal$coef <- NA
    res.subcal$discreteSL <- NULL
  }
  
  # Now bound normalized obsWeights after rescaling so they sum to sum(Delta.W)
  obsWeights <- Delta.W / res.subcal$pi
  obsWeights <- obsWeights / sum(obsWeights) * sum(Delta.W)
  ub <- sqrt(sum(Delta.W)) * log(sum(Delta.W)) / 5
  obsWeights <- tmle:::.bound(Delta.W/res.subcal$pi, c(0, ub))
  
  # call TMLE on fully-observed data
  argList <- list(...)
  argList$Y <- Y[Delta.W == 1]
  argList$A <- A[Delta.W == 1] 
  if (augmentW) {
    argList$W <- cbind(W[Delta.W == 1, , drop = FALSE], W.Q[Delta.W == 1, ], W.sub)
  } else {
    argList$W <- cbind(W[Delta.W == 1, , drop = FALSE], W.sub)
  }
  argList$Z <- Z[Delta.W == 1]
  argList$Delta <- Delta[Delta.W == 1]
  argList$obsWeights <- obsWeights[Delta.W == 1]
  argList$id <- id[Delta.W == 1]
  argList$family <- Q.family
  argList$verbose <- verbose
  if (rareOutcome){
    argList$Q.SL.library <- c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2") 
    argList$Q.discreteSL <- TRUE
    argList$V.Q <- 20
  }
  result <- list()
  result$tmle <- try(do.call(tmle, argList))
  if (inherits(result$tmle, "try-error")) {
    cat("Error calling tmle. Estimated sampling probabilites are being returned")
    result$tmle <- NULL
  }
  result$subcal <- res.subcal
  result$augW <- W.Q
  class(result) <- "subcalTMLE" # for tmleMSM can use same class
  return(result)
}


#-------------------------
# Function subcalTMLEMSM
# purpose: subset calibration TMLE for conditional effects
# Arguments:
#  	Y : outcome of interest (missingness allowed)
# 	A : binary treatment indicator
# 	W : matrix or data.frame of covariates measured on entire population
# 	Delta.W : Indicator of inclusion in subset with additional information
# 	W.sub : matrix or data.frame of covariates measured in subset population
# 			  (in same order as the corresponding rows in W. DO NOT include
# 			   rows for subjects not included in the subset)
# 	Z : Only required when the goal is estimating a controlled direct effect
# 		mediated by a binary variable,Z,  between treatment and outcome
# 	Delta : binary indicator that outcome Y is observed
# 	pi : optional vector of sampling probabilities
# 	piform : optional parametric regression model for estimating pi
# 	pi.SL.library : optional SL library specification for estimating pi (ignored when piform or pi is provided)
# 	V.pi : optional number of cross-validation folds for super learning (ignored when piform or pi is provided)
# 	pi.discreteSL: flag to indicate whether to use ensemble or discrete super learning (ignored when piform or pi is provided)
# 	condSet : variables to condition on when estimating pi. Covariates in W are always included. Can optionally also
# 			   condition on A and/or Y.
# 	id : optional indicator of independent units of observation
# 	Q.family : "gaussian" or "binomial", used to evaluate the TMLE
#   augmentW : TRUE to augment W with predicted values when A = 0 and A = 1
#			   from fitting the outcome regression conditional on A and W on all observations 
# 	rareOutcome: set to TRUE to automatically set V.Q = 20, Q.discreteSL = TRUE, and
# 		Q.SL.library = glm, glmnet, bart
# 	verbose : set to TRUE to print out informative messages
# 	...	other arugments passed to the tmleMSM function (see tmleMSM documentation in the tmle package)
# Return value
# 	object of class "subcalTMLE" (a list):
# 	 subcal - a list of details on estimating the conditional sampling probabilities
# 		pi - predicted probabilities
# 		coef - coefficients in the parametric model or SL ensemble weights
# 		type - "user-supplied regression formula" or "Super Learner"
# 		"discreteSL" - FALSE when ensemble SL used, TRUE when best single algorithm selected
# 	tmle - results returned from the call to the tmleMSM function with
# 		   observation weights equal to 1/pi. (class "tmleMSM", see tmle package documentation)
#-------------------------

subcalTMLEMSM <- function(Y, A, W, V, Delta.W, W.sub, Delta = rep(1, length(Y)),
                          pi = NULL, piform = NULL,
                          pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet", 
                                            "tmle.SL.dbarts.k.5"), 
                          V.pi = 10, pi.discreteSL = TRUE, condSetNames = c("A", "W"), 
                          id = NULL, Q.family = "gaussian", augmentW = TRUE, 
                          rareOutcome = FALSE, verbose = FALSE, ...) {
  if (!requireNamespace("tmle", versionCheck = ">=2.0", quietly = FALSE)) {
    stop("Loading required tmle package (>-2.0) failed", call. = FALSE)
  }
  
  if (is.null(id)) {
    id <- 1:length(Y)
  }
  
  if (is.vector(W.sub)) {
    W.sub <- as.matrix(W.sub)
    colnames(W.sub) <- "W.sub"
  }
  
  if (augmentW) {
    L <- list(...)
    
    n.effective <- sum(Delta)
    if (Q.family == "binomial") {
      n.effective <- min(c(table(Y[Delta == 1]) * 5, n.effective))
    }
    V.Q <- setV(n.effective)
    d.g <- data.frame(Y = Y, mget(condSetNames))
    W.g <- tmle:::estimateG(d = d.g, gform = "A ~ .", SL.library = "SL.glm", 
                            outcome = "A", message = "treatment assigment probability",
                            g1W = NULL, V = V.pi, discreteSL = TRUE,
                            id = id, obsWeights = rep(1, nrow(W)), verbose = verbose)
    W.Q <- tmle(Y, A, W, Delta = Delta, id = id, family = Q.family, 
                g1W = W.g$g1W, pDelta1 = cbind(A0 = rep(1,n), A1 = rep(1,n)),
                V.Q = V.Q, Q.SL.library = L$Q.SL.library,
                Q.discreteSL = TRUE)$Qinit$Q
    colnames(W.Q) <- c("W.Q0", "W.Q1")
  } else {
    W.Q <- NULL
  }
  
  V <- as.matrix(V)
  
  # Evaluate conditional sampling probabilities
  if (is.null(pi)) {
    validCondSetNames <- c("A", "V", "W", "Y")
    if (!(all(condSetNames %in% validCondSetNames))) {
      stop("condSetNames must be any combination of 'A', 'V', 'W', 'Y'")
    }
    if (any(condSetNames == "Y")) {
      if (any(is.na(Y))) {
        stop("Cannot condition on the outcome to evaluate sampling probabilities when some outcome values are missing")
      }
    }
    
    if (is.null(W.Q)) {
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames))
      colnames(d.pi)[-1] <- .getColNames(condSetNames, Wnames = colnames(W), Vnames = colnames(V))
    } else {
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames), W.Q)
      colnames(d.pi)[-c(1, tail(1:ncol(d.pi), n = 2L))] <- .getColNames(condSetNames, Wnames = colnames(W), Vnames = colnames(V))
    }
    res.subcal <- tmle:::estimateG(
      d = d.pi, g1W = NULL, gform = piform, SL.library = pi.SL.library,
      id = id, V = V.pi, message = "sampling weights", outcome = "Delta.W", discreteSL = pi.discreteSL,
      obsWeights = rep(1, nrow(W)), verbose = verbose
    )
    names(res.subcal)[1] <- "pi"
  } else {
    res.subcal <- list()
    res.subcal$pi <- pi
    res.subcal$type <- "User supplied values"
    res.subcal$coef <- NA
    res.subcal$discreteSL <- NULL
  }
  
  # Now bound normalized obsWeights after rescaling so they sum to sum(Delta.W)
  obsWeights <- Delta.W/res.subcal$pi
  obsWeights <- obsWeights / sum(obsWeights) * sum(Delta.W)
  ub <- sqrt(sum(Delta.W)) * log(sum(Delta.W)) / 5
  obsWeights <- tmle:::.bound(Delta.W/res.subcal$pi, c(0, ub))
  
  argList <- list(...)
  argList$Y <- Y[Delta.W == 1]
  argList$A <- A[Delta.W == 1] 
  if (augmentW) {
    argList$W <- cbind(W[Delta.W == 1, , drop = FALSE], W.Q[Delta.W == 1, ], W.sub)
  } else {
    argList$W <- cbind(W[Delta.W == 1, , drop = FALSE], W.sub)
  }
  argList$V <- V[Delta.W == 1, , drop = FALSE]
  argList$Delta <- Delta[Delta.W == 1]
  argList$obsWeights <- obsWeights[Delta.W == 1]
  argList$id <- id[Delta.W == 1]
  argList$family <- Q.family
  argList$verbose <- verbose
  if (rareOutcome){
    argList$Q.SL.library <- c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2") 
    argList$Q.discreteSL <- TRUE
    argList$V.Q <- 20
  }
  
  result <- list()
  result$tmleMSM <- try(do.call(tmleMSM, argList))		
  
  if (inherits(result$tmleMSM, "try-error")) {
    cat("Error calling tmleMSM. Estimated sampling probabilites are being returned")
    result$tmleMSM <- NULL
  }
  result$subcal <- res.subcal
  result$augW <- W.Q
  class(result) <- "subcalTMLE" # for tmleMSM can use same class
  return(result)
}

# utility functions for IPCW-TMLE ----------------------------------------------
# @param data the dataset (assumes in order and has is.complete variable)
# @param g_learner_lib the SL library for estimating the (treatment assignment) propensity score
# @param miss_learner_lib the SL library for estimating the probability of missingness
# @param q_learner_lib the SL library for estimating the outcome regression
# @param K the number of cross-validation folds for the Super Learner
# @return the TMLE estimators and SEs for the risk difference and odds ratio
run_tmle <- function(data, g_learner_lib = c("SL.glm"), miss_learner_lib = c("SL.glm"),
                     q_learner_lib = c("SL.glm"), K = 5) {
  O <- data
  O_cc <- O[complete.cases(O), ]
  # estimate propensity score of being missing
  # remove Ws (partially observed) and is.complete (delta) from SL
  miss_exclude_vars <- c(grep("_obs", names(O)), grep("is.complete", names(O)))
  # need to rename Y for SL to work
  dat_for_g_miss <- O[, -miss_exclude_vars]
  names(dat_for_g_miss)[names(dat_for_g_miss) == "Y"] <- "real_Y"
  # if (length(miss_learner_lib) == 1 & !is.list(miss_learner_lib)) {
  #   fitter <- eval(parse(text = miss_learner_lib))
  #   est_g_miss <- fitter(Y = O$is.complete,
  #                        X = dat_for_g_miss,
  #                        obsWeights = rep(1, nrow(O)),
  #                        newX = dat_for_g_miss,
  #                        family = "binomial")
  #   g_miss <- est_g_miss$pred
  # } else {
  est_g_miss <- SuperLearner::SuperLearner(
    Y = O$is.complete,
    X = dat_for_g_miss,
    SL.library = miss_learner_lib,
    cvControl = list(V = K, stratifyCV = TRUE),
    family = "binomial"
  )
  g_miss <- predict(est_g_miss)$pred
  # }
  
  miss_weights <- 1 / g_miss
  # normalized weights
  miss_weights_2 <- miss_weights / sum(miss_weights) * nrow(data)
  miss_weights <- miss_weights_2
  miss_weights_cc <- miss_weights[complete.cases(O)]
  # estimate propensity score of being assigned treatment
  tx_exclude_vars <- c(grep("X", names(O)), grep("Y", names(O)), grep("is.complete", names(O)))
  est_g_tx <- SuperLearner::SuperLearner(
    Y = O_cc$X,
    X = O_cc[, -tx_exclude_vars], # remove X, Y, is.complete (delta)
    obsWeights = miss_weights_cc,
    SL.library = g_learner_lib,
    cvControl = list(V = K, stratifyCV = TRUE),
    family = "binomial"
  )
  # estimate outcome regression
  y_exclude_vars <- c(grep("Y", names(O)), grep("is.complete", names(O)))
  if (length(unique(O_cc$Y)) == 2) {
    fam <- "binomial"
    cv_control <- list(V = K, stratifyCV = TRUE)
  } else {
    fam <- "gaussian"
    cv_control <- list(V = K)
  }
  est_q <- SuperLearner::SuperLearner(
    Y = O_cc$Y,
    X = O_cc[, -y_exclude_vars], # remove Y, is.complete (delta)
    SL.library = q_learner_lib,
    obsWeights = miss_weights_cc,
    cvControl = cv_control,
    family = fam
  )
  # update using TMLE
  w_for_tmle <- O_cc[, -tx_exclude_vars]
  a_for_tmle <- O_cc$X
  y_for_tmle <- O_cc$Y
  data_for_tmle <- data.frame(Y = y_for_tmle, X = a_for_tmle, w_for_tmle)
  # fluctuate the initial model
  final_q <- fluctuate_q(
    data = data_for_tmle, est_g = est_g_tx, est_q = est_q,
    weights = miss_weights_cc
  )
  # get TMLE estimators
  # q_1_star <- vector("numeric", length = nrow(O))
  # q_0_star <- vector("numeric", length = nrow(O))
  # q_1_star[O$is.complete == 1] <- final_q$q_1_star
  # q_0_star[O$is.complete == 1] <- final_q$q_0_star
  est_ey1_tmle <- mean(miss_weights_2[O$is.complete == 1] * final_q$q_1_star)
  est_ey0_tmle <- mean(miss_weights_2[O$is.complete == 1] * final_q$q_0_star)
  # est_ey1_tmle <- mean(O$is.complete * miss_weights * q_1_star, na.rm = TRUE)
  # est_ey0_tmle <- mean(O$is.complete * miss_weights * q_0_star, na.rm = TRUE)
  est_rd_tmle <- est_ey1_tmle - est_ey0_tmle
  est_or_tmle <- (est_ey1_tmle / (1 - est_ey1_tmle)) /
    (est_ey0_tmle / (1 - est_ey0_tmle))
  est_rr_tmle <- est_ey1_tmle / est_ey0_tmle
  # get standard errors
  # pred_g_all <- rep(0.5, length = nrow(O))
  # pred_g_all[O$is.complete == 1] <- predict(est_g_tx)$pred
  # q_star <- vector("numeric", length = nrow(O))
  # q_star[O$is.complete == 1] <- final_q$q_star
  est_eif_rd <- eif_rd(
    data = data_for_tmle, pred_q = final_q$q_star,
    pred_q_1 = final_q$q_1_star, pred_q_0 = final_q$q_0_star,
    pred_g = predict(est_g_tx)$pred, miss_weights = miss_weights[O$is.complete == 1],
    delta = O$is.complete[O$is.complete == 1], mu1 = est_ey1_tmle,
    mu0 = est_ey0_tmle
  )
  se_rd_tmle <- sqrt(var(est_eif_rd) / nrow(data_for_tmle))
  est_eif_or <- eif_or(
    data = data_for_tmle, pred_q = final_q$q_star,
    pred_q_1 = final_q$q_1_star, pred_q_0 = final_q$q_0_star,
    pred_g = predict(est_g_tx)$pred, miss_weights = miss_weights[O$is.complete == 1],
    delta = O$is.complete[O$is.complete == 1], mu1 = est_ey1_tmle,
    mu0 = est_ey0_tmle
  )
  se_or_tmle <- sqrt(var(est_eif_or) / nrow(data_for_tmle))
  est_eif_rr <- eif_rr(
    data = data_for_tmle, pred_q = final_q$q_star,
    pred_q_1 = final_q$q_1_star, pred_q_0 = final_q$q_0_star,
    pred_g = predict(est_g_tx)$pred, miss_weights = miss_weights[O$is.complete == 1],
    delta = O$is.complete[O$is.complete == 1], mu1 = est_ey1_tmle,
    mu0 = est_ey0_tmle
  )
  se_rr_tmle <- sqrt(var(est_eif_rr) / nrow(data_for_tmle))
  
  return(list(
    "results" = data.frame(
      "procedure" = "TMLE",
      "estimand" = c("RD", "OR", "RR"),
      "est" = c(est_rd_tmle, est_or_tmle, est_rr_tmle),
      "SE" = c(se_rd_tmle, se_or_tmle, se_rr_tmle)
    ),
    "coefs" = NA
  ))
}
# fluctuation-based updates for Q
# @param data the fully observed data
# @param est_q the initial estimator of E(Y | X = x, W = w, Delta = 1)
# @param est_g the estimator of E(X = 1 | W = w, Delta = 1)
# @param weights the IPCW weights
# @param tol the tolerance for approximately equal to zero (convergence)
# @param maxiter the maximum number of iterations
# @return the final updated estimator of E(Y | A = a, W = w, Delta = 1)
fluctuate_q <- function(data = NULL, est_q = NULL, est_g = NULL,
                        weights = rep(1, nrow(data)), tol = 1e-6,
                        maxiter = 100) {
  # initial estimators of the propensity and outcome regression
  g_w <- predict(est_g)$pred
  q_w <- predict(est_q)$pred
  q_1_w <- predict(est_q, newdata = data.frame(X = 1, subset(data, select = !(names(data) %in% c("X", "Y")))))$pred
  q_0_w <- predict(est_q, newdata = data.frame(X = 0, subset(data, select = !(names(data) %in% c("X", "Y")))))$pred
  Q <- data.frame(QAW = q_w, Q1W = q_1_w, Q0W = q_0_w)
  # compute covariates for fluctuation
  h_w <- data$X / g_w - (1 - data$X) / (1 - g_w)
  h_1_w <- 1 / g_w
  h_0_w <- (-1) / (1 - g_w)
  # initial fluctuation
  epsilon_fit <- glm(data$Y ~ -1 + h_w + offset(qlogis(Q$QAW)),
                     weights = weights, family = "quasibinomial"
  )
  epsilon <- coef(epsilon_fit)
  Q_j <- as.data.frame(plogis(apply(Q, 2, qlogis) + cbind(epsilon * h_w, epsilon * h_1_w, epsilon * h_0_w)))
  iter <- 1
  while ((abs(max(epsilon)) > tol) & (iter < maxiter)) {
    epsilon_fit <- glm(data$Y ~ -1 + h_w + offset(qlogis(Q_j$QAW)),
                       weights = weights, family = "quasibinomial"
    )
    epsilon <- coef(epsilon_fit)
    Q_j <- as.data.frame(plogis(apply(Q_j, 2, qlogis) + cbind(epsilon * h_w, epsilon * h_1_w, epsilon * h_0_w)))
    iter <- iter + 1
  }
  return(list("q_star" = Q_j$QAW, "q_1_star" = Q_j$Q1W, "q_0_star" = Q_j$Q0W))
}

# @param full_data_eif the EIF on the fully-observed data
# @param delta the indicator of being fully observed
# @param pred_miss the probability of being missing
# @param eif_preds the predictions from a regression of the full-data EIF on the observed data
get_obs_eif <- function(full_data_eif, delta, pred_miss, eif_preds) {
  full_data_eif2 <- vector("numeric", length = length(delta))
  full_data_eif2[delta == 1] <- full_data_eif
  full_data_eif2[delta == 0] <- 0
  eif <- delta / pred_miss * full_data_eif2 - (delta / pred_miss - 1) * eif_preds
  return(eif)
}
# influence function for the TMLE estimator
# @param data the entire dataset
# @param pred_q the initial estimator of E(Y | X = x, W = w, Delta = 1)
# @param pred_q_1 the initial estimator of E(Y | X = 1, W = w, Delta = 1)
# @param pred_q_0 the initial estimator of E(Y | X = 0, W = w, Delta = 1)
# @param pred_g the estimator of E(X | W = w, Delta = 1)
# @param pred_miss the estimator of E(Delta | W = w, X = x)
# @param delta the indicator of being fully observed
# @param mu1 the estimate of the outcome mean under treatment
# @param mu0 the estimate of the outcome mean under no treatment
eif_rd <- function(data = NULL, pred_q = NULL,
                   pred_q_1 = NULL, pred_q_0 = NULL, pred_g = NULL,
                   miss_weights = NULL,
                   delta = rep(1, nrow(data)),
                   mu1 = 0, mu0 = 0) {
  # g_1 <- pred_g  / miss_weights
  # g_0 <- (1 - pred_g)  / miss_weights
  g_1 <- pred_g
  g_0 <- 1 - pred_g
  eif <- miss_weights * ((data$X / g_1 - (1 - data$X) / g_0) * delta * (data$Y - pred_q) +
                           pred_q_1 - pred_q_0 - (mu1 - mu0))
  return(eif)
}
# influence function for the TMLE estimator
# @param data the entire dataset
# @param pred_q the initial estimator of E(Y | X = x, W = w, Delta = 1)
# @param pred_q_1 the initial estimator of E(Y | X = 1, W = w, Delta = 1)
# @param pred_q_0 the initial estimator of E(Y | X = 0, W = w, Delta = 1)
# @param pred_g the estimator of E(X | W = w, Delta = 1)
# @param pred_miss the estimator of E(Delta | W = w, X = x)
# @param delta the indicator of being fully observed
# @param mu1 the estimate of the outcome mean under treatment
# @param mu0 the estimate of the outcome mean under no treatment
eif_rr <- function(data = NULL, pred_q = NULL,
                   pred_q_1 = NULL, pred_q_0 = NULL, pred_g = NULL,
                   miss_weights = NULL,
                   delta = rep(1, nrow(data)),
                   mu1 = 0, mu0 = 0) {
  g_1 <- pred_g
  g_0 <- 1 - pred_g
  eif <- miss_weights * ((1 / mu1) * (data$X / g_1 * delta * (data$Y - pred_q) + pred_q_1 - mu1) -
                           (1 / mu0) * ((1 - data$X) / g_0 * delta * (data$Y - pred_q) + pred_q_0 - mu0))
  return(eif)
}
# @param data the entire dataset
# @param pred_q the initial estimator of E(Y | X = x, W = w, Delta = 1)
# @param pred_q_1 the initial estimator of E(Y | X = 1, W = w, Delta = 1)
# @param pred_q_0 the initial estimator of E(Y | X = 0, W = w, Delta = 1)
# @param pred_g the estimator of E(X | W = w, Delta = 1)
# @param pred_miss the estimator of E(Delta | W = w, X = x)
# @param delta the indicator of being fully observed
# @param mu1 the estimate of the outcome mean under treatment
# @param mu0 the estimate of the outcome mean under no treatment
eif_or <- function(data = NULL, pred_q = NULL,
                   pred_q_1 = NULL, pred_q_0 = NULL, pred_g = NULL,
                   miss_weights = NULL,
                   delta = rep(1, nrow(data)),
                   mu1 = 0, mu0 = 0) {
  g_1 <- pred_g
  g_0 <- 1 - pred_g
  eif <- miss_weights * (1 / (mu1 * (1 - mu1)) * (data$X / g_1 * delta * (data$Y - pred_q) + pred_q_1) -
                           1 / (mu0 * (1 - mu0)) * ((1 - data$X) / g_0 * delta * (data$Y - pred_q) + pred_q_0))
  return(eif)
}
