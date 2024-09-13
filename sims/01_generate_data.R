# Data generation

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 20230725       bdw       harmonize dataset creation into a single file
# 20230807       bdw       add missing data scenarios, match name to overleaf doc
# 20231114       chloe     add y scenarios 3.1--3.4
# 20231115       bdw       add y scenarios 4, 4.1
# 20231115       chloe     modified scenarios 3.1--3.4 to have outcome rates more comparable to other Y scenarios.
# 20231227       bdw       mild outcome and missingness misspecification
# 20240326       bdw       finalize missing data scenarios 2.2, 2.4: misspecified based on tree/interaction structure, with 40 or 80% missing data, missingness depends on Y
# 20240404       chloe     added y scenarios 1.15 and 4.15
# 20240424       bdw       finalize missing-data scenarios 2.1, 2.3; misspecified based on tree/interaction structure, 40 or 80% missing data, missingness does not depend on Y
# 20240501       eaj       Include missing data scenarios 2.5-2.9, MNAR
# 20240522       eaj       Adding X scenarios 1.5, 1.6 for higher correlation with unobserved variable
# 20240603       eaj       Adding M scenario 2.91

# generate an entire dataset ---------------------------------------------------
# @param N the sample size
# @param XScenario the covariate scenario
# @param YScenario the outcome regression model
# @param missScenario the missing-data mechanism
# @param lowcor low correlation
# @param midcor middle correlation
# @param highcor high correlation
# @param gencor general correlation
# @return a dataset, with covariates, outcome, missing-data indicator
gen_data <- function(N = 100, XScenario = 1, YScenario = 1, missScenario = 1, 
                     lowcor = 0.1, midcor = 0.4, highcor = 0.7, gencor = 0.2) {
  # generate covariates (everything fully observed)
  covariates <- XFunc(N = N, lowcor = lowcor, midcor = midcor, highcor = highcor, gencor = gencor, XScenario = XScenario)
  # generate outcome based on covariates (everything fully observed)
  Y <- YFunc(data = covariates, YScenario = YScenario)$Y
  df <- cbind.data.frame(Y = Y, covariates)
  # generate missing-data indicators
  miss <- missFunc(data = df, missScenario = missScenario)
  # create missing-data versions of confounders
  df$Ws_obs <- ifelse(miss == 1, NA, df$Ws)
  df$Ww_obs <- ifelse(miss == 1, NA, df$Ww)
  df$is.complete <- as.numeric(miss == 0)
  # return
  # return(list("df" = df, "obs_vars" = c("Y", "X", "Zs", "Zw", "As", "Aw", "Ws_obs", "Ww_obs", "is.complete")))
  return(df)
}
# outcome regression model -----------------------------------------------------
# @param data the covariate data (expects columns X, Zw, Zs, Ww, Ws, Uw, Us)
# @param YScenario the outcome regression model scenario of interest
# @return a vector of outcomes
YFunc <- function(data, YScenario = 1) {
  # specify regression model parameters
  # scenario 1, 1.1: simple outcome model
  if (YScenario == 1) {	# no treatment effect
    beta0 <- -2.4   
    betaX <- log(1.0)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZw2 <- 0
    betaZs2 <- 0
    betaWw2 <- 0
    betaWs2 <- 0
  }
  if (YScenario == 1.1) {	# conditional log OR of treatment = log(1.5)
    beta0 <- -2.4   
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZw2 <- 0
    betaZs2 <- 0
    betaWw2 <- 0
    betaWs2 <- 0
  }	
  if (YScenario == 1.15) {
    beta0 <- -5.1
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZw2 <- 0
    betaZs2 <- 0
    betaWw2 <- 0
    betaWs2 <- 0
  }
  if (YScenario == 1.16) {
    beta0 <- -3.9
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZw2 <- 0
    betaZs2 <- 0
    betaWw2 <- 0
    betaWs2 <- 0
  }
  if (YScenario == 1.17) {
    beta0 <- -3.4
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZw2 <- 0
    betaZs2 <- 0
    betaWw2 <- 0
    betaWs2 <- 0
  }
  # scenario 2: mild outcome model misspecification
  # achieve this in several ways:
  # 2, 2.1: exactly the same as 4, 4.1 except that the nonlinearities appear further in the tails of the distribution of X
  # so these scenarios are captured below
  # scenario 2.2 (Pam's 2.4), 2.3 (Pam's 2.5), 2.4 (Pam's 2.51), 2.5 (Pam's 2.52) from Pam
  if (YScenario == 2.2) {
    beta0 <- -2.4   
    betaX <- log(1.5)
    betaZw <- log(1.5)
    delta <- 0.5
    zeta <- 1.5
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
  }
  if (YScenario == 2.3) {
    beta0 <- -2.4  
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- 0
    delta <- 2.1
    zeta <- log(1.8)
    betaWw <- 0
    betaWs <- 0
    betaUw <- 0
    betaUs <- 0
  }
  if (YScenario == 2.4) {
    beta0 <- -2.4  
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- 0
    delta <- 4
    zeta <- log(1.2)
    betaWw <- 0
    betaWs <- 0.54
    betaUw <- 0
    betaUs <- 0
  }
  if (YScenario == 2.5) {
    beta0 <- -2.4  
    betaX <- log(1.5)
    betaZw <- log(1.2)
    betaZs <- 0
    delta <- 4
    zeta <- log(1.2)
    betaWw <- 0
    betaWs <- 0.54
    betaUw <- 0
    betaUs <- 0
  }
  # scenario 3: mild functional form misspecification (smooth functions of Z and/or W)
  if (YScenario == 3.1){ # Functional form misspecification of Zs, null X
    # beta0 <- -2.4   
    beta0 <- -10
    betaX <- 0 
    betaZw <- 0
    betaZs <- 0
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZs2 <- log(5.0)
    betaZw2 <- log(7.0)
    betaWs2 <- 0
    betaWw2 <- 0
  }
  if (YScenario == 3.2){ # Functional form misspecification of Zs, non-null X
    # beta0 <- -2.4   
    beta0 <- -10
    betaX <- log(0.5)
    betaZw <- 0
    betaZs <- 0
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUw <- 0
    betaUs <- 0
    betaZs2 <- log(5.0)
    betaZw2 <- log(7.0)
    betaWs2 <- 0
    betaWw2 <- 0
  }
  if (YScenario == 3.3){ # Functional form misspecification of Ws, null X
    # beta0 <- -2.4   
    beta0 <- -10
    betaX <- 0
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- 0
    betaWs <- 0
    betaUw <- 0
    betaUs <- 0
    betaZs2 <- 0
    betaZw2 <- 0
    betaWs2 <- log(5.0)
    betaWw2 <- log(7.0)
  }
  if (YScenario == 3.4){
    # beta0 <- -2.4   
    beta0 <- -10
    betaX <- log(0.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- 0
    betaWs <- 0
    betaUw <- 0
    betaUs <- 0
    betaZs2 <- 0
    betaZw2 <- 0
    betaWs2 <- log(5.0)
    betaWw2 <- log(7.0)
  }
  # scenario 4: complex functional form misspecification (including binary/categorical variables, interactions)
  if (YScenario %in% c(2, 4)) {
    beta0 <- -3
    betaX <- 0
    betaZw <- 0
    betaZs <- -0.4
    betaWw <- -0.6
    betaWs <- 0.5
    betaUw <- 0
    betaUs <- 0
    betaZw_low <- 0.1
    betaZw_high <- 0.8
    betaWwWs <- 1
    betaWsZs <- 3
    betaWsZw <- 1
  }
  if (YScenario %in% c(2.1, 4.1)) {
    beta0 <- -3   
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- -0.4
    betaWw <- -0.6
    betaWs <- 0.5
    betaUw <- 0
    betaUs <- 0
    betaZw_low <- 0.1
    betaZw_high <- 0.8
    betaWwWs <- 1
    betaWsZs <- 3
    betaWsZw <- 1
  }
  # incorrect outcome model, depends on Us
  if (YScenario == 2.11) {
    beta0 <- -2.5   
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUs <- -log(1.75)
    betaUw <- 0
  }
  # incorrect outcome model, depends on Us, rare
  if (YScenario == 2.17) {
    beta0 <- -3.56  
    betaX <- log(1.5)
    betaZw <- log(1.5)
    betaZs <- -log(1.3)
    betaWw <- log(1.5)
    betaWs <- -log(1.75)
    betaUs <- -log(1.75)
    betaUw <- 0
  }
  if (YScenario == 4.15) {
    beta0 <- -6.45
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- -0.4
    betaWw <- -0.6
    betaWs <- 0.5
    betaUw <- 0
    betaUs <- 0
    betaZw_low <- 0.1
    betaZw_high <- 0.8
    betaWwWs <- 1
    betaWsZs <- 3
    betaWsZw <- 1
  }
  if (YScenario == 4.16) {
    beta0 <- -4.9
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- -0.4
    betaWw <- -0.6
    betaWs <- 0.5
    betaUw <- 0
    betaUs <- 0
    betaZw_low <- 0.1
    betaZw_high <- 0.8
    betaWwWs <- 1
    betaWsZs <- 3
    betaWsZw <- 1
  }
  if (YScenario == 4.17) {
    beta0 <- -4.1
    betaX <- log(1.5)
    betaZw <- 0
    betaZs <- -0.4
    betaWw <- -0.6
    betaWs <- 0.5
    betaUw <- 0
    betaUs <- 0
    betaZw_low <- 0.1
    betaZw_high <- 0.8
    betaWwWs <- 1
    betaWsZs <- 3
    betaWsZw <- 1
  }
  if (YScenario %in% c(1, 1.1, 1.15, 1.16, 1.17, 2.11, 2.17)) {
    Xbeta <- beta0 + betaX * data$X + betaZw * data$Zw + betaZs * data$Zs + 
      betaWw * data$Ww + betaWs * data$Ws
  } else if (YScenario %in% c(2, 2.1)) {
    # note these are the same as 4, 4.1 below but are further out in the tails of the X distribution
    Zw_low <- as.numeric(data$Zw < -3)
    Zw_high <- as.numeric(data$Zw > 3)
    Zs_low <- as.numeric(data$Zs < -3)
    Xbeta <- beta0 + betaX * data$X + betaWw * data$Ww + betaZw_low * Zw_low + betaZw_high * Zw_high + 
      betaWwWs * data$Ww * data$Ws + betaZs * Zs_low + betaWs * data$Ws + betaWsZs * data$Ws * Zs_low +
      betaWsZw * data$Ws * Zw_high + betaUw * data$Uw + betaUs * data$Us
  } else if (YScenario %in% c(2.2, 2.3, 2.4, 2.5)) {
    ws_extra <- delta * I(abs(data$Ws) > zeta) * (abs(data$Ws) - zeta)
    if (YScenario %in% c(2.2, 2.5)) {
      ws_extra <- ws_extra * sign(data$Ws)
    }
    Xbeta <- beta0 + betaX * data$X + betaZw * data$Zw + 
      betaZs * data$Zs + betaWw * data$Ww + betaWs * data$Ws + 
      ws_extra 
  } else if (YScenario %in% c(3.1, 3.2, 3.3, 3.4)) {
    Xbeta <- beta0 + betaX * data$X + betaZw * data$Zw + betaZs* data$Zs + betaWw * data$Ww + betaWs * data$Ws + betaUw * data$Uw + betaUs * data$Us +
      betaZw2 * exp(data$Zw) + betaZs2 * (data$Zs^2) + betaWw2 * exp(data$Ww) + betaWs2 * (data$Ws^2)
  } else if (YScenario %in% c(4, 4.1, 4.15, 4.16, 4.17)) {
    Zw_low <- as.numeric(data$Zw < -0.5)
    Zw_high <- as.numeric(data$Zw > 2)
    Zs_low <- as.numeric(data$Zs < -1)
    Xbeta <- beta0 + betaX * data$X + betaWw * data$Ww + betaZw_low * Zw_low + betaZw_high * Zw_high + 
      betaWwWs * data$Ww * data$Ws + betaZs * Zs_low + betaWs * data$Ws + betaWsZs * data$Ws * Zs_low +
      betaWsZw * data$Ws * Zw_high + betaUw * data$Uw + betaUs * data$Us
  } else {
    # not yet implemented
  }
  Y <- as.numeric(runif(nrow(data)) < stats::plogis(Xbeta))
  return(list("Y" = Y, "beta" = c(betaX, betaZs, betaWs)))
} #End Yfunc

# missing data model -----------------------------------------------------------
# @param data the covariate data (expects columns Y, X, Zw, Zs, Ww, Ws, Uw, Us)
# @param missScenario the missing-data scenario of interest
# @return a vector of missing data indicators
missFunc<-function(data, missScenario = 1, pMiss = 0, return_prob = FALSE) {
  ### Generates a missing inidicator vector: assumes the set of confounders (W) are either complete or missing
  if (missScenario == 1)  { ### Scenario 1 MAR, pMiss about 40% missingness 
    alpha0 <- -.67   
    alphaX <- log(1.5)
    alphaZw <- log(2.5)
    alphaZs <- log(2)
    alphaUs <- 0
    alphaY <- 0
  }
  
  if (missScenario == 1.1) {
    alpha0 <- -.67   
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- 0
    alphaY <- log(2.5)
  }
  
  # if (missScenario == 2.1)  { ### Scenario 2 mispec MAR, pMiss about .4, unobserved covar U 
  #   alpha0 <- -.78   
  #   alphaX <-  log(1.5)
  #   alphaZw <- log(2.5)
  #   alphaZs <- log(2)
  #   alphaUs <- log(3)
  #   alphaY <-  log(2)
  # }
  # if (missScenario == 2.3)  { ### Scenario 3 mispec MAR, pMiss about .4, unobserved covar U 
  #   alpha0 <- -.78   
  #   alphaX <-  log(1.5)
  #   alphaZw <- log(2.5)
  #   alphaZs <- log(2)
  #   alphaUs <- log(3)
  #   alphaY <- log(2)
  # }
  # misspecified MAR, missingness does not depend on Y, 40 or 80% missing data
  if (missScenario %in% c(2.1, 2.3)) { 
    alpha0 <- -1.5
    if (missScenario == 2.3) {
      alpha0 <- alpha0 + 4.3
    }
    alphaX <- 1
    alphaZw <- 0
    alphaZs <- -1
    alphaUs <- 0
    alphaZw_1 <- -3
    alphaZw_2 <- -2.5
    alphaZw_3 <- 2
    alphaZw_4 <- 3
    alphaZsZw <- 1
    alphaZsX <- 1
    alphaZwX <- 1
    alphaY <- 0
    alphaZsY <- 0
    alphaXY <- 0
  }
  # misspecified MAR, missingness depends on Y, 40 or 80% missing data
  if (missScenario %in% c(2.2, 2.4)) { 
    alpha0 <- -0.9
    if (missScenario == 2.4) {
      alpha0 <- alpha0 + 2.2
    }
    alphaX <- 1
    alphaZw <- 0
    alphaZs <- -2
    alphaUs <- 0
    alphaZw_high <- 2
    alphaZw_low <- -0.9
    alphaZsZw <- 0
    alphaZsX <- 3
    alphaZwX <- 0
    alphaY <- 0.2
    alphaZsY <- -3
    alphaXY <- 0
  }
  if (missScenario == 2.5) { #MNAR, unobserved covariate, 40% missing, same as 1.1 otherwise
    alpha0 <- -.97   
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- log(2.5)
    alphaUw <- 0
    alphaY <- log(2.5)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.6) { #MNAR, depends directly on missing covariates, 40% missing, same as 1.1 otherwise 
    alpha0 <- -.67   
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- 0
    alphaUw <- 0
    alphaY <- log(2.5)
    alphaWs <- log(2.5) #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- log(2.5) #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.7) { #MNAR, unobserved covariate, 80% missing, same as 1.1 otherwise 
    # alpha0 <- -2.2 
    alpha0 <- 1.28 # -.67 + 1.95   
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- log(2.5)
    alphaUw <- 0
    alphaY <- log(2.5)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.8) { #MNAR, depends directly on missing covariates, 80% missing, same as 1.1 otherwise
    alpha0 <- -.67 +2.3  
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- 0
    alphaUw <- 0
    alphaY <- log(2.5)
    alphaWs <- log(2.5) #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- log(2.5) #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.9) { #MNAR, unobserved covariate, 40% missing, same as 1.1 otherwise, alphaUs effect increased, others decreased
    alpha0 <- -.67   
    alphaX <- log(0.5)
    alphaZw <- log(0.5)
    alphaZs <- log(0.5)
    alphaUs <- log(3.5)
    alphaUw <- log(1)
    alphaY <- log(2.5)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.91) { #MNAR, unobserved covariate, 40% missing, same as 1.1 otherwise, alphaUs and alphaUw high, others decreased
    alpha0 <- -.67   
    alphaX <- log(0.5)
    alphaZw <- log(0.5)
    alphaZs <- log(0.5)
    alphaUs <- log(2.5)
    alphaUw <- log(2.5)
    alphaY <- log(2.5)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.92) { #MNAR, unobserved covariate, 40% missing, same as 1.1 otherwise, alphaUs and alphaUw high, others decreased
    alpha0 <- -.67   
    alphaX <- log(1)
    alphaZw <- log(1)
    alphaZs <- log(1)
    alphaUs <- log(2.5)
    alphaUw <- log(2.5)
    alphaY <- log(1)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 2.93) { #MNAR, unobserved covariate, 40% missing, same as 1.1 otherwise, alphaUs and alphaUw high, others decreased
    alpha0 <- -.67   
    alphaX <- log(1)
    alphaZw <- log(1)
    alphaZs <- log(1)
    alphaUs <- log(2.5)
    alphaUw <- log(2.5)
    alphaY <- log(2.5)
    alphaWs <- 0 #Ws, Ww have mean 0, but this goes through a logit function.  May impact overall prevalence of missingness.
    alphaWw <- 0 #Higher values of Ws+Ww are more likely to be missing.  Imputations will be lower than they should be.
  }
  if (missScenario == 3)  { ### Missing scenario 1, but with 80% missing instead. 
    alpha0 <- -.67+2.3   
    alphaX <- log(1.5)
    alphaZw <- log(2.5)
    alphaZs <- log(2)
    alphaUs <- 0
    alphaY <- 0
  }
  if (missScenario == 3.1) { # Scenario 1.1, but with 80% missing instead.
    alpha0 <- -.67+1.75
    alphaX <- log(2.5)
    alphaZw <- log(1.5)
    alphaZs <- log(1.5)
    alphaUs <- 0
    alphaY <- log(2.5)
  }
  
  if (missScenario %in% c(1, 1.1, 3, 3.1)) {
    XbetaM <- alpha0 + alphaX * data$X + alphaZw * data$Zw + alphaZs * data$Zs + 
      alphaUs * data$Us + alphaY * data$Y
  } else if (missScenario %in% c(2.5, 2.6, 2.7, 2.8, 2.9, 2.91, 2.92, 2.93)){ #MNAR value-based or unobserved covariate
    XbetaM <- alpha0 + alphaX * data$X + alphaZw * data$Zw + alphaZs * data$Zs + 
      alphaUs * data$Us + alphaUw * data$Uw + alphaY * data$Y + alphaWs * data$Ws + alphaWw * data$Ww
  } else if (missScenario %in% c(2.2, 2.4)) {
    Zw_low <- as.numeric(data$Zw < -0.5)
    Zw_high <- as.numeric(data$Zw > 1)
    Zs_low <- as.numeric(data$Zs < -1)
    XbetaM <- alpha0 + alphaX * data$X + alphaY * data$Y +
      alphaZw_low * Zw_low + alphaZw_high * Zw_high +
      alphaZsZw * Zs_low * Zw_low + alphaZs * Zs_low +
      alphaZsY * data$Y * Zs_low + alphaZsX * data$X * Zs_low +
      alphaXY * data$X * data$Y + alphaZwX * data$X * Zw_high
  } else {
    threshold_1 <- qnorm(.2)
    threshold_2 <- qnorm(.4)
    threshold_3 <- qnorm(.6)
    threshold_4 <- qnorm(.8)
    Zw_1 <- as.numeric(data$Zw < threshold_1)
    Zw_2 <- as.numeric(threshold_1 >= data$Zw & data$Zw < threshold_2)
    Zw_3 <- as.numeric(threshold_3 < data$Zw & data$Zw <= threshold_4)
    Zw_4 <- as.numeric(data$Zw > threshold_4)
    Zs_low <- as.numeric(data$Zs < -1)
    XbetaM <- alpha0 + alphaX * data$X + alphaY * data$Y +
      alphaZw_1 * Zw_1 + alphaZw_2 * Zw_2 + alphaZw_3 * Zw_3 + alphaZw_4 * Zw_4 +
      alphaZsZw * Zs_low * Zw_1 + alphaZs * Zs_low +
      alphaZsY * data$Y * Zs_low + alphaZsX * data$X * Zs_low +
      alphaXY * data$X * data$Y + alphaZwX * data$X * Zw_4
  }
  miss <- as.numeric(runif(nrow(data)) < stats::plogis(XbetaM))
  if (return_prob) {
    return(stats::plogis(XbetaM))
  } else {
    return(miss)
  }
} ### End missFunc

# generate de novo covariates --------------------------------------------------
# @param N the sample size
# @param lowcor low correlation
# @param midcor middle correlation
# @param highcor high correlation
# @param gencor general correlation
# @param XScenario the scenario for generating covariates
# @return a dataset with the covariates
XFunc <- function(N = 100, lowcor = 0.1, midcor = 0.4, 
                  highcor = 0.7, gencor = 0.2, XScenario = 1) {
  if (XScenario %in% c(1, 1.5, 1.6)) {
    corMat <- matrix(gencor,nrow=9,ncol=9)
    corMat[2,1] <- corMat[1,2] <- midcor  ###cor (X,Zs)
    corMat[3,1] <- corMat[1,3] <- lowcor  ###cor (X,Zw)
    corMat[4,1] <- corMat[1,4] <- midcor  ###cor (X,Ws)
    corMat[5,1] <- corMat[1,5] <- lowcor  ###cor (X,Ww)
    corMat[6,1] <- corMat[1,6] <- midcor  ###cor (X,Us)
    corMat[7,1] <- corMat[1,7] <- lowcor  ###cor (X,Uw)
    corMat[8,1] <- corMat[1,8] <- midcor  ###cor (X,As)
    corMat[9,1] <- corMat[1,9] <- lowcor  ###cor (X,Aw)
    
    if (XScenario == 1.5) { #moderate correlation between unobserved variable and missing variable 
      corMat[6,4] <- corMat[4,6] <- midcor  ###cor (Ws,Us)
      corMat[7,5] <- corMat[5,7] <- midcor  ###cor (Ww,Uw)
    }
    if (XScenario == 1.6) { #high correlation between unobserved variable and missing variable 
      corMat[6,4] <- corMat[4,6] <- highcor  ###cor (Ws,Us)
      corMat[7,5] <- corMat[5,7] <- highcor  ###cor (Ww,Uw)
    }
    varMat <- corMat
    diag(varMat) <- 1
    
    XZWUA <- mvtnorm::rmvnorm(N, rep(0,9), varMat)
    
    latentX <- XZWUA[,1]
    Zs <- XZWUA[,2]
    Zw <- XZWUA[,3]
    Ws <- XZWUA[,4]
    Ww <- XZWUA[,5]
    Us <- XZWUA[,6]
    Uw <- XZWUA[,7]
    As <- XZWUA[,8]
    Aw <- XZWUA[,9]
    
    ### Treatment covariate of interest. 
    ### Not 100% if this is the best treatment assignment mechanism for scenario 1.
    ### Open to suggestions. 
    ### Maybe better to try a direct X propensity to treat model?
    pX <- 0.4
    X <- ifelse(latentX < quantile(latentX, p = pX), 1, 0)
  } # End XScenario1
  if (XScenario == 1.1) {
    corMat <- matrix(gencor,nrow=9,ncol=9)
    corMat[2,1] <- corMat[1,2] <- midcor  ###cor (X,Zs)
    corMat[3,1] <- corMat[1,3] <- lowcor  ###cor (X,Zw)
    corMat[4,1] <- corMat[1,4] <- midcor  ###cor (X,Ws)
    corMat[5,1] <- corMat[1,5] <- lowcor  ###cor (X,Ww)
    corMat[6,1] <- corMat[1,6] <- midcor  ###cor (X,Us)
    corMat[7,1] <- corMat[1,7] <- lowcor  ###cor (X,Uw)
    corMat[8,1] <- corMat[1,8] <- midcor  ###cor (X,As)
    corMat[9,1] <- corMat[1,9] <- lowcor  ###cor (X,Aw)
    
    corMat[4,8] <- corMat[8,4] <- 0.80  ###cor (Ws,As) #establishes strongly correlated auxillary variables for Ws, Ww, the variables with missing data
    corMat[5,8] <- corMat[8,5] <- 0.80  ###cor (Ww,Aw)
    
    varMat <- corMat
    diag(varMat) <- 1
    
    XZWUA <- mvtnorm::rmvnorm(N, rep(0,9), varMat)
    
    latentX <- XZWUA[,1]
    Zs <- XZWUA[,2]
    Zw <- XZWUA[,3]
    Ws <- XZWUA[,4]
    Ww <- XZWUA[,5]
    Us <- XZWUA[,6]
    Uw <- XZWUA[,7]
    As <- XZWUA[,8]
    Aw <- XZWUA[,9]
    
    pX <- 0.4
    X <- ifelse(latentX < quantile(latentX, p = pX), 1, 0)
  } # End XScenario1.1
  data <- data.frame(X = X, Zs = Zs, Zw = Zw, Ws = Ws, Ww = Ww, Us = Us, Uw = Uw, As = As, Aw = Aw)
  return(data)
}
# obtain plasmode covariates ---------------------------------------------------
# coming soon!
