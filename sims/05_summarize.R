# summarize the results of the base-case simulation in run_simple_outcome_simple_MAR.bat (xscenario 1, yscenario 1.1, mscenario 1.1), print out to a latex table 

# this file is meant to be run on its own (to get all LaTeX tables) or "source"d as part of .Rmd reports

# load required packages -------------------------------------------------------
library("dplyr")
library("tidyr")
library("stringr")
library("optparse") # for getting command-line arguments
library("knitr") # for LaTeX tables
library("kableExtra") # for LaTeX tables

# Y, X, missing scenarios and estimators ---------------------------------------
# a single scenario run in run_simple_outcome_simple_MAR.bat
xscenarios <- "1"
yscenarios <- "1.1"
mscenarios <- "1.1"
if (!exists("estimators")) {
  estimators <- c("cc_oracle", "cc_population", "cc_noW", "cc", "ipw", "gr", "mice", "xgb", "rf", "tmle_m", "tmle_mto")
}
if (!exists("ns")) {
  ns <- c(10000)
}
if (!exists("estimand")) {
  estimand <- "cOR"
}
if (!exists("nreps")) {
  nreps <- 2500
  max_nreps <- 2500
}
if (!exists("tmle_nreps")) {
  tmle_nreps <- 2500
}
if (!exists("seeds")) {
  seeds <- 1
}
if (!exists("rd_mult")) {
  rd_mult <- 1
}
if (nreps == 1000) {
  seeds <- 1
}
all_scenarios <- expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[1], yscenario = yscenarios[1], seed = seeds)
unique_scenarios <- all_scenarios
if (!exists("plasmode")) {
  plasmode <- FALSE
}
if (!exists("write_out_tables")) {
  write_out_tables <- TRUE
}
if (!exists("robust")) {
  robust <- TRUE
}
# set up directory structrue ---------------------------------------------------
main_dir <- getwd()
proj_root <- paste0(main_dir)
setwd(proj_root)
production_code_dir <- paste0(main_dir, "/")
production_output_dir <- paste0(proj_root, "/Results/")
if (!dir.exists(production_output_dir)) {
  dir.create(production_output_dir, recursive = TRUE)
}

code_subdir <- NULL
output_subdir <- NULL

code_dir <- paste0(production_code_dir, code_subdir)
source(paste0(code_dir, "00_utils.R"))
output_dir <- paste0(production_output_dir, output_subdir)
latex_output_dir <- paste0(output_dir, "/latex_files/")
if (!dir.exists(latex_output_dir)) {
  dir.create(latex_output_dir, recursive = FALSE)
}

# read in results --------------------------------------------------------------
# the number of digits to round true values to
suffix <- paste0(ifelse(isTRUE(plasmode), "_plasmode.rds", ".rds"))
by_vars <- c("estimand", "yscenario")
if (!isTRUE(plasmode)) {
  by_vars <- c(by_vars, "xscenario")
}
round_digits <- 3
truths_oracle <- readRDS(file = paste0(production_output_dir, "true_values_oracle", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

truths_pop <- readRDS(file = paste0(production_output_dir, "true_values_pop", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

if (!isTRUE(plasmode)) {
  truths_oracle <- truths_oracle %>% select(-mscenario)
  truths_pop <- truths_pop %>% select(-mscenario)
}

colnames(truths_pop)[which(colnames(truths_pop)=="truth")] <- "truth_pop"
colnames(truths_oracle)[which(colnames(truths_oracle)=="truth")] <- "truth_oracle"
truths <- merge(truths_pop, truths_oracle,
                by=by_vars) %>% 
  group_by(pick(!contains("truth"))) %>% 
  slice(1) %>% 
  ungroup()

if (!isTRUE(plasmode)) {
  truths$yscenario <- as.numeric(as.character(truths$yscenario))
  truths$xscenario <- as.numeric(as.character(truths$xscenario))
}

# for the continuous covariate Zs
truths_oracle_zs <- readRDS(file = paste0(production_output_dir, "true_values_oracle_zs", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

truths_pop_zs <- readRDS(file = paste0(production_output_dir, "true_values_pop_zs", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

if (!isTRUE(plasmode)) {
  truths_oracle_zs <- truths_oracle_zs %>% select(-mscenario)
  truths_pop_zs <- truths_pop_zs %>% select(-mscenario)
}

colnames(truths_pop_zs)[which(colnames(truths_pop_zs)=="truth")] <- "truth_pop"
colnames(truths_oracle_zs)[which(colnames(truths_oracle_zs)=="truth")] <- "truth_oracle"
truths_zs <- merge(truths_pop_zs, truths_oracle_zs,
                by=by_vars) %>% 
  group_by(pick(!contains("truth"))) %>% 
  slice(1) %>% 
  ungroup()

if (!isTRUE(plasmode)) {
  truths_zs$yscenario <- as.numeric(as.character(truths_zs$yscenario))
  truths_zs$xscenario <- as.numeric(as.character(truths_zs$xscenario))
}

# for the continuous covariate Ws
truths_oracle_ws <- readRDS(file = paste0(production_output_dir, "true_values_oracle_ws", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

truths_pop_ws <- readRDS(file = paste0(production_output_dir, "true_values_pop_ws", suffix)) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "RD" ~ "mRD",
    estimand == "logRR" ~ "mlogRR"
  ))

if (!isTRUE(plasmode)) {
  truths_oracle_ws <- truths_oracle_ws %>% select(-mscenario)
  truths_pop_ws <- truths_pop_ws %>% select(-mscenario)
}

colnames(truths_pop_ws)[which(colnames(truths_pop_ws)=="truth")] <- "truth_pop"
colnames(truths_oracle_ws)[which(colnames(truths_oracle_ws)=="truth")] <- "truth_oracle"
truths_ws <- merge(truths_pop_ws, truths_oracle_ws,
                   by=by_vars) %>% 
  group_by(pick(!contains("truth"))) %>% 
  slice(1) %>% 
  ungroup()

if (!isTRUE(plasmode)) {
  truths_ws$yscenario <- as.numeric(as.character(truths_ws$yscenario))
  truths_ws$xscenario <- as.numeric(as.character(truths_ws$xscenario))
}


# hash table relating procedure names to nice procedure names
all_est_procedures <- c("cc_oracle", "cc_population", "cc", "cc_noW", 
                        "IPW", "iptw", "GR_Vanilla", "GR_X", "RRZ",
                        "MICE", "XGB", "RF", 
                        "TMLE", "TMLE-M", "TMLE-MTO", "tmle_m", "tmle_mto",
                        "ipcw-tmle_m", "ipcw-tmle_mto",
                        "ipcw-a-tmle_m", "ipcw-a-tmle_mto", 
                        "r-ipcw-tmle_m", "r-ipcw-tmle_mto",
                        "r-ipcw-a-tmle_m", "r-ipcw-a-tmle_mto")
all_est_procedures_fct <- factor(all_est_procedures,
                                 levels = all_est_procedures,
                                 labels = c("Oracle model", "Population model", "Complete-case", "Confounded model", 
                                            "IPW", "IPTW", "Raking (vanilla)", "Raking (X only)", "RRZ",
                                            "MICE", "MI-XGB", "MI-RF", 
                                            "IPCW-TMLE-M", "IPCW-TMLE-M", "IPCW-TMLE-MTO", 
                                            "IPCW-TMLE-M", "IPCW-TMLE-MTO",
                                            "IPCW-TMLE-M", "IPCW-TMLE-MTO",
                                            "IPCW-a-TMLE-M", "IPCW-a-TMLE-MTO", 
                                            "r-IPCW-TMLE-M", "r-IPCW-TMLE-MTO",
                                            "r-IPCW-a-TMLE-M", "r-IPCW-a-TMLE-MTO"))
procedure_hash_table <- data.frame(
  "procedure" = all_est_procedures,
  "nice_procedure" = all_est_procedures_fct
)

all_output_list <- NULL
all_zs_list <- NULL
all_ws_list <- NULL 
for (i in 1:nrow(all_scenarios)) {
  this_scenario <- all_scenarios[i, ]
  if (isTRUE(plasmode)) {
    this_output_dir <- paste0(output_dir, "plasmode/", this_scenario$yscenario,
                              "_n", this_scenario$n, "/")
  } else {
    this_output_dir <- paste0(output_dir, "m", this_scenario$mscenario, 
                              "_y", this_scenario$yscenario, 
                              "_x", this_scenario$xscenario, 
                              "_s", this_scenario$seed, "/")  
  }
  all_output_files <- list.files(this_output_dir, pattern = paste0("output_m", this_scenario$mscenario,
                                                                   "_y", this_scenario$yscenario,
                                                                   "_x", this_scenario$xscenario,
                                                                   "_s", this_scenario$seed))
  all_coef_files <- list.files(this_output_dir, pattern = paste0("coefs_m", this_scenario$mscenario,
                                                                 "_y", this_scenario$yscenario,
                                                                 "_x", this_scenario$xscenario,
                                                                 "_s", this_scenario$seed))
  # both "gr" and "rak" refer to vanilla raking; only keep one
  if (any(grepl("gr", all_output_files)) & any(grepl("rak", all_output_files))) {
    all_output_files <- all_output_files[!grepl("rak", all_output_files)]
    all_coef_files <- all_coef_files[!grepl("rak", all_coef_files)]
  }
  # both "ipcw-tmle_m" and "tmle_m" refer to the same thing; only keep one
  if (any(grepl("ipcw-tmle_m", all_output_files)) & any(grepl("tmle_m", all_output_files) & !grepl("ipcw", all_output_files))) {
    all_output_files <- all_output_files[!(grepl("tmle", all_output_files) & !grepl("ipcw", all_output_files))]
    all_coef_files <- all_coef_files[!(grepl("tmle", all_coef_files) & !grepl("ipcw", all_coef_files))]
  }
  if (nreps != max_nreps) {
    all_output_files <- all_output_files[grepl("1000", all_output_files)]
    all_coef_files <- all_coef_files[grepl("1000", all_coef_files)]
  } 
  if (tmle_nreps != max_nreps) {
    all_output_files <- c(all_output_files[!grepl("1000", all_output_files)], all_output_files[grepl("1000_est_tmle_m", all_output_files)])
    all_coef_files <- c(all_coefs_files[!grepl("1000", all_coefs_files)], all_coefs_files[grepl("1000_est_tmle_m", all_coefs_files)])
  }
  # if there are two versions of a file, e.g., "cc" and "cc_nreps_2500", keep only "nreps_2500" version
  minus_nreps <- gsub("_nreps_[0-9]*\\.", ".", all_output_files)
  minus_nreps_coef <- gsub("_nreps_[0-9]*\\.", ".", all_coef_files)
  if (length(minus_nreps) > 0) {
    use_last <- !(all(grepl("nreps", all_output_files)) | plasmode)
    all_output_files <- all_output_files[!duplicated(minus_nreps, fromLast = use_last)]
    all_coef_files <- all_coef_files[!duplicated(minus_nreps_coef, fromLast = use_last)]  
  }
  
  this_output_list <- lapply(as.list(seq_len(length(all_output_files))), function(j) {
    read.csv(file = paste0(this_output_dir, all_output_files[j]), stringsAsFactors = FALSE) %>% 
      mutate(seed = this_scenario$seed, .before = "n") %>% 
      select(-matches("nreps"))
  })
  all_output_list <- c(all_output_list, this_output_list)
  
  this_zs_list <- lapply(as.list(seq_len(length(all_coef_files))), function(j) {
    read.csv(file = paste0(this_output_dir, all_coef_files[j]), stringsAsFactors = FALSE) %>% 
      mutate(seed = this_scenario$seed, .before = "n") %>% 
      select(-matches("nreps")) %>% 
      filter(grepl("Zs", Variable)) %>% 
      rename(est = Estimate)
  })
  all_zs_list <- c(all_zs_list, this_zs_list)
  
  this_ws_list <- lapply(as.list(seq_len(length(all_coef_files))), function(j) {
    read.csv(file = paste0(this_output_dir, all_coef_files[j]), stringsAsFactors = FALSE) %>% 
      mutate(seed = this_scenario$seed, .before = "n") %>% 
      select(-matches("nreps")) %>% 
      filter(grepl("Ws", Variable)) %>% 
      rename(est = Estimate)
  })
  all_ws_list <- c(all_ws_list, this_ws_list)
}

# compile output and transform, if necessary -----------------------------------
# output for treatment assignment
all_output <- do.call(rbind, all_output_list) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "logRR" ~ "mlogRR",
    estimand == "RD" ~ "mRD",
    estimand == "OR" ~ "mOR",
    estimand == "RR" ~ "mRR",
    estimand == "cOR" ~ "cOR"
    # estimand != "canonical_parameter" & estimand != "logOR" & estimand != "logRR" & estimand != "RD" ~ estimand
  )) 

# Need to take log of TMLE OR (marginal and conditional) and RR output
# Tried doing this in format above and it did not work.
# estimand_log <- c("mOR","mRR","cOR") # Use once scale of cOR is confirmed (or not)
estimand_log <- c("mOR","mRR")
all_output$est_new <- ifelse(all_output$estimand %in% estimand_log,log(all_output$est),all_output$est)
all_output$est <- all_output$est_new

all_output <- all_output %>%
  mutate(estimand = case_when(
    estimand == "cOR" ~ "clogOR", # Note that the initial output is on the log scale, but not labeled as such.
    estimand == "mRR" ~ "mlogRR",
    estimand == "mOR" ~ "mlogOR",
    estimand == "mlogRR" ~ "mlogRR",
    estimand == "mlogOR" ~ "mlogOR",
    estimand == "clogOR" ~ "clogOR",
    estimand == "mRD" ~ "mRD"
    ),
         # Now truth and estimates on same scale (log for OR and RR, natural for RD), so CI scale is now identity.
         ci_scale = "identity") %>% 
  # Adding on truths
  left_join(truths, by = by_vars) %>%
  left_join(procedure_hash_table, by = "procedure") %>% 
  mutate(across(ends_with("scenario"), as.character)) %>% 
  filter(procedure != "iptw", procedure != "RRZ") 

# output for Zs
all_output_zs <- do.call(rbind, all_zs_list) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "logRR" ~ "mlogRR",
    estimand == "RD" ~ "mRD",
    estimand == "OR" ~ "mOR",
    estimand == "RR" ~ "mRR",
    estimand == "cOR" ~ "cOR"
    # estimand != "canonical_parameter" & estimand != "logOR" & estimand != "logRR" & estimand != "RD" ~ estimand
  )) 

# Need to take log of TMLE OR (marginal and conditional) and RR output
# Tried doing this in format above and it did not work.
# estimand_log <- c("mOR","mRR","cOR") # Use once scale of cOR is confirmed (or not)
estimand_log <- c("mOR","mRR")
all_output_zs$est_new <- ifelse(all_output_zs$estimand %in% estimand_log,log(all_output_zs$est),all_output_zs$est)
all_output_zs$est <- all_output_zs$est_new

all_output_zs <- all_output_zs %>%
  mutate(estimand = case_when(
    estimand == "cOR" ~ "clogOR", # Note that the initial output is on the log scale, but not labeled as such.
    estimand == "mRR" ~ "mlogRR",
    estimand == "mOR" ~ "mlogOR",
    estimand == "mlogRR" ~ "mlogRR",
    estimand == "mlogOR" ~ "mlogOR",
    estimand == "clogOR" ~ "clogOR",
    estimand == "mRD" ~ "mRD"
  ),
  # Now truth and estimates on same scale (log for OR and RR, natural for RD), so CI scale is now identity.
  ci_scale = "identity") %>% 
  # Adding on truths
  left_join(truths_zs, by = by_vars) %>%
  left_join(procedure_hash_table, by = "procedure") %>% 
  mutate(across(ends_with("scenario"), as.character)) %>% 
  filter(procedure != "iptw", procedure != "RRZ") %>%
  select(-Variable)


# output for Ws
all_output_ws <- do.call(rbind, all_ws_list) %>% 
  mutate(estimand = case_when(
    estimand == "canonical_parameter" ~ "clogOR",
    estimand == "logOR" ~ "mlogOR",
    estimand == "logRR" ~ "mlogRR",
    estimand == "RD" ~ "mRD",
    estimand == "OR" ~ "mOR",
    estimand == "RR" ~ "mRR",
    estimand == "cOR" ~ "cOR"
    # estimand != "canonical_parameter" & estimand != "logOR" & estimand != "logRR" & estimand != "RD" ~ estimand
  )) 

# Need to take log of TMLE OR (marginal and conditional) and RR output
# Tried doing this in format above and it did not work.
# estimand_log <- c("mOR","mRR","cOR") # Use once scale of cOR is confirmed (or not)
estimand_log <- c("mOR","mRR")
all_output_ws$est_new <- ifelse(all_output_ws$estimand %in% estimand_log,log(all_output_ws$est),all_output_ws$est)
all_output_ws$est <- all_output_ws$est_new

all_output_ws <- all_output_ws %>%
  mutate(estimand = case_when(
    estimand == "cOR" ~ "clogOR", # Note that the initial output is on the log scale, but not labeled as such.
    estimand == "mRR" ~ "mlogRR",
    estimand == "mOR" ~ "mlogOR",
    estimand == "mlogRR" ~ "mlogRR",
    estimand == "mlogOR" ~ "mlogOR",
    estimand == "clogOR" ~ "clogOR",
    estimand == "mRD" ~ "mRD"
  ),
  # Now truth and estimates on same scale (log for OR and RR, natural for RD), so CI scale is now identity.
  ci_scale = "identity") %>% 
  # Adding on truths
  left_join(truths_ws, by = by_vars) %>%
  left_join(procedure_hash_table, by = "procedure") %>% 
  mutate(across(ends_with("scenario"), as.character)) %>% 
  filter(procedure != "iptw", procedure != "RRZ") %>%
  select(-Variable)

# compute performance metric: nominal (ASE) coverage ---------------------------
# bdw: switched on 20240304 to not change the true value for TMLE
all_output_fixedtruth <- all_output 

add_cover_oracle_init <- all_output_fixedtruth %>% 
  mutate(bias_oracle_init = est - truth_oracle,
         nominal_cover_oracle_init = cover(mu = truth_oracle, est = est, 
                                           SE = SE, scale = ci_scale))

add_cover_pop_init <- add_cover_oracle_init %>% 
  mutate(bias_pop_init = est - truth_pop,
         nominal_cover_pop_init = cover(mu = truth_pop, est = est, SE = SE, 
                                        scale = ci_scale))

all_output_fixedtruth_zs <- all_output_zs 

add_cover_oracle_init_zs <- all_output_fixedtruth_zs %>% 
  mutate(bias_oracle_init = est - truth_oracle,
         nominal_cover_oracle_init = cover(mu = truth_oracle, est = est, 
                                           SE = SE, scale = ci_scale))

add_cover_pop_init_zs <- add_cover_oracle_init_zs %>% 
  mutate(bias_pop_init = est - truth_pop,
         nominal_cover_pop_init = cover(mu = truth_pop, est = est, SE = SE, 
                                        scale = ci_scale))

all_output_fixedtruth_ws <- all_output_ws 

add_cover_oracle_init_ws <- all_output_fixedtruth_ws %>% 
  mutate(bias_oracle_init = est - truth_oracle,
         nominal_cover_oracle_init = cover(mu = truth_oracle, est = est, 
                                           SE = SE, scale = ci_scale))

add_cover_pop_init_ws <- add_cover_oracle_init_ws %>% 
  mutate(bias_pop_init = est - truth_pop,
         nominal_cover_pop_init = cover(mu = truth_pop, est = est, SE = SE, 
                                        scale = ci_scale))

# note all NAs are from seed = 2, where we don't have ASEs

# compute performance metrics: oracle (ESE) coverage, MAD coverage -------------
empirical_ses <- add_cover_pop_init %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(ESE = sd(est, na.rm = TRUE), .groups = "drop")

mad <- add_cover_pop_init %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(MAD = mad(est, na.rm = TRUE), .groups = "drop")

oracle_coverage <- add_cover_pop_init %>% 
  left_join(empirical_ses,
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>% 
  left_join(mad, 
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>%
  mutate(oracle_cover_pop_init = cover(mu = truth_pop, est = est, SE = ESE,
                                       scale = ci_scale),
         oracle_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = ESE,
                                          scale = ci_scale),
         mad_cover_pop_init = cover(mu = truth_pop, est = est, SE = MAD, scale = ci_scale),
         mad_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = MAD, scale = ci_scale),
         power_init = power(mu = 0, est = est, SE = SE, scale = ci_scale),
         oracle_power_init = power(mu = 0, est = est, SE = ESE, scale = ci_scale),
         mad_power_init = power(mu = 0, est = est, SE = MAD, scale = ci_scale))

empirical_ses_zs <- add_cover_pop_init_zs %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(ESE = sd(est), .groups = "drop")

mad_zs <- add_cover_pop_init_zs %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(MAD = mad(est, na.rm = TRUE), .groups = "drop")

oracle_coverage_zs <- add_cover_pop_init_zs %>% 
  left_join(empirical_ses_zs,
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>% 
  left_join(mad_zs, 
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>%
  mutate(oracle_cover_pop_init = cover(mu = truth_pop, est = est, SE = ESE,
                                       scale = ci_scale),
         oracle_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = ESE,
                                          scale = ci_scale),
         mad_cover_pop_init = cover(mu = truth_pop, est = est, SE = MAD, scale = ci_scale),
         mad_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = MAD, scale = ci_scale),
         power_init = power(mu = 0, est = est, SE = SE, scale = ci_scale),
         oracle_power_init = power(mu = 0, est = est, SE = ESE, scale = ci_scale),
         mad_power_init = power(mu = 0, est = est, SE = MAD, scale = ci_scale))

empirical_ses_ws <- add_cover_pop_init_ws %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(ESE = sd(est), .groups = "drop")

mad_ws <- add_cover_pop_init_ws %>%
  filter(seed == 1) %>%
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>%
  summarize(MAD = mad(est, na.rm = TRUE), .groups = "drop")

oracle_coverage_ws <- add_cover_pop_init_ws %>% 
  left_join(empirical_ses_ws,
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>% 
  left_join(mad_ws, 
            by = c("yscenario", "mscenario", "xscenario", "nice_procedure", "estimand")) %>%
  mutate(oracle_cover_pop_init = cover(mu = truth_pop, est = est, SE = ESE,
                                       scale = ci_scale),
         oracle_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = ESE,
                                          scale = ci_scale),
         mad_cover_pop_init = cover(mu = truth_pop, est = est, SE = MAD, scale = ci_scale),
         mad_cover_oracle_init = cover(mu = truth_oracle, est = est, SE = MAD, scale = ci_scale),
         power_init = power(mu = 0, est = est, SE = SE, scale = ci_scale),
         oracle_power_init = power(mu = 0, est = est, SE = ESE, scale = ci_scale),
         mad_power_init = power(mu = 0, est = est, SE = MAD, scale = ci_scale))

# compute performance metric: prop. of sims that finished without error and returned reasonable results ---------------------------
reasonable_threshold <- log(10)
num_finished_and_reasonable <- add_cover_pop_init %>% 
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand, seed) %>% 
  summarize(prop_complete = n() / nreps * 100,
            prop_reasonable = sum(abs(est) < reasonable_threshold & abs(SE) < reasonable_threshold, na.rm = TRUE) / nreps * 100,
            .groups = "drop")

# check monte-carlo IDs that failed to converge (no output) or have large point estimates; return
unreasonable_ids <- add_cover_pop_init %>% 
  mutate(reasonable = as.numeric(abs(est) < reasonable_threshold & abs(SE) < reasonable_threshold)) %>% 
  filter(reasonable == 0) %>% 
  select(yscenario, mscenario, xscenario, seed, procedure, mc_id)

all_ids <- 1:nreps
unique_estimators <- unique(add_cover_pop_init$procedure)
all_id_mat <- do.call(rbind, lapply(as.list(1:length(unique_estimators)), function(i) {
  tibble::tibble(
    mc_id = rep(all_ids, nrow(all_scenarios)),
    yscenario = rep(all_scenarios$yscenario, each = nreps),
    xscenario = rep(all_scenarios$xscenario, each = nreps),
    mscenario = rep(all_scenarios$mscenario, each = nreps),
    seed = rep(all_scenarios$seed, each = nreps),
    procedure = unique_estimators[i]
  )
})) %>% 
  filter(seed == 1, !(procedure == "XGB"))
completed <- add_cover_pop_init %>% 
  group_by(mc_id, yscenario, xscenario, mscenario, seed, procedure) %>% 
  slice(1)
uncomplete_ids <- all_id_mat %>% filter(seed == 1) %>% 
  anti_join(completed %>% filter(seed == 1), by = c("mc_id", "yscenario", "xscenario", "mscenario", "seed", "procedure"))
all_bad_ids <- bind_rows(unreasonable_ids, uncomplete_ids) %>% 
  arrange(yscenario, mscenario, xscenario, seed, procedure, mc_id)

# obtain performance summaries by averaging over MC id -------------------------
# summaries get to average over seed
# obtain estimate of monte-carlo variance by taking sd over seed
initial_summary <- get_summaries(oracle_coverage, num_finished_and_reasonable)

# Monte-Carlo variances
mc_vars <- initial_summary %>% 
  select(-seed) %>% 
  group_by(yscenario, mscenario, xscenario, nice_procedure, estimand) %>% 
  summarize(across(where(is.double),
                   .fns = ~ sd(.x, na.rm = TRUE)), .groups = "drop")

mc_vars_overall <- mc_vars %>% 
  group_by(estimand) %>% 
  summarize(across(where(is.double), 
                   .fns = ~ max(.x, na.rm = TRUE)), .groups = "drop")

# save the summary over all estimators, estimands
plasmode_txt <- ifelse(plasmode, "_plasmode", "")
saveRDS(initial_summary, file = paste0(output_dir, "summaries_all", plasmode_txt, ".rds"))
saveRDS(mc_vars, file = paste0(output_dir, "mc_vars", plasmode_txt, ".rds"))
saveRDS(mc_vars_overall, file = paste0(output_dir, "mc_vars_overall", plasmode_txt, ".rds"))
write.csv(all_bad_ids, file = paste0(output_dir, "problematic_ids", plasmode_txt, ".csv"))

# for Zs
initial_summary_zs <- get_summaries(oracle_coverage_zs, num_finished_and_reasonable)

# save the summary over all estimators, estimands
saveRDS(initial_summary_zs, file = paste0(output_dir, "summaries_all_zs", plasmode_txt, ".rds"))

# for Ws
initial_summary_ws <- get_summaries(oracle_coverage_ws, num_finished_and_reasonable)

# save the summary over all estimators, estimands
saveRDS(initial_summary_ws, file = paste0(output_dir, "summaries_all_ws", plasmode_txt, ".rds"))

# create LaTeX tables ---------------------------------------------------------- 
summary_tib <- get_summary_tib(initial_summary, rd_mult = rd_mult,
                               reasonable_threshold = reasonable_threshold, threshold_scale = "log")

if (!plasmode) {
  # internal tests
  test_tib <- initial_summary %>% 
    group_by(yscenario, mscenario, xscenario, estimand) %>% 
    filter(nice_procedure != "Confounded model") %>% 
    summarize(max_census_bias = max(mean_bias_pop_truth),
              max_oracle_bias = max(mean_bias_oracle_truth),
              max_se_diff = max(ASE - ESE), .groups = "drop")
  
  # check GR vs MICE
  gr_vs_mice <- initial_summary %>% 
    filter(nice_procedure %in% c("Raking (vanilla)", "MICE"), seed == 1) %>% 
    pivot_wider(names_from = nice_procedure,
                values_from = mean_bias_pop_truth:rmse_oracle) %>% 
    mutate(census_bias_diff = abs(`mean_bias_pop_truth_Raking (vanilla)`) - abs(`mean_bias_pop_truth_MICE`),
           oracle_bias_diff = abs(`mean_bias_oracle_truth_Raking (vanilla)`) - abs(`mean_bias_oracle_truth_MICE`),
           ese_diff = `ESE_Raking (vanilla)` - `ESE_MICE`,
           ase_diff = `ASE_Raking (vanilla)` - `ASE_MICE`,
           census_ase_cover_diff = `nominal_coverage_pop_Raking (vanilla)` - `nominal_coverage_pop_MICE`,
           census_ese_cover_diff = `oracle_coverage_pop_Raking (vanilla)` - `oracle_coverage_pop_MICE`,
           oracle_ase_cover_diff = `nominal_coverage_oracle_Raking (vanilla)` - `nominal_coverage_oracle_MICE`,
           oracle_ese_cover_diff = `oracle_coverage_oracle_Raking (vanilla)` - `oracle_coverage_oracle_MICE`) 
  gr_vs_mice %>% 
    group_by(estimand) %>% 
    summarize(max_census_bias_diff = max(census_bias_diff),
              max_oracle_bias_diff = max(oracle_bias_diff),
              max_ese_diff = max(ese_diff),
              max_ase_diff = max(ase_diff))
}


summary_tib_zs <- get_summary_tib(initial_summary_zs, rd_mult = rd_mult,
                                  reasonable_threshold = reasonable_threshold, threshold_scale = "log")

summary_tib_ws <- get_summary_tib(initial_summary_ws, rd_mult = rd_mult,
                                  reasonable_threshold = reasonable_threshold, threshold_scale = "log")

# write out to LaTeX tables ----------------------------------------------------
this_estimand <- mc_vars_overall %>% 
  mutate(estimand = gsub("log", "", estimand)) %>% 
  filter(estimand == !!estimand)
max_continuous_mce <- max(this_estimand[!grepl("coverage", names(this_estimand))][-1])
max_binary_mce <- max(this_estimand[grepl("coverage", names(this_estimand))])
# tables broken down by scenario and estimand
all_scenarios_estimands <- unique_scenarios %>% 
  mutate(cOR = "cOR", mOR = "mOR", mRD = "mRD", mRR = "mRR") %>% 
  pivot_longer(cOR:mRR, names_to = "estimand") %>% 
  # mutate(clogOR = "clogOR", mlogOR = "mlogOR", mRD = "mRD", mlogRR = "mlogRR") %>% 
  # pivot_longer(clogOR:mlogRR, names_to = "estimand") %>% 
  select(-value)
if (all(is.na(summary_tib$`M scenario`))) {
  all_scenarios_estimands$xscenario <- ""
  all_scenarios_estimands$mscenario <- ""
  unique_scenarios$xscenario <- ""
  unique_scenarios$mscenario <- ""
  summary_tib$`M scenario` <- ""
  summary_tib$`X scenario` <- ""
  summary_tib_zs$`M scenario` <- ""
  summary_tib_zs$`X scenario` <- ""
  summary_tib_ws$`M scenario` <- ""
  summary_tib_ws$`X scenario` <- ""
}
robust_txt <- ifelse(robust, "_robust", "")
if (write_out_tables) {
  for (i in 1:nrow(all_scenarios_estimands)) {
    # separate tables for oracle and population truth (keeping oracle and nominal coverage)
    this_scenario <- all_scenarios_estimands[i, ]
    this_summary_tib <- summary_tib %>% 
      filter(`X scenario` == this_scenario$xscenario,
             `Y scenario` == this_scenario$yscenario,
             `M scenario` == this_scenario$mscenario,
             `Estimand` == this_scenario$estimand) %>% 
      select(-`Y scenario`, -`X scenario`, -`M scenario`,
             -Estimand)
    oracle_census_equal <- (this_scenario$yscenario == 1 | this_scenario$yscenario == 1.1 | this_scenario$yscenario == 1.15)
    # table for oracle truth
    oracle_kable <- create_report_kable(this_summary_tib, oracle_census_equal = oracle_census_equal,
                                        estimand = estimand, oracle = TRUE, matched = FALSE,
                                        this_scenario = this_scenario, ns = ns, nreps = nreps,
                                        tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                        max_continuous_mce = max_continuous_mce,
                                        true_val = this_summary_tib$truth_oracle[1], 
                                        robust = robust)
    oracle_kable %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_oracletruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
    
    # table for population truth
    pop_kable <- create_report_kable(this_summary_tib, oracle_census_equal = oracle_census_equal,
                                     estimand = estimand, oracle = FALSE, matched = FALSE,
                                     this_scenario = this_scenario, ns = ns, nreps = nreps,
                                     tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                     max_continuous_mce = max_continuous_mce,
                                     true_val = this_summary_tib$truth_pop[1], 
                                     robust = robust)
    pop_kable %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_poptruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
    
    # table for oracle truth, Zs
    oracle_kable_zs <- create_report_kable(summary_tib_zs, oracle_census_equal = oracle_census_equal,
                                           estimand = estimand, oracle = TRUE, matched = FALSE,
                                           this_scenario = this_scenario, ns = ns, nreps = nreps,
                                           tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                           max_continuous_mce = max_continuous_mce,
                                           true_val = summary_tib_zs$truth_oracle[1], 
                                           robust = robust)
    oracle_kable_zs %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_zs_oracletruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
    
    # table for population truth, Zs
    pop_kable_zs <- create_report_kable(summary_tib_zs, oracle_census_equal = oracle_census_equal,
                                     estimand = estimand, oracle = FALSE, matched = FALSE,
                                     this_scenario = this_scenario, ns = ns, nreps = nreps,
                                     tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                     max_continuous_mce = max_continuous_mce,
                                     true_val = summary_tib_zs$truth_pop[1], 
                                     robust = robust)
    pop_kable %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_zs_poptruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
    # table for oracle truth, Ws
    oracle_kable_ws <- create_report_kable(summary_tib_ws, oracle_census_equal = oracle_census_equal,
                                           estimand = estimand, oracle = TRUE, matched = FALSE,
                                           this_scenario = this_scenario, ns = ns, nreps = nreps,
                                           tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                           max_continuous_mce = max_continuous_mce,
                                           true_val = summary_tib_ws$truth_oracle[1], 
                                           robust = robust)
    oracle_kable_ws %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_ws_oracletruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
    
    # table for population truth, Ws
    pop_kable_ws <- create_report_kable(summary_tib_ws, oracle_census_equal = oracle_census_equal,
                                        estimand = estimand, oracle = FALSE, matched = FALSE,
                                        this_scenario = this_scenario, ns = ns, nreps = nreps,
                                        tmle_nreps = tmle_nreps, max_binary_mce = max_binary_mce,
                                        max_continuous_mce = max_continuous_mce,
                                        true_val = summary_tib_ws$truth_pop[1], 
                                        robust = robust)
    pop_kable %>% 
      kableExtra::save_kable(file = paste0(output_dir, "/latex_files/",
                                           "summary_table_ws_poptruth_", this_scenario$estimand,
                                           robust_txt,
                                           "_m", this_scenario$mscenario,
                                           "_y", this_scenario$yscenario,
                                           "_x", this_scenario$xscenario, ".tex"))
  }
}
