# Main analysis file
cat("Starting simulation\n")

# load required packages -------------------------------------------------------
library("mvtnorm") # for generating data
library("pROC")
library("survey") # for IPW, calibration
library("sandwich") # for robust SEs for glms
library("mice") # for MI
library("SuperLearner") # for TMLE
library("glmnet") # for TMLE
library("xgboost") # for TMLE
library("ranger") # for TMLE
library("dbarts") # for TMLE
library("optparse") # for getting command-line arguments
library("future.apply") # for parallelization
library("dplyr") # for piping, mutate
library("readr")
library("tmle") # version 2.0.0
library("mixgb")
library("marginaleffects")

main_dir <- getwd()
proj_root <- paste0(main_dir)
setwd(proj_root)
production_code_dir <- paste0(main_dir, "/")
production_output_dir <- paste0(proj_root, "/Results/")
if (!dir.exists(production_output_dir)) {
  dir.create(production_output_dir, recursive = TRUE)
}

# get arguments from the command line ------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--nreps-total", default = 2500, help = "Number of simulation replicates in total")
parser <- add_option(parser, "--n", default = 100, help = "Sample size")
parser <- add_option(parser, "--yscenario", default = 1, help = "Outcome regression model")
parser <- add_option(parser, "--mscenario", default = 1, help = "Missing data model")
parser <- add_option(parser, "--xscenario", default = 1, help = "Covariates and treatment assignment model")
parser <- add_option(parser, "--estimator", default = "cc_oracle", help = "The estimator(s) to run")
parser <- add_option(parser, "--use-cached-datasets", default = 0, help = "Should we use cached datasets or create a new dataset each time?")
parser <- add_option(parser, "--data-only", default = 0, help = "Should we only generate datasets?")
parser <- add_option(parser, "--seed", default = 1, help = "Change the seed (already depends on y, m, x scenarios)")
parser <- add_option(parser, "--testing", default = 0, help = "Are we testing?")
parser <- add_option(parser, "--subdir", default = "", help = "Subdirectory to use if testing")
parser <- add_option(parser, "--plasmode", default = 0, help = "Is this a plasmode simulation?")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
# make argument for directory: main results vs testing results

code_subdir <- NULL
data_subdir <- NULL
output_subdir <- NULL
nreps_total <- args$nreps_total
n <- args$n
yscenario <- args$yscenario
mscenario <- args$mscenario
xscenario <- args$xscenario
estimator <- args$estimator
use_cached_datasets <- args$use_cached_datasets
data_only <- args$data_only
seed <- args$seed
plasmode <- args$plasmode

# check if binary outcome or not
if (yscenario %in% c("1", "1.1")) {
  fam <- "binomial"
} else {
  fam <- "binomial"
}
# if plasmode and yscenario = SH_90day or SH_365day, then set rare_outcome = TRUE, otherwise FALSE
if (plasmode & (grepl("SH_90day", yscenario) | grepl("SH_365day", yscenario))) {
  rare_outcome <- TRUE
} else {
  rare_outcome <- FALSE
}
# print out the args
cat("\nArguments: \n",
    "Number of replicates: ", nreps_total, "\n",
    "Sample size: ", n, "\n",
    "Outcome regression scenario: ", yscenario, "\n",
    "Missing-data scenario: ", mscenario, "\n",
    "Covariates scenario: ", xscenario, "\n",
    "Estimator: ", estimator, "\n",
    ifelse(use_cached_datasets == 1, "Using cached datasets", "Generating new datasets"), "\n",
    "Plasmode:", plasmode)
code_dir <- paste0(production_code_dir, code_subdir)
if (plasmode) {
  data_dir <- paste0(proj_root, "Data/plasmode data sets/")
  output_dir <- paste0(production_output_dir, output_subdir, "plasmode/")
  this_model <- ifelse(grepl("glm_", yscenario), "glm", "tree")
  this_outcome <- gsub("glm_", "", gsub("tree_", "", yscenario))
  this_tx <- "numEpiType"
} else {
  data_dir <- paste0(production_output_dir, data_subdir)
  output_dir <- paste0(production_output_dir, output_subdir)
  this_model <- NULL
  this_outcome <- NULL
  this_tx <- "X"
}

source(paste0(code_dir, "00_utils.R"))
source(paste0(code_dir, "01_generate_data.R"))
source(paste0(code_dir, "02_methods_design.R"))
source(paste0(code_dir, "02_methods_mi.R"))
source(paste0(code_dir, "02_methods_xgb.R"))
source(paste0(code_dir, "02_methods_tmle.R"))
source(paste0(code_dir, "03_estimate.R"))

est_vec <- unlist(strsplit(estimator, ";", fixed = TRUE))
est_for_save <- paste0(est_vec, collapse = "_")

# set up directory to save results to and filename
if (plasmode) {
  this_output_dir <- paste0(output_dir, yscenario, "_n", n, "/")
  this_data_dir <- paste0(data_dir, yscenario, "/", "N", n, "/")
  filename_prefix <- paste0("plasmode_data")
  filename_suffix <- ".Rdata"
  read_func <- function(file) {
    load(file)
    # 'file' contains a dataset named 'data'
    data
  }
} else {
  this_output_dir <- paste0(output_dir, "m", mscenario, "_y", yscenario, 
                            "_x", xscenario, "_s", seed, "/")
  this_data_dir <- paste0(data_dir, "m", mscenario, "_y", yscenario, 
                          "_x", xscenario, "_s", seed, "/datasets/")    
  filename_prefix <- paste0("data_m", mscenario,"_y", yscenario,
                            "_x", xscenario, "_n", n, "_id")
  filename_suffix <- ".rds"
  read_func <- readRDS
}

if (!dir.exists(this_data_dir)) {
  dir.create(this_data_dir, recursive = TRUE)
}
if (!dir.exists(this_output_dir)) {
  dir.create(this_output_dir, recursive = TRUE)
}
# run the simulation -----------------------------------------------------------
# read in cached datasets, if we're using these
cat("Loading datasets\n")
if (use_cached_datasets & !plasmode) {
  cached_datasets <- lapply(
    as.list(1:nreps_total), function(i) {
      read_func(file = paste0(
        this_data_dir, filename_prefix, i, filename_suffix
      ))
    }
  )
} else if (plasmode) {
  cached_datasets <- TRUE
  cached_dataset <- read_func(paste0(this_data_dir, filename_prefix, 1, filename_suffix))
} else {
  cached_datasets <- NULL
}
# set up Super Learner library (only used if "TMLE" is an estimator)
xgb_tune_params <- list(max_depth = c(1, 3), shrinkage = c(1e-2, 1e-1), ntrees = 500)
num_n <- as.numeric(n)
if (is.na(num_n)) {
  num_n <- nrow(cached_dataset)
}
min_node_sizes <- round(seq.int(num_n / 100, num_n / 10, length.out = 4))
rf_tune_params <- list(num.trees = 500, min.node.size = min_node_sizes,
                       verbose = FALSE)
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params, 
                               detailed_names = TRUE, name_prefix = "xgb")
rf_learners <- create.Learner("SL.ranger", tune = rf_tune_params, 
                              detailed_names = TRUE, name_prefix = "rf")
continuous_learner_lib <- c("SL.glm", xgb_learners$names, rf_learners$names)
binary_learner_lib <- c("SL.glm", xgb_learners$names, rf_learners$names)
simple_lib <- "SL.glm"
tmle_args <- list("g_lib" = binary_learner_lib, 
                  "miss_lib" = binary_learner_lib, 
                  "q_lib" = continuous_learner_lib, "K" = 10,
                  "phase1_covars" = c("Y", "X", "Zs", "Zw", "As", "Aw"),
                  "phase2_covars" = c("Ws_obs", "Ww_obs"))
if (plasmode) {
  tmle_args$phase1_covars <- c("Y", "numEpiType", "sex", "AgeAtIndex",
                               "AgeIndexsq", "AgeAtIndex_cat", "anx_dx_priorYr",
                               "aud_dx_priorYr", "Charlson_cat", "priorSH",
                               "MHIP_prior5yr")
  tmle_args$phase2_covars <- c("IndexTrtPHQ8_score_cat_obs", "IndexTrtPHQ_item9_score_cat_obs")
}
n_imp <- 25
n_imp_raking <- 10
mice_args <- list("n_imp" = n_imp, "maxit" = 20)
raking_args <- list("NimpRaking" = n_imp_raking)

# get the current seed
if (plasmode) {
  ns <- c("5000", "15000", "30000", "full")
  models <- c("glm", "tree")
  outcomes <- c("SH_90day", "SH_365day", "SH_HOSP_365day", "SH_HOSP_1826day")
  current_seed <- round(which(this_outcome == outcomes) * 10) + round(which(this_model == models) * 100) + 
    round(which(n == ns) * 1000)
  num_cores <- parallel::detectCores() - 6
} else {
  current_seed <- round(yscenario * 10) + round(mscenario * 100) + round(xscenario * 100) + 
    round(n * 1000) + round((seed - 1) * 51)
  num_cores <- parallel::detectCores()
}

# set up parallelization (platform-agnostic)
future::plan(multisession)
options(future.globals.maxSize = +Inf)
set.seed(current_seed)
seeds <- future_lapply(as.list(seq_len(nreps_total)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
# run the simulation
cat("\nRunning simulation\n")
start <- Sys.time()
output_list <- future.apply::future_lapply(
  X = as.list(seq_len(nreps_total)), FUN = function(i) {
    tryCatch(
      investigate_performance_once(
        mc_id = i, n = n, XScenario = xscenario, fam = fam,
        YScenario = yscenario, missScenario = mscenario,
        lowcor = 0.2, midcor = 0.4, highcor = 0.7, gencor = 0.2,
        outcome_name = "Y", tx_name = this_tx, missing_indicator = "is.complete",
        estimators = est_vec, tmle_args = tmle_args, mi_args = mice_args,
        raking_args = raking_args, cached_datasets = cached_datasets,
        data_only = data_only, plasmode = plasmode,
        data_dir = this_data_dir, filename_prefix = filename_prefix,
        filename_suffix = filename_suffix, read_func = read_func, 
        rare_outcome = rare_outcome
      ), error = function(e) {
        message(conditionMessage(e))
        print(paste0("Error occurred when running Monte-Carlo iteration ", i))
        traceback()
        list("results" = NULL, "coefs" = NULL)
      }
    )
  }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
  future.conditions = "message",
  future.globals = structure(TRUE, add = c(xgb_learners$names, rf_learners$names))
)
end <- Sys.time()
cat("\nElapsed time: ", format(end - start), "\n")
# write results ----------------------------------------------------------------
if (data_only) {
  lapply(as.list(1:length(output_list)), function(i) {
    saveRDS(output_list[[i]], file = paste0(
      this_data_dir, "/data_m", mscenario,"_y", yscenario,
      "_x", xscenario, "_n", n, "_id", i, ".rds" 
    ))
  })
} else {
  output <- do.call(rbind, lapply(output_list, function(x) x$results))
  output <- output %>% 
    mutate(yscenario = yscenario, mscenario = mscenario, xscenario = xscenario, seed = seed, n = n,
           .before = "procedure") %>% 
    select(yscenario, mscenario, xscenario, seed, n, procedure, mc_id, estimand, est, SE)
  coefs <- do.call(rbind, lapply(output_list, function(x) {
    if (!is.null(x$coefs)) {
      x$coefs %>%
        mutate(yscenario = yscenario, mscenario = mscenario, xscenario = xscenario, seed = seed, n = n,
              .before = "Variable")
    }
  }))
  # save each estimator separately (so that we can update one at a time if need be)
  for (i in 1:length(est_vec)) {
    this_output <- output %>%
      filter(procedure == get_nice_procedure(est_vec[i]))
    these_coefs <- coefs %>%
      filter(procedure == get_nice_procedure(est_vec[i]))
    # save as csvs
    readr::write_csv(this_output, file = paste0(
      this_output_dir, "/output_m", mscenario,"_y", yscenario,
      "_x", xscenario, "_s", seed, "_n", n, "_est_", est_vec[i], "_nreps_", nreps_total, ".csv"
    ))
    readr::write_csv(these_coefs, file = paste0(
      this_output_dir, "/coefs_m", mscenario,"_y", yscenario,
      "_x", xscenario, "_s", seed, "_n", n, "_est_", est_vec[i], "_nreps_", nreps_total, ".csv"
    ))
  }
}
cat("\nSimulation complete!\n")
