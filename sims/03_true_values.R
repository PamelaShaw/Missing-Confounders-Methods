# NOTES
# Chloe 11/20/2023: updated to compute log true values and both oracle and population-level true values.

# compute true values (marginal, conditional) for each combination of outcome regression and covariates
library("tidyr")
library("dplyr")
library("future.apply")
library("ranger")
# These are the limiting values that the estimation procedures are targeting
# For parametric approaches: the canonical parameter from the regression analysis (e.g., log odds ratio)
# For all approaches: marginal parameters: risk difference, relative risk, odds ratio
this_os <- .Platform$OS.type
if (grepl("windows", this_os, ignore.case = TRUE)) {
  main_dir <- "G:/"
} else {
  main_dir <- ""
}
proj_root <- paste0(main_dir, "CTRHS/Sentinel/Y14_Task_Orders_2022/Subset Calibration/")
production_code_dir <- paste0(proj_root, "R code/")
production_output_dir <- paste0(proj_root, "Results/")

all_production_files <- gsub(".R", "", list.files(production_code_dir, pattern = ".R"))
get_latest_version <- function(all_files = NULL, file = "00_utils") {
  these_files <- all_files[grepl(file, all_files)]
  dates_chr <- gsub(paste0(file, "_"), "", these_files)
  dates <- as.Date(dates_chr, format = "%Y%m%d")
  today <- Sys.Date()
  latest <- dates_chr[which.min(today - dates)]
  return(latest)
}

utils_version <- get_latest_version(all_files = all_production_files,
                                    file = "00_utils")
source(paste0(production_code_dir, "00_utils_", utils_version, ".R"))
data_version <- get_latest_version(all_files = all_production_files,
                                   file = "01_generate_data")
source(paste0(production_code_dir, "01_generate_data_", data_version, ".R"))

# compute true values for fully-synthetic scenarios ----------------------------

xscenarios <- c("1", "1.1", "1.5", "1.6")
yscenarios <- c("1", "1.1", "3.1", "3.2", "3.3", "3.4", "4", "4.1", "1.15", "1.16", "1.17", "4.15", "4.16", "4.17", "2.11", "2.17")
mscenarios <- c("1", "1.1", "2.2", "2.4", "3", "3.1", "2.5", "2.6", "2.7", "2.8", "2.9", "2.1",  "2.3")
# the list of "all_scenarios" depends on which aim this is
aim_1_scenarios <- rbind(
  # Scenarios where both the outcome and missing-data models are correctly specified
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[1:2], yscenario = yscenarios[1:2]),
  # from 04 Jan 2024: suppress X scenario 1.1 for now
  # expand.grid(xscenario = xscenarios[2], mscenario = mscenarios, yscenario = yscenarios[1:2]),
  # Scenarios where only the outcome model is misspecified
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[2], yscenario = yscenarios[3:8])
)
aim_2_scenarios <- rbind(
  # Scenarios where only the missing-data model is misspecified MAR
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[c(3:4, 12:13)], yscenario = yscenarios[2]),
  # Scenarios where both the outcome and missing-data model are misspecified MAR
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[c(3:4, 12:13)], yscenario = yscenarios[8]),
  # Scenarios with higher proportion of missing data
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[5:6], yscenario = yscenarios[c(2, 6, 8)]),
  # Scenarios with a rare outcome
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[c(2, 6)], yscenario = yscenarios[c(9:14)]),
  # Scenarios with MNAR missing-data model; MNAR-value
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[c(8, 10)], yscenario = yscenarios[c(2, 8, 11, 14)]),
  # Scenarios with MNAR missing-data model; MNAR-unobserved
  expand.grid(xscenario = xscenarios[1], mscenario = mscenarios[c(7, 9)], yscenario = yscenarios[c(2, 11, 15, 16)])
)
all_xym <- rbind(aim_1_scenarios, aim_2_scenarios)
# set lowcor, midcor, highcor, gencor to be the same values used to generate datasets in the simulation file
lowcor <- 0.2
midcor <- 0.4
highcor <- 0.7
gencor <- 0.2

# for how we did this before (not in parallel) see 03_true_values_20231211.R
future::plan(multisession)
set.seed(20230802)
output_list <- future.apply::future_lapply(
  X = as.list(seq_len(nrow(all_xym))), FUN = function(i) {
    get_true_values(all_xym = all_xym, iter = i, sample_size = 1e6,
                    nreps = 10, lowcor = lowcor, midcor = midcor, gencor = gencor,
                    highcor = highcor, y_only = FALSE)
  }, future.seed = TRUE
)
all_output <- do.call(rbind.data.frame, output_list)

# true values for treatment variable
true_values_oracle_long <- all_output %>% 
  filter(type == "oracle", variable == "X") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long <- all_output %>% 
  filter(type == "census", variable == "X") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long, file = paste0(production_output_dir, "true_values_oracle.rds"))
saveRDS(true_values_pop_long, file = paste0(production_output_dir, "true_values_pop.rds"))

# true values for continuous variable that is always observed
true_values_oracle_long_zs <- all_output %>% 
  filter(type == "oracle", variable == "Zs") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long_zs <- all_output %>% 
  filter(type == "census", variable == "Zs") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long_zs, file = paste0(production_output_dir, "true_values_oracle_zs.rds"))
saveRDS(true_values_pop_long_zs, file = paste0(production_output_dir, "true_values_pop_zs.rds"))

# true values for continuous variable that is sometimes missing
true_values_oracle_long_ws <- all_output %>% 
  filter(type == "oracle", variable == "Ws") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long_ws <- all_output %>% 
  filter(type == "census", variable == "Ws") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long_ws, file = paste0(production_output_dir, "true_values_oracle_ws.rds"))
saveRDS(true_values_pop_long_ws, file = paste0(production_output_dir, "true_values_pop_ws.rds"))

# check many replications at N = 10K
set.seed(20230802)
output_list_10K <- future.apply::future_lapply(
  X = as.list(seq_len(nrow(all_xym))), FUN = function(i) {
    get_true_values(all_xym = all_xym, iter = i, sample_size = 1e4,
                    nreps = 250, lowcor = lowcor, midcor = midcor, gencor = gencor,
                    highcor = highcor, y_only = FALSE)
  }, future.seed = TRUE
)
all_output_10K <- do.call(rbind.data.frame, output_list_10K)

# true values for treatment variable
true_values_oracle_long_10K <- all_output_10K %>% 
  filter(type == "oracle", variable == "X") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long_10K <- all_output_10K %>% 
  filter(type == "census", variable == "X") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long_10K, file = paste0(production_output_dir, "true_values_oracle_10K.rds"))
saveRDS(true_values_pop_long_10K, file = paste0(production_output_dir, "true_values_pop_10K.rds"))

# true values for continuous variable that is always observed
true_values_oracle_long_zs_10K <- all_output_10K %>% 
  filter(type == "oracle", variable == "Zs") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long_zs_10K <- all_output_10K %>% 
  filter(type == "census", variable == "Zs") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long_zs_10K, file = paste0(production_output_dir, "true_values_oracle_zs_10K.rds"))
saveRDS(true_values_pop_long_zs_10K, file = paste0(production_output_dir, "true_values_pop_zs_10K.rds"))

# true values for continuous variable that is sometimes missing
true_values_oracle_long_ws_10K <- all_output_10K %>% 
  filter(type == "oracle", variable == "Ws") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")
true_values_pop_long_ws_10K <- all_output_10K %>% 
  filter(type == "census", variable == "Ws") %>% 
  select(-type, -variable) %>% 
  group_by(xscenario, yscenario, mscenario, estimand) %>% 
  summarize(truth = mean(value), .groups = "drop")

saveRDS(true_values_oracle_long_ws_10K, file = paste0(production_output_dir, "true_values_oracle_ws_10K.rds"))
saveRDS(true_values_pop_long_ws_10K, file = paste0(production_output_dir, "true_values_pop_ws_10K.rds"))

# check against each other
true_values_oracle_long %>% 
  left_join(true_values_oracle_long_10K %>% 
              rename(truth_10K = truth), by = c("xscenario", "yscenario", "mscenario", "estimand")) %>% 
  mutate(diff = truth - truth_10K) %>% 
  pull(diff) %>% 
  summary()
true_values_pop_long %>% 
  left_join(true_values_pop_long_10K %>% 
              rename(truth_10K = truth), by = c("xscenario", "yscenario", "mscenario", "estimand")) %>% 
  mutate(diff = truth - truth_10K) %>% 
  pull(diff) %>% 
  summary()

# compute true values for plasmode scenarios -----------------------------------
yscenarios <- apply(expand.grid(model = c("glm", "tree"), outcome = c("_SH_HOSP_1826day", "_SH_90day", "_SH_365day")), 1, paste0, collapse = "")
xscenarios <- c("")
mscenarios <- c("")
ns <- c("full")
all_xym <- rbind(
  expand.grid(xscenario = xscenarios, mscenario = mscenarios, yscenario = yscenarios,
              n = ns)
)
plasmode_folder <- paste0(proj_root, "Data/plasmode data sets/")
nreps_total <- 10

future::plan(sequential)
set.seed(20240119)
seeds <- future_lapply(as.list(seq_len(nrow(all_xym))), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = 20240119)
output_list <- future.apply::future_lapply(
  X = as.list(seq_len(nrow(all_xym))), FUN = function(i) {
    get_true_values(all_xym = all_xym, iter = i, sample_size = 1e6,
                    nreps = nreps_total, y_only = FALSE, plasmode = TRUE,
                    plasmode_folder = paste0(proj_root, "Data/plasmode data sets/"))
  }, future.seed = seeds
)
all_output <- do.call(rbind.data.frame, output_list)
saveRDS(all_output, file = paste0(production_output_dir, "all_true_values_plasmode.rds"))
# average over the 10 replications
true_values_plasmode <- all_output %>% 
  group_by(estimand, type, variable, yscenario) %>% 
  summarize(truth = mean(value, na.rm = TRUE), .groups = "drop")
saveRDS(true_values_plasmode %>% 
          filter(type == "oracle", variable == "numEpiType") %>% 
          select(-type, -variable), 
        file = paste0(production_output_dir, "true_values_oracle_plasmode.rds"))
saveRDS(true_values_plasmode %>% 
          filter(type == "census", variable == "numEpiType") %>% 
          select(-type, -variable), 
        file = paste0(production_output_dir, "true_values_pop_plasmode.rds"))

saveRDS(true_values_plasmode %>% 
          filter(type == "oracle", grepl("AgeAtIndex", variable)) %>% 
          select(-type, -variable), 
        file = paste0(production_output_dir, "true_values_oracle_zs_plasmode.rds"))
saveRDS(true_values_plasmode %>% 
          filter(type == "census", grepl("AgeAtIndex", variable)) %>% 
          select(-type), 
        file = paste0(production_output_dir, "true_values_pop_zs_plasmode.rds"))

saveRDS(true_values_plasmode %>% 
          filter(type == "oracle", grepl("item9", variable)) %>% 
          select(-type), 
        file = paste0(production_output_dir, "true_values_oracle_ws_plasmode.rds"))
saveRDS(true_values_plasmode %>% 
          filter(type == "census", grepl("item9", variable)) %>% 
          select(-type), 
        file = paste0(production_output_dir, "true_values_pop_ws_plasmode.rds"))

