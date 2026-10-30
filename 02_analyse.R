###############################################################################
### Analysis of ONE simulated data set (step by step, same functions as 03)
### - Section 0: Settings, data and true values
### - Section 1: Data Preparation
### - Section 2: Summary Statistics
### - Section 3: Estimation of Weights
### - Section 4: Fit models (LWYY and Ghosh-Lin)
### - Section 5: Mean number of falls (nonparametric and model-based)
### - Section 6: Collect results
###############################################################################


## Load necessary libraries
library(tidyverse)   # Data manipulation and visualization
library(magrittr)    # Pipe operators and utilities
library(survival)    # LWYY (coxph)
library(mets)        # Ghosh-Lin (recreg)
library(ipw)         # Inverse Probability Weighting (ipwtm)
library(splines)     # ns() in the weight models

# functions (the same as in 03_multiple_analyse_cluster.R)
source("functions.R")


#####################################################
#----------- Section 0: Data and true values -------#
#####################################################
sim_log        <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

# choose the data set by its parameters, not by scenario number (numbers change with the design)
current_params <- sim_log %>%
  filter(scenario_name == "2_Informative", effect == "low", drop_beta_arm == -1, mortality_level == "none", sim_run == 1)
if (nrow(current_params) != 1) stop("data set not found (or not unique) in 00_simulation_log.csv")
current_file <- current_params$file_path
Odat_all     <- readRDS(current_file)

# true values per method (one row per causal scenario)
truth <- true_estimands %>%
  filter(Scenario_Name == current_params$scenario_name[1],
         Mortality     == current_params$mortality_level[1],
         True_Beta_Arm == current_params$true_beta_arm[1])
if (nrow(truth) != 1) stop("true values not found (or not unique) in 00_true_estimands.")

cat("data set:", current_file, "\n")
cat("conditional DGP parameter (true_beta_arm):", current_params$true_beta_arm[1], "\n")
cat("True_Beta_LWYY:   ", truth$True_Beta_LWYY, "\n")
cat("True_Beta_GL:     ", truth$True_Beta_GL, "\n")
cat("True_LogRateRatio:", truth$True_LogRateRatio, "(ratio of mean falls)\n")


#####################################################
#----------- Section 1: Data Preparation -----------#
#####################################################
# columns from the data generation:
# observed, n_falls, n_falls_full, drop, death, death_full, cum_falls_lag, train_rate_lag, frailty_t
# prep_odat() (functions.R) adds log_base_falls (log(1 + falls in the year before)), tstop (time for the
# weight models) and log_falls_lag (log(1 + observed falls up to week t-1))
Odat_all <- prep_odat(Odat_all)

Odat <- Odat_all %>% filter(observed)                            # observed data in the study


#######################################################
#----------- Section 2: Summary Statistics -----------#
#######################################################
# observed data (with informative dropout)
summary_stat_obs <- Odat %>%
  group_by(arm, ID) %>%
  summarise(falls = sum(n_falls), weeks = n(), dropped = max(drop), died = max(death), .groups = "drop") %>%
  group_by(arm) %>%
  summarise(total_pts      = n(),
            total_dropouts = sum(dropped),
            dropout_prop   = mean(dropped),
            deaths_obs     = sum(died),
            ttl_events_arm = sum(falls),
            ttl_time_arm   = sum(weeks),
            event_rate_arm = 52 * ttl_events_arm / ttl_time_arm,   # annualised event rate
            avg_evts       = ttl_events_arm / total_pts,
            avg_time       = ttl_time_arm / total_pts)
print(summary_stat_obs)

# proportion of all falls (complete data) that is observed in the study
summary_stat_full <- Odat_all %>%
  group_by(arm) %>%
  summarise(falls_full        = sum(n_falls_full),
            falls_observed    = sum(n_falls, na.rm = TRUE),
            prop_observed     = falls_observed / falls_full,
            falls_per_pt_full = falls_full / n_distinct(ID))
print(summary_stat_full)


##########################################################
#----------- Section 3: Estimation of Weights -----------#
##########################################################
# Stabilised weights with ipwtm, estimated separately per arm.
# - numerator  : baseline covariates + time (the same covariates are in the outcome model)
# - denominator: numerator + proxies (sym, sym_train) or true frailty (oracle)
# - ipwtm is fitted on weeks without death (dropout only possible if the person survives the week)
# - lag_ipw(): the weight for week t is the ipwtm weight of week t-1, because the falls of
#   week t are observed before a dropout at the end of week t
# Unstabilised weights (ipwtm without numerator) are needed for the mean number of falls (Section 5).

Odat_w <- compute_weights(Odat)



# check weights for plausibility
cat("\nStabilised weights (models) and unstabilised weights (MCF):\n")
print(summary(Odat_w[, c("w_sym", "w_symtrain", "w_oracle", "wm_sym", "wm_symtrain", "wm_oracle")]))

# mean of the stabilised weights per arm and week should be ~1
weight_check <- Odat_w %>%
  group_by(arm, t) %>%
  summarise(across(c(w_sym, w_symtrain, w_oracle), mean), .groups = "drop") %>%
  group_by(arm) %>%
  summarise(across(c(w_sym, w_symtrain, w_oracle), ~ paste(round(range(.x), 2), collapse = " - ")))
cat("\nRange of weekly mean weights per arm (should be close to 1):\n"); print(weight_check)

ess <- Odat_w %>%
  group_by(arm) %>%
  summarise(across(c(w_sym, w_symtrain, w_oracle), ~ sum(.x)^2 / sum(.x^2) / n()))
cat("\nEffective sample size relative to n person-weeks:\n"); print(ess)


###########################################################
#----------- Section 4: Fit models -----------------------#
###########################################################
# fit_models() (functions.R), the same 14 models as in 03:
# - counting-process data with prep_counting(): falls of the death / dropout week are kept,
#   weights are carried along per week
# - LWYY (coxph): death = censoring
# - Ghosh-Lin (recreg): death = terminal event, without death GL = LWYY
#   cens.model = ~ strata(arm): censoring Kaplan-Meier per arm (dropout differs between the arms)
#   control = gl_ctrl (tol 1e-12, Armijo line search): with the default tol 1e-6 recreg stops too early
# - each model: full data (oracle M0), naive, IPCW sym, sym_train, oracle, sym and sym_train with truncated weights
# NOTE: whether recreg handles row-specific (time-varying) weights correctly is not yet
#       verified -> compare with the nonparametric mean number of falls in Section 5
fits <- fit_models(Odat_all, Odat_w)

print(summary(fits[["LWYY 3: IPCW sym"]]))             # look at single models if needed



#############################################################
#----- Section 5: Mean number of falls ---------------------#
#############################################################
# mcf_estimate() (functions.R), as in 03:
# - nonparametric: mu_a(52) = sum_t S_a(t-1) * dR_a(t), weighted Nelson-Aalen + weighted Kaplan-Meier per arm
#   (weighted_mcf with the unstabilised weights wm_*)
# - model-based: standardised prediction of every LWYY / GL fit (model_mcf)
# true values: LWYY -> limit of the LWYY fit on complete data (MCF_LWYY_*, as True_Beta_LWYY),
#              GL and nonparametric MCF -> real world with death
mcf_est <- mcf_estimate(Odat_all, Odat_w, fits)

mean_falls_table <- tibble(Model            = names(mcf_est),
                           MCF_Control      = sapply(mcf_est, `[`, 1),
                           MCF_Intervention = sapply(mcf_est, `[`, 2)) %>%
  mutate(lwyy                  = startsWith(Model, "MCF-LWYY"),   # LWYY: limit on complete data
         True_MCF_Control      = if_else(lwyy, truth$MCF_LWYY_Control,      truth$MCF_Control),
         True_MCF_Intervention = if_else(lwyy, truth$MCF_LWYY_Intervention, truth$MCF_Intervention),
         True_LogRatio         = if_else(lwyy, truth$True_Beta_LWYY,        truth$True_LogRateRatio),
         Bias_MCF_Control      = MCF_Control - True_MCF_Control,
         Bias_MCF_Intervention = MCF_Intervention - True_MCF_Intervention,
         Log_Ratio             = log(MCF_Intervention / MCF_Control),
         Bias                  = Log_Ratio - True_LogRatio) %>%
  dplyr::select(-lwyy)
print(mean_falls_table)

# nonparametric MCF with mets (as in 03): mean falls per arm and robust SE of the log ratio
# check: MCF_*_mets = weighted_mcf (unweighted identical, weighted ~1e-5 because of tie_offset)
mcf_mets(Odat_all, Odat_w) %>%
  left_join(mean_falls_table %>% dplyr::select(Model, MCF_Control, MCF_Intervention), by = "Model") %>%
  mutate(diff_control      = MCF_Control_mets - MCF_Control,
         diff_intervention = MCF_Intervention_mets - MCF_Intervention) %>%
  print()

# time-varying effect: mean falls per quarter (as MCF-Q1..Q4 in 03), here for IPCW sym
mcf_q_sym <- weighted_mcf(Odat_w, "wm_sym", breaks = c(0, 13, 26, 39, 52)) %>%
  pivot_wider(names_from = arm, values_from = MCF, names_prefix = "arm") %>%
  mutate(Log_Ratio     = log(arm1 / arm0),
         True_LogRatio = unlist(truth[paste0("LogRR_Q", period)]))   # true log ratio per quarter
print(mcf_q_sym)



#############################################################
# ----- Section 6: Collect results -------------------------#
#############################################################
# effect estimates with the true value of each method (as in 03)
comparison_table <- bind_rows(lapply(names(fits), function(m) {
  true_beta <- if (startsWith(m, "LWYY")) truth$True_Beta_LWYY else truth$True_Beta_GL
  extract_model_results(fits[[m]], m, true_beta)
})) %>%
  mutate(scenario_id = current_params$scenario_id[1],      # for the loop over all files
         sim_run     = current_params$sim_run[1],
         file_path   = current_file)

print(comparison_table)

