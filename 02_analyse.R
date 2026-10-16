###############################################################################
### Analysis of ONE simulated data set
### - Section 0: Settings, data and true values
### - Section 1: Data Preparation
### - Section 2: Summary Statistics
### - Section 3: Estimation of Weights
### - Section 4: Fit models (LWYY and Ghosh-Lin)
### - Section 5: Mean number of falls (nonparametric, weighted)
### - Section 6: Collect results
###############################################################################


## Load necessary libraries
library(tidyverse)   # Data manipulation and visualization
library(magrittr)    # Pipe operators and utilities
library(survival)    # LWYY (coxph)
library(mets)        # Ghosh-Lin (recreg)
library(ipw)         # Inverse Probability Weighting (ipwtm)
library(splines)     # ns() in the weight models

# functions
source("functions.R")   

# Constants
t_end <- 52

#####################################################
#----------- Section 0: Data and true values -------#
#####################################################
current_file <- "./results/Odat_scen48_run0001.rds"
Odat_all <- readRDS(current_file)

sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

current_params <- sim_log %>% filter(file_path == current_file)

# true values per method (one row per causal scenario)
truth <- true_estimands %>%
  filter(Scenario_Name == current_params$scenario_name[1],
         Mortality     == current_params$mortality_level[1],
         True_Beta_Arm == current_params$true_beta_arm[1])
if (nrow(truth) != 1) stop("true values not found (or not unique) in 00_true_estimands.")

cat("conditional DGP parameter (true_beta_arm):", current_params$true_beta_arm[1], "\n")
cat("True_Beta_LWYY:   ", truth$True_Beta_LWYY, "\n")
cat("True_Beta_GL:     ", truth$True_Beta_GL, "\n")
cat("True_LogRateRatio:", truth$True_LogRateRatio, "(ratio of mean falls)\n")


#####################################################
#----------- Section 1: Data Preparation -----------#
#####################################################
# columns from the data generation:
# observed, n_falls, n_falls_full, drop, death, death_full, cum_falls_lag, train_rate_lag, frailty_t
Odat_all <- Odat_all %>%
# log_base_falls = log(1 + falls in the year before)
# The difference becomes smaller as the number of iterations increases, if the values are too high, recreg does not converge
  mutate(log_base_falls = log1p(base_risk),                       # log(1 + falls in the year before)
         tstop         = t,                                     # time variable for the weight models
         fall_rate_lag  = cum_falls_lag / pmax(t - 1, 1)) %>%    # observed fall rate up to week t-1
  arrange(ID, t)

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
cat("\nStabilised weights (used in models):\n")                
print(summary(Odat_w[, c("w_sym", "w_symtrain", "w_oracle")]))

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

# counting-process data for LWYY and Ghosh-Lin with prep_counting(): falls of the death / dropout week
# are kept, weights are carried along per week
cp      <- prep_counting(Odat_w, "observed")        # observed data incl. weights
cp_full <- prep_counting(Odat_all, "full")          # complete data (oracle M0)

# --- LWYY (death = censoring) ---------------------------------------------------
f_lwyy <- Surv(start, stop, status == 1) ~ arm + age + sex + log_base_falls + cluster(ID)

lwyy_full     <- coxph(f_lwyy, data = cp_full, ties = "breslow")
lwyy_naive    <- coxph(f_lwyy, data = cp,      ties = "breslow")
lwyy_sym      <- coxph(f_lwyy, data = cp,      ties = "breslow", weights = w_sym)
lwyy_symtrain <- coxph(f_lwyy, data = cp,      ties = "breslow", weights = w_symtrain)
lwyy_oracle   <- coxph(f_lwyy, data = cp,      ties = "breslow", weights = w_oracle)

# --- Ghosh-Lin (death = terminal event) -----------------------------------------
# fitted in all scenarios, without death GL = LWYY
# NOTE: whether recreg handles row-specific (time-varying) weights correctly is not yet
#       verified -> compare with the nonparametric mean number of falls in Section 5
f_gl <- Event(start, stop, status) ~ arm + age + sex + log_base_falls + cluster(ID)

# cens.model = ~ strata(arm):  recreg weights the persons who died with the censoring
#   probabilities G_c (Kaplan-Meier from the rows with status 0, ?recreg: "Formula for a stratified Cox model
#   without covariates used to estimate censoring probabilities"). Default ~1 = one Kaplan-Meier for both arms;
#   dropout differs strongly between the arms -> one Kaplan-Meier per arm, as for the IPCW models (per arm)
#   and as in the example of ?recreg ("censoring stratified after treatment").
#   The IPCW weights (weights = ...) are a different correction: informative dropout given the history.
gl_full     <- recreg(f_gl, data = cp_full, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm))
gl_naive    <- recreg(f_gl, data = cp,      cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm))
gl_sym      <- recreg(f_gl, data = cp,      cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), weights = cp$w_sym)
gl_symtrain <- recreg(f_gl, data = cp,      cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), weights = cp$w_symtrain)
gl_oracle   <- recreg(f_gl, data = cp,      cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), weights = cp$w_oracle)



#############################################################
#----- Section 5: Mean number of falls (nonparametric) -----#
#############################################################
# MCF for LWYY
#   mu_a(52) = sum_t S_a(t-1) * dR_a(t), weighted Nelson-Aalen + weighted Kaplan-Meier per arm
#   true values: MCF_Control, MCF_Intervention, True_LogRateRatio
mcf_list <- list(
  "0: Full data (oracle M0)" = weighted_mcf(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full"),
  "1: Naive (Unweighted)"    = weighted_mcf(Odat_w %>% mutate(w1 = 1), "w1"),
  "3: IPCW sym"              = weighted_mcf(Odat_w, "wm_sym"),
  "4: IPCW sym_train"        = weighted_mcf(Odat_w, "wm_symtrain"),
  "5: IPCW oracle"           = weighted_mcf(Odat_w, "wm_oracle")
)


# baseline covariates, one row per person (for the standardised model-based MCFs)
base_cov <- Odat_all %>% filter(t == 1) %>% dplyr::select(ID, age, sex, log_base_falls)

gl_list   <- list("0: Full data (oracle M0)" = gl_full,   "1: Naive (Unweighted)" = gl_naive,
                  "3: IPCW sym" = gl_sym,   "4: IPCW sym_train" = gl_symtrain,   "5: IPCW oracle" = gl_oracle)
lwyy_list <- list("0: Full data (oracle M0)" = lwyy_full, "1: Naive (Unweighted)" = lwyy_naive,
                  "3: IPCW sym" = lwyy_sym, "4: IPCW sym_train" = lwyy_symtrain, "5: IPCW oracle" = lwyy_oracle)


# table: nonparametric, GL- and LWYY-based MCF with the matching true values
# function to compute mcf lwyy and gl see function.R
mean_falls_table <- bind_rows(
  imap_dfr(mcf_list, function(m, name)
    tibble(Model = paste("MCF", name), MCF_Control = m$MCF[m$arm == 0], MCF_Intervention = m$MCF[m$arm == 1])),
  imap_dfr(gl_list, function(fit, name) {
    mu <- model_mcf(fit, base_cov, t_end)
    tibble(Model = paste("MCF-GL", name), MCF_Control = mu[1], MCF_Intervention = mu[2])
  }),
  imap_dfr(lwyy_list, function(fit, name) {
    mu <- model_mcf(fit, base_cov, t_end)
    tibble(Model = paste("MCF-LWYY", name), MCF_Control = mu[1], MCF_Intervention = mu[2])
  })
) %>%
  mutate(no_death              = startsWith(Model, "MCF-LWYY"),   # LWYY: hypothetical world without death
         True_MCF_Control      = if_else(no_death, truth$MCF_Control_NoDeath,       truth$MCF_Control),
         True_MCF_Intervention = if_else(no_death, truth$MCF_Intervention_NoDeath,  truth$MCF_Intervention),
         True_LogRatio         = if_else(no_death, truth$True_LogRateRatio_NoDeath, truth$True_LogRateRatio),
         Bias_MCF_Control      = MCF_Control - True_MCF_Control,
         Bias_MCF_Intervention = MCF_Intervention - True_MCF_Intervention,
         Log_Ratio             = log(MCF_Intervention / MCF_Control),
         Bias                  = Log_Ratio - True_LogRatio) %>%
  dplyr::select(-no_death)
print(mean_falls_table)





#############################################################
# ----- Section 6: Collect results -------------------------#
#############################################################
comparison_table <- bind_rows(
  extract_model_results(lwyy_full,     "LWYY 0: Full data (oracle M0)", truth$True_Beta_LWYY),
  extract_model_results(lwyy_naive,    "LWYY 1: Naive (Unweighted)",    truth$True_Beta_LWYY),
  extract_model_results(lwyy_sym,      "LWYY 3: IPCW sym",              truth$True_Beta_LWYY),
  extract_model_results(lwyy_symtrain, "LWYY 4: IPCW sym_train",        truth$True_Beta_LWYY),
  extract_model_results(lwyy_oracle,   "LWYY 5: IPCW oracle",           truth$True_Beta_LWYY),
  extract_model_results(gl_full,       "GL 0: Full data (oracle M0)",   truth$True_Beta_GL),
  extract_model_results(gl_naive,      "GL 1: Naive (Unweighted)",      truth$True_Beta_GL),
  extract_model_results(gl_sym,        "GL 3: IPCW sym",                truth$True_Beta_GL),
  extract_model_results(gl_symtrain,   "GL 4: IPCW sym_train",          truth$True_Beta_GL),
  extract_model_results(gl_oracle,     "GL 5: IPCW oracle",             truth$True_Beta_GL)
) %>%
  mutate(scenario_id = current_params$scenario_id[1],      # for the loop over all files
         sim_run     = current_params$sim_run[1],
         file_path   = current_file)

print(comparison_table)

