###############################################################################
### Analysis for Scenario 1 with the Time-Varying Covariate Measured Weekly ###
### - Section 1: Data Preparation
### - Section 2: Summary Statistics
### - Section 3: Estimation of Weights
### - Section 4: Fit LWYY and NB Models
### - Section 5: Bootstrap for three IPW approaches to calculate variance 
### - Section 6: Collect and Save Results
###############################################################################

## Load necessary libraries
library(ipw)         # Inverse Probability Weighting
library(MASS)        # For fitting negative binomial models
library(tidyverse)   # Data manipulation and visualization
library(magrittr)    # Pipe operators and utilities
library(survival)    # Recurrent event analysis
library(mets)        # For Gosh-Lin Model
library(splines)     # For time-dependent weights in ipw

# functions
source("functions.R")

# Constants
n_pts <- 2778        # Total number of patients (both arms)
nboot <- 1000        # Number of bootstrap samples
trunc_level <- 0.01  # IPW truncation: weights below 1st / above 99th percentile are set to these percentiles

# Load simulated data
current_file <- "./results/Odat_scen14_run001.rds"
Odat_all <- readRDS(current_file)   

# load log data and true estimands
sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

# später schleife über sim_log
current_params <- sim_log %>% 
  filter(file_path == current_file)
true_beta_dgp <- current_params$true_beta_arm[1] # true beta from data generation process

# select true marginal effect
true_beta_marginal <- true_estimands %>%
  filter(
    Scenario_Name == current_params$scenario_name[1],
    Mortality     == current_params$mortality_level[1],
    True_Beta_Arm == current_params$true_beta_arm[1]
  ) %>%
  pull(True_LogRateRatio)

if (length(true_beta_marginal) != 1) {
  stop("not found in 00_true_estimands.")
}

# true value without death (for LWYY)
true_beta_nodeath <- true_estimands %>%
  filter(
    Scenario_Name == current_params$scenario_name[1],
    Mortality     == current_params$mortality_level[1],
    True_Beta_Arm == current_params$true_beta_arm[1]
  ) %>%
  pull(True_LogRateRatio_NoDeath)

cat("conditional DGP-Parameter (true_beta_arm):", true_beta_dgp, "\n")
cat("marginal true beta (True_LogRateRatio):    ", true_beta_marginal, "\n")
cat("marginal true beta without death:          ", true_beta_nodeath, "\n")

# method depends on mortality: LWYY only for low, Ghosh-Lin only for high mortality
method <- ifelse(current_params$mortality_level[1] == "low", "LWYY", "GL")
cat("Method:", method, "\n")



#####################################################
#----------- Section 1: Data Preparation -----------#
#####################################################
# observed data in the study
Odat <- Odat_all %>% filter(observed)

# complete data for oracle M0 (no dropout, no random censoring, only death or t.end)
Odat_full <- Odat_all %>%
  filter(t != 0) %>%
  mutate(tstart = t - 1, tstop = t) %>%
  as.data.frame()



Odat_prepared <- Odat %>%
  group_by(ID) %>%
  mutate(
    ttl_events_pts = sum(Yrep, na.rm = TRUE), # Sum of reported events per patient
    ttl_time_pts = max(t)                     # Maximum follow-up time
  ) %>%
  ungroup()


#######################################################
#----------- Section 2: Summary Statistics -----------#
#######################################################

# Summary stats for the actual observed data (with informative dropout)
summary_stat_obs <- Odat_prepared %>%
  group_by(ID) %>%
  slice(1) %>%  # Select the first record for each patient
  ungroup() %>%
  group_by(arm) %>%                            # In each arm, calculate:
  mutate(
    total_pts = n(),                           # Total number of patients 
    total_dropouts = sum(ever_dropped),        # Number of patients who dropped out
    ttl_events_arm = sum(ttl_events_pts),      # Total events in the group
    ttl_time_arm = sum(ttl_time_pts),          # Total follow-up time in the group
    event_rate_arm = 52 * ttl_events_arm / ttl_time_arm,  # Annualized event rate (52 weeks)
    avg_evts = ttl_events_arm / total_pts,     # Average events per patient
    dropout_prop = total_dropouts / total_pts, # Proportion of dropouts
    avg_time = ttl_time_arm / total_pts        # Average follow-up time per patient
  ) %>%
  slice(1) %>%  # Retain one row per treatment arm
  ungroup() %>%
  dplyr::select(arm, total_pts, total_dropouts, dropout_prop, ttl_events_arm, ttl_time_arm, event_rate_arm, avg_evts, avg_time)

print(summary_stat_obs)

# proportion of all falls (complete data) that is observed in the study
summary_stat_full <- Odat_all %>%
  group_by(arm) %>%
  summarise(
    falls_full = sum(Yfull),                        # all falls
    falls_observed = sum(Yobs[observed]),           # falls observed in study
    prop_observed = falls_observed / falls_full,
    falls_per_pt_full = falls_full / n_distinct(ID)
  )

print(summary_stat_full)


##########################################################
#----------- Section 3: Estimation of Weights -----------#
##########################################################

# Prepare counting process data for weight estimation
# tstart and tstop are required
Odat_counting <- Odat_prepared %>%
  group_by(ID) %>% 
  filter(t != 0) %>%  # Exclude the baseline (t = 0)
  mutate(
    tstart = t - 1,                # Start time of each interval 
    tstop = t,                     # Stop time of each interval
    dropout_event = ifelse(!is.na(dropout_time) & dropout_time == tstop, 1, 0), # Event indicator for Dropout
    cum_events = lag(cumsum(Yrep)),               # reported falls (lag = until week t-1)
    cum_events = replace_na(cum_events, 0),
    cumTrain_lag = lag(CumTrain),
    cumTrain_lag = replace_na(cumTrain_lag, 0),
    cumContact_lag = lag(CumContact),             # answered calls until t-1 (only control)
    cumContact_lag = replace_na(cumContact_lag, 0),
    fall_rate_lag    = cum_events     / pmax(tstart, 1),  # fall rate per week
    train_rate_lag   = cumTrain_lag   / pmax(tstart, 1),
    contact_rate_lag = cumContact_lag / pmax(tstart, 1)
  ) %>%
  ungroup() %>% 
  as.data.frame()   

# split data for weights calculation

Odat_arm0 <- subset(Odat_counting, arm == 0)
Odat_arm1 <- subset(Odat_counting, arm == 1)

# Three different weights:
#  - falls:  only (reported) falls as proxy
#  - proxy:  falls + Proxy for Frailty (Treatment group: Training, control group: answered calls)
#  - oracle: true Frailty (not possible, only to see what's possible)
#
# ns(tstop, df = 3) : later rates are more important, because more information

# Control group: falls
ipw_arm0_falls <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) * fall_rate_lag,
  id = ID,
  tstart = tstart,  # is ignored with family "binomial"
  timevar = tstop,
  type = "cens",
  data = Odat_arm0
)

# Control group: proxy (reported calls and falls)
ipw_arm0_proxy <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) * (fall_rate_lag + contact_rate_lag),
  id = ID,
  tstart = tstart,
  timevar = tstop,
  type = "cens",
  data = Odat_arm0
)

# Control group: oracle
ipw_arm0_oracle <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) + frailty,
  id = ID,
  tstart = tstart,
  timevar = tstop,
  type = "cens",
  data = Odat_arm0
)

# Training group: falls
ipw_arm1_falls <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) * fall_rate_lag,
  id = ID,
  tstart = tstart,
  timevar = tstop,
  type = "cens",
  data = Odat_arm1
)

# Training group: proxy (falls and training)
ipw_arm1_proxy <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) * (fall_rate_lag + train_rate_lag),
  id = ID,
  tstart = tstart,
  timevar = tstop,
  type = "cens",
  data = Odat_arm1
)

# Training group: oracle
ipw_arm1_oracle <- ipwtm(
  exposure = dropout_event,
  family = "binomial",
  link = "logit",
  numerator = ~ age + sex + ns(tstop, df = 3),
  denominator = ~ age + sex + ns(tstop, df = 3) + frailty,
  id = ID,
  tstart = tstart,
  timevar = tstop,
  type = "cens",
  data = Odat_arm1
)

# truncation as extra function, because in sumulation as variable. Not possible in trunc in ipw
# evtl trunc_level auf 0.005
Odat_arm0$w_falls       <- trunc_weights(ipw_arm0_falls$ipw.weights, trunc_level)
Odat_arm0$w_proxy       <- trunc_weights(ipw_arm0_proxy$ipw.weights, trunc_level)
Odat_arm0$w_oracle      <- trunc_weights(ipw_arm0_oracle$ipw.weights, trunc_level)

Odat_arm1$w_falls       <- trunc_weights(ipw_arm1_falls$ipw.weights, trunc_level)
Odat_arm1$w_proxy       <- trunc_weights(ipw_arm1_proxy$ipw.weights, trunc_level)
Odat_arm1$w_oracle      <- trunc_weights(ipw_arm1_oracle$ipw.weights, trunc_level)


Odat_counting_weights <- bind_rows(Odat_arm0, Odat_arm1) %>%
  arrange(ID, tstop) %>%
  mutate(weights_cens = w_proxy) %>%    
  as.data.frame()

# check ipw weights for plausibility
cat("\nDistribution of IPCW-weights before truncation:\n")
print(summary(data.frame(
  w_falls  = c(ipw_arm0_falls$ipw.weights,  ipw_arm1_falls$ipw.weights),
  w_proxy  = c(ipw_arm0_proxy$ipw.weights,  ipw_arm1_proxy$ipw.weights),
  w_oracle = c(ipw_arm0_oracle$ipw.weights, ipw_arm1_oracle$ipw.weights)
)))
cat("\nDistribution of IPCW-weights after truncation (used in models):\n")
print(summary(Odat_counting_weights[, c("w_falls", "w_proxy", "w_oracle")]))
cat("Number of weights > 10 (proxy):", sum(Odat_counting_weights$w_proxy > 10), "\n")
cat("Number of weights > 20 (proxy):", sum(Odat_counting_weights$w_proxy > 20), "\n")



###########################################################
#----------- Section 4: Fit models -----------------------#
###########################################################
# mortality low: LWYY, mortality high: Ghosh-Lin
# M0 (complete data, oracle) added for both
# ties = "breslow": because of possible ties (with poi~falls normally no ties left)

if (method == "LWYY") {

  # not include week of death/dropout (status =2 or 3)
  # only status 0/1 left
  set.seed(1)   
  Odat_lwyy <- Odat_counting_weights %>% filter(status %in% c(0, 1)) %>%
    expand_falls("status", "Yrep")
  Odat_lwyy_full <- Odat_full %>% filter(status_full %in% c(0, 1)) %>%
    expand_falls("status_full", "Yfull")  

  # LWYY on complete data (oracle)
  fit_full <- coxph(Surv(tstart, tstop, status_full) ~ arm + age + sex + cluster(ID),
                    data = Odat_lwyy_full, ties = "breslow", robust = TRUE)
  # LWYY without IPW
  fit_naive <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                     data = Odat_lwyy, ties = "breslow", robust = TRUE)
  # LWYY + IPW (falls)
  fit_falls <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                     data = Odat_lwyy, weights = w_falls, ties = "breslow", robust = TRUE)
  # LWYY + IPW (proxy) 
  fit_weighted <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                        data = Odat_lwyy, weights = w_proxy, ties = "breslow", robust = TRUE)
  # LWYY + IPW (oracle)
  fit_oracle <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                      data = Odat_lwyy, weights = w_oracle, ties = "breslow", robust = TRUE)
}

if (method == "GL") {

  # data in the form recreg needs 
  set.seed(1)  
  Odat_gl <- prep_gl_data(expand_falls(Odat_counting_weights, "status", "Yrep"),
                               status_col = "status")
  Odat_gl_full <- prep_gl_data(expand_falls(Odat_full, "status_full", "Yfull"),
                               status_col = "status_full")

  # Ghosh-Lin on complete data (oracle)
  fit_full <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                     data = Odat_gl_full, cause = 1, death.code = 3, cens.code = 0)
  # Ghosh-Lin without IPW
  fit_naive <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                      data = Odat_gl, cause = 1, death.code = 3, cens.code = 0)
  # Ghosh-Lin + IPW (falls)
  fit_falls <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                      data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                      weights = Odat_gl$w_falls)
  # Ghosh-Lin + IPW (proxy)
  fit_weighted <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                         data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                         weights = Odat_gl$w_proxy)
  # Ghosh-Lin + IPW (oracle)
  fit_oracle <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                       data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                       weights = Odat_gl$w_oracle)
}


#############################################################
# -----Calculate mean number of falls--------------------------#
############################################################
# Use G computation to compute mean number of falls
# only for LWYY (survfit needs a coxph model) TODO: Gosh-Lin

if (method == "LWYY") {

  # one row per patient
  baseline_pts <- Odat_counting_weights %>% distinct(ID, age, sex)
  exposed   <- baseline_pts %>% mutate(arm = 1)     # Training group
  unexposed <- baseline_pts %>% mutate(arm = 0)     # Control group
  
  fits <- list("0: Full data (oracle M0)" = fit_full,
               "1: Naive (Unweighted)"    = fit_naive,
               "2: IPCW falls"            = fit_falls,
               "3: IPCW proxy"            = fit_weighted,
               "4: IPCW oracle"           = fit_oracle)
  
  mean_falls_table <- map_dfr(names(fits), function(m) {
    pred_exposed   <- survfit(fits[[m]], newdata = exposed, se.fit = FALSE)
    pred_unexposed <- survfit(fits[[m]], newdata = unexposed, se.fit = FALSE)
    
    data.frame(
      time           = pred_exposed$time,
      mean_exposed   = rowMeans(pred_exposed$cumhaz),
      mean_unexposed = rowMeans(pred_unexposed$cumhaz)
    ) %>%
      filter(time <= 52) %>%
      tail(1) %>%
      mutate(Model      = paste(method, m),
             difference = mean_exposed - mean_unexposed,
             rate_ratio = mean_exposed / mean_unexposed)
  })
  
  print(mean_falls_table)
}


#############################################################
# -----Extract results        --------------------------#
############################################################


# combine both results in table
# fit_* for LWYY or GL, M0 added, method in model name
comparison_table <- bind_rows(
  extract_model_results(fit_full,     paste(method, "0: Full data (oracle M0)"), true_beta_marginal, true_beta_nodeath),
  extract_model_results(fit_naive,    paste(method, "1: Naive (Unweighted)"),    true_beta_marginal, true_beta_nodeath),
  extract_model_results(fit_falls,    paste(method, "2: IPCW falls"),            true_beta_marginal, true_beta_nodeath),
  extract_model_results(fit_weighted, paste(method, "3: IPCW proxy"),            true_beta_marginal, true_beta_nodeath),
  extract_model_results(fit_oracle,   paste(method, "4: IPCW oracle"),           true_beta_marginal, true_beta_nodeath)
)

# print
print(comparison_table)




#############################################################
# -----------------Plots-----------------------------------#
############################################################
# Cumulated falls over time
Odat %>%
  group_by(ID) %>%
  mutate(cum_y = cumsum(Yrep)) %>%          
  group_by(arm, t) %>%
  summarise(mean_cum_y = mean(cum_y), .groups = "drop") %>%
  ggplot(aes(x = t, y = mean_cum_y, color = factor(arm))) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("0" = "#D55E00", "1" = "#0072B2"),
    labels = c("0" = "Kontrolle", "1" = "Training")
  ) +
  labs(
    title = "Durchschnittlich kumulierte Stürze pro Patient",
    x = "Beobachtungswoche (t)",
    y = "Kumulierte Stürze (Mittelwert)",
    color = "Gruppe"
  ) +
  theme_minimal(base_size = 12)

# mean cumulative falls, complete vs. observed data
# per week: falls / patients at risk, then summed up over the weeks
bind_rows(
  Odat_all %>% filter(t > 0) %>% mutate(falls = Yfull, data = "vollständig"),
  Odat     %>% filter(t > 0) %>% mutate(falls = Yrep,  data = "beobachtet")
) %>%
  group_by(data, arm, t) %>%
  summarise(rate = sum(falls) / n(), .groups = "drop") %>%
  group_by(data, arm) %>%
  arrange(t) %>%
  mutate(mcf = cumsum(rate)) %>%
  ggplot(aes(x = t, y = mcf, color = factor(arm), linetype = data)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c("0" = "#D55E00", "1" = "#0072B2"),
    labels = c("0" = "Kontrolle", "1" = "Training")
  ) +
  labs(
    title = "Mittlere kumulierte Stürze: vollständige vs. beobachtete Daten",
    x = "Woche (t)",
    y = "Stürze pro Patient",
    color = "Gruppe",
    linetype = "Daten"
  ) +
  theme_minimal(base_size = 12)

# Time in study
Odat %>%
  group_by(arm, t) %>%
  summarise(n_active = n(), .groups = "drop") %>%
  group_by(arm) %>%
  mutate(prop_active = n_active / first(n_active)) %>%
  ggplot(aes(x = t, y = prop_active, color = factor(arm))) +
  geom_step(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(
    values = c("0" = "#D55E00", "1" = "#0072B2"),
    labels = c("0" = "Kontrolle", "1" = "Training")
  ) +
  labs(
    title = "Verbleib der Patienten im Studienverlauf",
    x = "Beobachtungswoche (t)",
    y = "Anteil aktiver Patienten",
    color = "Gruppe"
  ) +
  theme_minimal(base_size = 12)

# some individual dropouts
set.seed(42)
sample_ids <- sample(unique(Odat$ID), 20)

Odat %>%
  filter(ID %in% sample_ids) %>%
  ggplot(aes(x = t, y = factor(ID), group = ID)) +
  geom_line(color = "grey70", linewidth = 0.8) +
  geom_point(data = . %>% filter(Yrep >= 1), aes(color = "Sturz"), size = 2) +
  geom_point(data = . %>% filter(status == 2), aes(color = "Dropout"), shape = 4, size = 3, stroke = 1.5) +
  scale_color_manual(values = c("Sturz" = "#D55E00", "Dropout" = "black")) +
  labs(
    title = "Individuelle Verläufe von 20 zufälligen Patienten",
    x = "Beobachtungswoche (t)",
    y = "Patienten-ID",
    color = "Ereignis"
  ) +
  theme_minimal(base_size = 12)
