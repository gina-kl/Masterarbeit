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
library(mets)

# Constants
n_pts <- 2778        # Total number of patients (both arms)
nboot <- 1000        # Number of bootstrap samples

# Load simulated data
current_file <- "./results/Odat_scen08_run001.rds"
Odat <- readRDS("./results/Odat_scen08_run001.rds")

# load log data
sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
# später schleife über sim_log
current_params <- sim_log %>% 
  filter(file_path == current_file)
true_beta <- current_params$true_beta_arm[1]


#####################################################
#----------- Section 1: Data Preparation -----------#
#####################################################


Odat_prepared <- Odat %>%
  group_by(ID) %>%
  mutate(
    ttl_events_pts = sum(Yobs, na.rm = TRUE), # Sum of events per patient
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
    dropout_event = ifelse(!is.na(dropout_time) & dropout_time == tstop, 1, 0) # Event indicator for Dropout
  ) %>%
  ungroup() %>% 
  as.data.frame()                  # ipw package strictly need data.frame

# Estimate inverse probability of censoring weights (IPCW)
# Here we model the probability of dropping out based on observed covariates
ipwtm_cens <- ipwtm(
  exposure = dropout_event,        # Dependent variable: dropout indicator
  family = "binomial",             # Logistic regression
  link = "logit",                  
  numerator = ~ arm + age + base_risk,                     # Model for numerator (baseline covariates only)
  denominator = ~ arm + age + base_risk + CumTrain,        # Full model including time-varying covariates (CumTrain)
  id = ID,                         # Cluster identifier
  tstart = tstart,                 # Start of time interval
  timevar = tstop,                 # End of time interval
  type = "cens",                   # Type of weights: Censoring weights (IPCW)
  data = Odat_counting             
)

# Incorporate weights into the counting process data
Odat_counting_weights <- Odat_counting %>%
  mutate(weights_cens = ipwtm_cens$ipw.weights) %>%  # Add estimated weights
  group_by(ID) %>%
  mutate(
    cum_events = lag(cumsum(Yobs)),                  # Cumulative events up to previous time point
    cum_events = replace_na(cum_events, 0)           # Replace NA with 0 for the first row
  ) %>%
  ungroup()

Odat_counting_weights <- as.data.frame(Odat_counting_weights) # as data.frame for coxph


###########################################################
#----------- Section 4: Fit LWYY Model -------------------#
###########################################################

# LWYY + IPW (Marginal Cox Model for Recurrent Events)

cox_weighted <- coxph(
  Surv(tstart, tstop, Yobs) ~ arm + age + base_risk + cluster(ID), 
  data = Odat_counting_weights,  
  #cluster = ID,                                
  weights = weights_cens,
  robust = TRUE                                 
)

#  LWYY without IPW
cox_unweighted <- coxph(
  Surv(tstart, tstop, Yobs) ~ arm + age + base_risk + cluster(ID), 
  data = Odat_counting_weights,  
  #cluster = ID,                                
  #weights = weights_cens,
  robust = TRUE                                 
)

# # Gosh-Lin without IPW
# # Problem mit Dropout! zu zensierung?
# gl_unweighted <- recreg(
#   Event(tstart, tstop, status) ~ arm + age + base_risk + cluster(ID), 
#   data = Odat_counting_weights,
#   cause = 1,       # recurrent event 
#   death.code = 3   # death
# )
# 
# # Gosh-Lin + IPW
# gl_weighted <- recreg(
#   Event(tstart, tstop, status) ~ arm + age + base_risk, 
#   data = Odat_counting_weights,
#   cause = 1,
#   death.code = 3,
#   weights = Odat_counting_weights$weights_cens, 
#   id = Odat_counting_weights$ID
# )

# Extract summary statistics
summary_cox <- summary(cox_weighted)
coef_cox_weighted <- summary_cox$coefficients["arm", "coef"]            # Log Hazard Ratio for treatment
se_cox_weighted_naive <- summary_cox$coefficients["arm", "robust se"]   # Robust standard error
p_val_cox_weighted <- summary_cox$coefficients["arm", "Pr(>|z|)"]       # P-Value
rej_cox_weighted_naive <- ifelse(p_val_cox_weighted < 0.05, 1, 0)       # Reject Null Hypothesis (1 = Yes)

# bias
#bias_unweighted <- coef_unweighted - true_beta
bias_weighted   <- coef_cox_weighted - true_beta


cat("\n--- LWYY + IPCW Results ---\n")
cat("Treatment Effect (Log HR):", round(coef_cox_weighted, 4), "\n")
cat("Robust SE:              ", round(se_cox_weighted_naive, 4), "\n")
cat("P-Value:                ", round(p_val_cox_weighted, 4), "\n")
cat("Significant (alpha=0.05)?", ifelse(rej_cox_weighted_naive == 1, "Yes", "No"), "\n")





#############################################################
# -----------------Plots-----------------------------------#
############################################################
# Cumulated falls over time
Odat %>%
  group_by(ID) %>%
  mutate(cum_y = cumsum(Yobs)) %>%
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
  geom_point(data = . %>% filter(Yobs == 1), aes(color = "Sturz"), size = 2) +
  geom_point(data = . %>% filter(status == 2), aes(color = "Dropout"), shape = 4, size = 3, stroke = 1.5) +
  scale_color_manual(values = c("Sturz" = "#D55E00", "Dropout" = "black")) +
  labs(
    title = "Individuelle Verläufe von 20 zufälligen Patienten",
    x = "Beobachtungswoche (t)",
    y = "Patienten-ID",
    color = "Ereignis"
  ) +
  theme_minimal(base_size = 12)
