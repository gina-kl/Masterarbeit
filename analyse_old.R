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
library(splines)

# functions
source("functions.R")

# Constants
n_pts <- 2778        # Total number of patients (both arms)
nboot <- 1000        # Number of bootstrap samples

# Load simulated data
current_file <- "./results/Odat_scen21_run001.rds"
Odat <- readRDS(current_file)

# load log data and true estimands
sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

# später schleife über sim_log
current_params <- sim_log %>% 
  filter(file_path == current_file)
true_beta_dgp <- current_params$true_beta_arm[1] # true beta from data generation process

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

cat("conditional DGP-Parameter (true_beta_arm):", true_beta_dgp, "\n")
cat("marginal true beta (True_LogRateRatio):    ", true_beta_marginal, "\n")



#####################################################
#----------- Section 1: Data Preparation -----------#
#####################################################
# cut off after death
Odat <- Odat %>%
  group_by(ID) %>%
  mutate(
    cutoff_t = pmin(dropout_time, death_time, real_followtime, na.rm = TRUE)
  ) %>%
  filter(t <= cutoff_t) %>%
  ungroup() %>%
  dplyr::select(-cutoff_t)



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
    dropout_event = ifelse(!is.na(dropout_time) & dropout_time == tstop, 1, 0), # Event indicator for Dropout
    cum_events = lag(cumsum(Yobs)),               # lag = time point before    
    cum_events = replace_na(cum_events, 0),
    cumTrain_lag = lag(CumTrain),                
    cumTrain_lag = replace_na(cumTrain_lag, 0)  
  ) %>%
  ungroup() %>% 
  as.data.frame()   
    
# split data for weights calculation

Odat_arm0 <- subset(Odat_counting, arm == 0)
Odat_arm1 <- subset(Odat_counting, arm == 1)

# IPW weights for control group  
ipw_arm0 <- ipwtm(
  exposure = dropout_event,        
  family = "binomial",             
  link = "logit",                  
  numerator = ~ age +sex,                                   
  denominator = ~ age + sex + cum_events,                    
  id = ID,                         
  tstart = tstart,                 
  timevar = tstop,                 
  type = "cens",  
  data = Odat_arm0              
)

# IPW weights for training group
ipw_arm1 <- ipwtm(
  exposure = dropout_event,        
  family = "binomial",             
  link = "logit",                  
  numerator = ~ age + sex ,                                   
  denominator = ~ age + sex + cum_events + cumTrain_lag,         
  id = ID,                         
  tstart = tstart,                 
  timevar = tstop,                 
  type = "cens", 
  data = Odat_arm1              
)

Odat_arm0$weights_cens <- ipw_arm0$ipw.weights
Odat_arm1$weights_cens <- ipw_arm1$ipw.weights


Odat_counting_weights <- bind_rows(Odat_arm0, Odat_arm1) %>%
  arrange(ID, tstop) %>%
  as.data.frame()

# check ipw wights for plausibality 
cat("\nDistribtion of IPCW-wights:\n")
print(summary(Odat_counting_weights$weights_cens))
cat("Number of weights > 10:", sum(Odat_counting_weights$weights_cens > 10), "\n")
cat("Number of weights > 20:", sum(Odat_counting_weights$weights_cens > 20), "\n")



###########################################################
#----------- Section 4: Fit LWYY Model -------------------#
###########################################################
# update status for LWYY and Gosh

Odat_counting_weights <- Odat_counting_weights %>%
  mutate(
    status_gl = ifelse(status == 2, 0, status),  # 0 = cens/Dropout, 1 = Fall, 3 = Death
    status_cox = ifelse(status == 1, 1, 0)       # 0 = cens/Dropout/Death, 1 = Fall
      )


# LWYY + IPW 
cox_weighted <- coxph(
  Surv(tstart, tstop, status_cox) ~ arm + age + sex  + cluster(ID), 
  data = Odat_counting_weights,  
  weights = weights_cens,
  robust = TRUE                                 
)

# LWYY without IPW
cox_unweighted <- coxph(
  Surv(tstart, tstop, status_cox) ~ arm + age + sex + cluster(ID), 
  data = Odat_counting_weights,  
  robust = TRUE                                 
)

#######################

Odat_counting_weights_gl <- Odat_counting_weights %>%
  mutate(
    # because of discrete simulation times (error)
    tstop_jitter = ifelse(status_gl %in% c(0, 3), tstop + 0.001, tstop)
  )
# Ghosh-Lin + IPW
gl_weighted <- recreg(
  Event(tstart, tstop_jitter, status_gl) ~ arm + age + sex + cluster(ID),
  data = Odat_counting_weights_gl,
  cause = 1,
  cens.code = 0,
  death.code = 3,
  weights = Odat_counting_weights_gl$weights_cens  
)

# Ghosh-Lin without IPW
gl_unweighted <- recreg(
  Event(tstart, tstop_jitter, status_gl) ~ arm + age + sex + cluster(ID),
  data = Odat_counting_weights_gl,
  cause = 1,
  cens.code = 0,
  death.code = 3
)


#############################################################
# -----Calculate mean number of falls--------------------------#
############################################################
# Use G computation to compute mean number of falls
# or with recurrent_marginal in mets package?

exposed <- Odat_counting_weights %>% mutate(arm = 1)       # Training group
unexposed <- Odat_counting_weights %>% mutate(arm = 0)     # Control group

# predict both outcomes for each patient
pred_exposed <- survfit(cox_weighted, newdata = exposed) 
pred_unexposed <- survfit(cox_weighted, newdata = unexposed)

# compute marginal effect
mcf_marginal <- data.frame(
  time = pred_exposed$time,
  mean_exposed = rowMeans(pred_exposed$cumhaz),
  mean_unexposed = rowMeans(pred_unexposed$cumhaz)
)

mcf_marginal%>%
  filter(time <= 52) %>%
  tail(1) %>%
  mutate(
    difference = mean_exposed - mean_unexposed,
    rate_ratio = mean_exposed / mean_unexposed
  )

pred_exposed_unweighted <- survfit(cox_unweighted, newdata = exposed) 
pred_unexposed_unweighted <- survfit(cox_unweighted, newdata = unexposed)

# compute marginal effect
mcf_marginal_unweighted <- data.frame(
  time = pred_exposed_unweighted$time,
  mean_exposed = rowMeans(pred_exposed_unweighted$cumhaz),
  mean_unexposed = rowMeans(pred_unexposed_unweighted$cumhaz)
)

mcf_marginal_unweighted%>%
  filter(time <= 52) %>%
  tail(1) %>%
  mutate(
    difference = mean_exposed - mean_unexposed,
    rate_ratio = mean_exposed / mean_unexposed
  )


#############################################################
# -----Extract results        --------------------------#
############################################################


# combine both results in table
comparison_table <- bind_rows(
  extract_model_results(cox_unweighted, "1: Naive (Unweighted)", true_beta = true_beta_marginal),
  extract_model_results(cox_weighted,   "2: LWYY + IPCW",        true_beta = true_beta_marginal)
)

# print
print(comparison_table)


comparison_table_dgp <- bind_rows(
  extract_model_results(cox_unweighted, "1: Naive (Unweighted)", true_beta = true_beta_dgp),
  extract_model_results(cox_weighted,   "2: LWYY + IPCW",        true_beta = true_beta_dgp)
)
cat("\nZum Vergleich -- Bias gegen konditionalen DGP-Parameter (nicht empfohlen als Hauptmass):\n")
print(comparison_table_dgp)


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
