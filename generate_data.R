############################################################################
### Scenario 1: Recurrent Event Probability Independent of Event History ###
############################################################################
# Adjustment to SURE-FIT ###################################################



# Parse command-line arguments
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("The analysis tag needs to be provided! Exiting...\n")
# }
# 
# seed <- as.numeric(args[1])
seed <- 123

 # Load required libraries
 library(simcausal)  # To simulate time-varying data#
 library(ipw)        # To estimate the weights for IPCW
 library(MASS)       # To fit the negative binomial model
 library(tidyverse)
 library(magrittr)
 library(survival)
 library(reshape2)

 # Simulation parameters
 n_pts <- 2778  # Total number of patients from both arms
 t.end <- 52  # Study duration in weeks (1 year)
 entry_time <- runif(n_pts, 0, 12)  # Entry time for each patient
 exp_followtime <- rexp(n_pts, rate = 0.001) + 1  # min 1 week follow-up
 max_followtime <- t.end - entry_time # maximal follow-up time
 real_followtime <- pmin(exp_followtime, max_followtime) %>% round()  # real follow up time
 sum(exp_followtime < max_followtime) / n_pts # Independent Dropout proportion


# Initialize an empty DAG
D <- DAG.empty()

# Define nodes at baseline (t = 0)
D <- D +
  node("Frailty", t=0, distr = "rnorm", mean=0, sd=1) + # unobserved motivation/variable
  node("B3", t = 0, distr = "runif", min = 65, max = 95) +    # age
  node("L1", t = 0, distr = "rnorm", mean = 50, sd = 10) +    # Baseline Fall-Risk
  node("A1", t = 0, distr = "rbern", prob = 0.5) +            # 2 arms, 1 = Training, 0=no training
  node("Dropout", t = 0, distr = "rbern", prob = 0) +  # no dropouts at day 1
  node("CumTrain", t = 0, distr = "rconst", const = 0) + # start with 0 trainings
  node("Yobs", t = 0, distr = "rconst", const = 0) +
  node("Train", t=0, distr = "rconst", const =0 )

# Define nodes for subsequent time points (t > 0)
D <- D +
  node("A1", t = 1:t.end, distr = "rbern",
       prob = ifelse(A1[t-1] == 1, 1, 0)) + # Treatment assignment remains the same as baseline (time-independent)
  node("Dropout", t = 1:t.end, distr = "rbern",
       prob = ifelse(Dropout[t-1] == 1, 1,
                     plogis(-4 + 0.03 * B3[0] - 0.6 * Frailty[0]))) + # informative dropout (Patients who are more seriously ill are more likely to drop out)
  node("Train", t = 1:t.end, distr = "rbern",
       prob = ifelse(A1[0] == 0 | Dropout[t] == 1, 0, # control group or dropout
                     plogis(
                       -0.5                        # baseline probability for training (~38%)
                       - 0.02 * B3[0]              # baseline: older patients train less
                       + 0.50 * Frailty[0]         # unobserved motivation
                       + 0.05 * CumTrain[t-1]      # history of trainings
                     ))) +
  node("CumTrain", t = 1:t.end, distr = "rconst",  # update cumulative trainings
       const = CumTrain[t-1] + Train[t]) +
  node("Yobs", t = 1:t.end, distr = "rbern",
       prob = plogis(-3 + 0.02 * B3[0] + 0.4 * Frailty[0] - 0.05 * CumTrain[t])) # recurrent events (fall-rate goes down with more training)





# Lock the DAG and simulate data
IDAG <- set.DAG(D)
Odat <- sim(DAG = IDAG, n = n_pts, wide = F, rndseed = seed)  # Set wide = F for counting process formatted data


# The data were simulated with each patient followed for t.end,
# the longest possible follow-up time across all patients.
# However, due to staggered recruitment, actual follow-up times vary and are less than t.end.
# filter(t <= real_followtime[ID]) ensures each patient's history is retained only up to their actual follow-up time.
Odat %<>%
  group_by(ID) %>%
  filter(t <= real_followtime[ID]) %>%
  mutate(
    arm = A1[1],               # treatment at day 0
    age = B3[1],               # age at baseline
    base_risk = L1[1],         # Baseline fall-risk

    # Informative dropouts
    ever_dropped = max(Dropout), # indicator if dropout
    temp_drop_time = ifelse(Dropout == 1, t, NA),
    dropout_time = ifelse(!all(is.na(temp_drop_time)), min(temp_drop_time, na.rm = TRUE), NA)
  ) %>%

  # censoring after dropout (ipw)
  filter(is.na(dropout_time) | t <= dropout_time) %>%

  mutate(
    status = ifelse(!is.na(dropout_time) & t == dropout_time, 2, Yobs) #(0=regular censored, 1 =recurrent event, 2=informative dropout)
  ) %>%

  ungroup() %>%
  select(-temp_drop_time) # drop helper variable

# save the raw data
saveRDS(Odat, file=paste0('./results/Odat_1_', seed, ".rds"))

