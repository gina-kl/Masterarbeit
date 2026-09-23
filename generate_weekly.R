############################################################################################
#------- Data generation ---------------------------------------------------------------##
###########################################################################################
# Load packages
library(simcausal)
library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)

# Function to simulate data
simulate_scenario <- function(
    seed = 123,
    n_pts = 2778,                    # number of patients
    t.end = 52,                      # maximal time (one year in weeks)
    true_beta_arm = -0.4,            # TRUE treatment effect on recurrent falls
    frailty_sd = 1.0,                # Variance of Frailty
    drop_intercept = -5,           # Basis-Dropout-Rate
    drop_beta_arm = 0.0,             # Unbalanced Dropout (0 = balanced, >0 = Training group drop out less, 0> = Training group drop out more)
    cens_intercept = -6.9,           # Baseline Censoring rate (exp(-6.9) ~ 0.001)
    cens_age = 0.02,                 # Influence of age on censoring
    cens_frailty = 0.5,              # Influence of frailty on censoring
    train_intercept = -0.5,          # Baseline train rate
    death_intercept = -9,           # Baseline death rate (low)
    Y_intercept = -3                 # Baseline fall rate
) {
  # DAG initialise
  D <- DAG.empty()
  
  # Baseline (t = 0)
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +    # unobserved frailty (higher frailty -> less training, more falls, more dropout)
    node("age", t = 0,distr = "runif", min = 65, max = 80) +        # age 
    node("sex", t = 0, distr = "rbern", prob = 0.5) +                       # Binary variable for sex, with equal probability for male or female
    node("L1", t = 0, distr = "rnorm", mean = 50, sd = 10) +                # Baseline Fall-Risk
    node("A1", t = 0, distr = "rbern", prob = 0.5) +                        # 2 arms, 1 = Training, 0 = no training
    node("Dropout", t = 0, distr = "rbern", prob = 0) +                     # no dropouts at day 0
    node("CumTrain", t = 0, distr = "rconst", const = 0) +                  # start with 0 trainings
    node("CumYobs", t = 0, distr = "rconst", const = 0) +                   # Cumulated falls at baseline
    node("Yobs", t = 0, distr = "rconst", const = 0) + 
    node("Train", t = 0, distr = "rconst", const = 0) + 
    node("Death", t = 0, distr = "rconst", const = 0)
  
  # t > 0
  D <- D +
    # Treatment status (constant)
    node("A1", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[t-1] == 1, 1, 0)) + 
    
    # Informative dropout (dependent on baseline age, frailty, and treatment arm)
    node("Dropout", t = 1:t.end, distr = "rbern",
         prob = ifelse(Dropout[t-1] == 1, 1,
                       plogis(drop_intercept + 0.03 * age[0] + 0.6 * Frailty[0] + drop_beta_arm * A1[0]))) + 
    
    # Training participation (driven by age, frailty, and past training)
    node("Train", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 0 | Dropout[t] == 1, 0, 
                       plogis(train_intercept - 0.02 * age[0] - 0.50 * Frailty[0]))) +
    
    # Cumulative falls up to time t-1
    node("CumYobs", t = 1:t.end, distr = "rconst",
         const = CumYobs[t-1] + Yobs[t-1]) +
    
    # Cumulative number of trainings
    node("CumTrain", t = 1:t.end, distr = "rconst", 
         const = CumTrain[t-1] + Train[t]) +
    
    # Death event (influenced by age, frailty, number of falls)
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       plogis(death_intercept + 0.05 * age[0] + 0.5 * Frailty[0]))) +
    
    # Recurrent event: True treatment effect (true_beta_arm) and frailty (0.4 * Frailty)
    node("Yobs", t = 1:t.end, distr = "rbern",
         prob = plogis(Y_intercept + 0.02 * age[0] + 0.4 * Frailty[0] + true_beta_arm * A1[0]))
  
  # Simulate data
  IDAG <- set.DAG(D)
  Odat <- sim(DAG = IDAG, n = n_pts, wide = FALSE, rndseed = seed)
  
  Odat %<>%
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = age[1],
      base_risk = L1[1],
      frailty = Frailty[1],
      
      # Censoring time per patient (influenced by age, frailty)
      censor_rate = exp(cens_intercept + cens_age * age + cens_frailty * frailty),
      exp_followtime = rexp(1, rate = censor_rate[1]) + 1,
      entry_time = 0,
      max_followtime = t.end - entry_time, 
      real_followtime = round(pmin(exp_followtime, max_followtime))
    ) %>%
    filter(t <= real_followtime) %>%
    mutate(
      ever_dropped = max(Dropout), 
      ever_died = max(Death), 
      
      temp_drop_time = ifelse(Dropout == 1, t, NA),
      dropout_time = ifelse(!all(is.na(temp_drop_time)), min(temp_drop_time, na.rm = TRUE), NA),
      
      temp_death_time = ifelse(Death == 1, t, NA),
      death_time = ifelse(!all(is.na(temp_death_time)), min(temp_death_time, na.rm = TRUE), NA),
      
      end_of_obs = pmin(dropout_time, death_time, na.rm = TRUE)
    ) %>%
    filter(is.na(dropout_time) | t <= dropout_time) %>%
    mutate(
      # 0 = cen, 1 = fall, 2 = Dropout, 3 = death
      status = case_when(
        !is.na(death_time) & t == death_time ~ 3,
        !is.na(dropout_time) & t == dropout_time ~ 2,
        TRUE ~ Yobs
      )
    ) %>%
    ungroup() %>%
    select(-temp_drop_time, -temp_death_time, -end_of_obs)
  
  return(Odat)
}

############################################################################################
#-------Simulate true values ---------------------------------------------------------------##
###########################################################################################
#For true values, simulate a large number of patients without dropouts or censoring for each arm

simulate_true <- function(n_pts_true = 50000, 
                          arm_value = 0, 
                          t.end = 52,
                          true_beta_arm = -0.4, 
                          frailty_sd = 1.0, 
                          seed=123,
                          death_intercept = -10,
                          Y_intercept = -3  ) {
  D <- DAG.empty()
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +
    node("age", t = 0, distr = "runif", min = 65, max = 95) +
    node("sex", t = 0, distr = "rbern", prob = 0.5) +
    node("A1", t = 0, distr = "rconst", const = arm_value) +
    node("CumYobs", t = 0, distr = "rconst", const = 0) +
    node("Yobs", t = 0, distr = "rconst", const = 0) +
    node("Death", t = 0, distr = "rconst", const = 0)
  
  D <- D +
    node("A1", t = 1:t.end, distr = "rconst", const = arm_value) + 
    node("CumYobs", t = 1:t.end, distr = "rconst", const = CumYobs[t-1] + Yobs[t-1]) +
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1, plogis(death_intercept + 0.05 * age[0] + 0.5 * Frailty[0]))) +
    node("Yobs", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 0, plogis(Y_intercept + 0.02 * age[0] + 0.4 * Frailty[0] + true_beta_arm * A1[0])))
  
  sim(DAG = set.DAG(D), n = n_pts_true, wide = FALSE, rndseed = seed) 
}


############################################################################################
#-------Scenarios & Parallelisation (foreach) ---------------------------------------------------------------##
###########################################################################################

n_sim <- 1                    # Number of simulations per scenario

# Grid for scenarios
base_scenarios <- expand_grid(
  n_pts = 2778,
  true_beta_arm = c(0.0, -0.4), # 0=same event-rate (H0), -0.4 = different event-rate (H1)
  frailty_sd = c(0.5, 1.5),
  drop_beta_arm = c(0.0, -1.5),
  cens_frailty = c(0.0, 0.5)
) |> 
  mutate(scenario_id = row_number()) # Unique ID 1 to 8

# Simulate true values
# unique causal scenarios
unique_causal_scenarios <- base_scenarios %>% 
  distinct(true_beta_arm, frailty_sd, .keep_all = TRUE)

# use parallel workers
no_cores <- parallel::detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# for each parallel loop for all true sceanrios
true_estimands_df <- foreach(i = 1:nrow(unique_causal_scenarios), 
                             .combine = bind_rows,
                             .packages = c("simcausal", "dplyr")) %dopar% {
                               
                               params <- unique_causal_scenarios[i, ]
                               
                               # unique seed for each scenario
                               scenario_seed <- 888000 + i
                               
                               # controll group
                               dat_ctrl <- simulate_true(
                                 n_pts_true = 50000, 
                                 arm_value = 0, 
                                 true_beta_arm = params$true_beta_arm, 
                                 frailty_sd = params$frailty_sd,
                                 seed = scenario_seed 
                               )
                               mcf_ctrl <- sum(dat_ctrl$Yobs) / 50000
                               
                               # Training group
                               dat_int <- simulate_true(
                                 n_pts_true = 50000, 
                                 arm_value = 1, 
                                 true_beta_arm = params$true_beta_arm, 
                                 frailty_sd = params$frailty_sd,
                                 seed = scenario_seed 
                               )
                               mcf_int <- sum(dat_int$Yobs) / 50000
                               
                               # true causal effect
                               rate_ratio <- mcf_int / mcf_ctrl
                               
                               # results
                               data.frame(
                                 scenario_id = params$scenario_id,
                                 True_Beta_Arm = params$true_beta_arm,
                                 Frailty_SD = params$frailty_sd,
                                 MCF_Control = mcf_ctrl,
                                 MCF_Intervention = mcf_int,
                                 True_RateRatio = rate_ratio,
                                 True_LogRateRatio = log(rate_ratio)
                               )
                             }

# stop cluster
stopCluster(cl)

# save true values in csv
if(!dir.exists("./results")) dir.create("./results")
write_csv(true_estimands_df, "./results/00_true_estimands.csv")

cat("True values are simulated and saved.\n")


#test_Odat <- simulate_scenario(n_pts = 100) # to check error in DAG


# Multiply by number of simulation runs
scenarios <- expand_grid(
  base_scenarios,
  sim_run = 1:n_sim
) |> 
  mutate(seed = 123000 + row_number())

if(!dir.exists("./results")) {
  dir.create("./results")
} else {
  # Clean previous results
  file.remove(list.files("./results", pattern = "\\.rds$", full.names = TRUE))
}

# Setup parallel worker
cat("Set up a parallel backend\n")
no_cores <- parallel::detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

cat("Start simulation with foreach...\n")

# Run simulations in parallel
results_log <- foreach(i = 1:nrow(scenarios), 
                       .combine = bind_rows,
                       .packages = c("simcausal", "tidyverse", "magrittr")) %dopar% {
                         
                         p <- scenarios[i, ]      # Current parameters
                         
                         # Generate data
                         sim_data <- simulate_scenario(
                           seed = p$seed,
                           n_pts = p$n_pts,
                           true_beta_arm = p$true_beta_arm,
                           frailty_sd = p$frailty_sd,
                           drop_beta_arm = p$drop_beta_arm,
                           cens_frailty = p$cens_frailty
                         )
                         
                         # Save individual dataset
                         filename <- sprintf("./results/Odat_scen%02d_run%03d.rds", 
                                             p$scenario_id, p$sim_run)
                         saveRDS(sim_data, file = filename)
                         
                         # Append path to parameter row for logging
                         p$file_path <- filename
                         return(p)
                       }

stopCluster(cl)

# Save logbook
write_csv(results_log, "./results/00_simulation_log.csv")
cat("All data records were successfully generated in parallel\n")

Odat <- readRDS("./results/Odat_scen16_run001.rds")

# check for mortality and dropout
check_mortality <- Odat %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%  # last observation week
  ungroup() %>%
  summarise(
    Total_Patients = n(),
    Died = sum(status == 3),
    Dropped_Out = sum(status == 2),
    Death_Rate_Percent = round(mean(status == 3) * 100, 2),
    Dropout_Rate_Percent = round(mean(status == 2) * 100, 2)
  )

print(check_mortality)

