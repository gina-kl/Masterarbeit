############################################################################################
######### Data generation #################################################################
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
    drop_intercept = -4.0,           # Basis-Dropout-Rate
    drop_beta_arm = 0.0,             # Unbalanced Dropout (0 = balanced, >0 = Training group drop out less, 0> = Training group drop out more)
    cens_intercept = -6.9,           # Baseline Censoring rate (exp(-6.9) ~ 0.001)
    cens_age = 0.02,                 # Influence of age on censoring
    cens_frailty = 0.5,              # Influence of frailty on censoring
    train_intercept = -0.5,          # Baseline train rate
    death_intercept = -7,            # Baseline death rate
    Y_intercept = -3                 # Baseline fall rate
) {
  # DAG initialise
  D <- DAG.empty()
  
  # Baseline (t = 0)
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +     # unobserved frailty (higher frailty -> less training, more falls, more dropout)
    node("B3", t = 0, distr = "runif", min = 65, max = 95) +               # age
    node("L1", t = 0, distr = "rnorm", mean = 50, sd = 10) +               # Baseline Fall-Risk
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
                       plogis(drop_intercept + 0.03 * B3[0] + 0.6 * Frailty[0] + drop_beta_arm * A1[0]))) + 
    
    # Training participation (driven by age, frailty, and past training)
    node("Train", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 0 | Dropout[t] == 1, 0, 
                       plogis(train_intercept - 0.02 * B3[0] - 0.50 * Frailty[0] + 0.05 * CumTrain[t-1]))) +
    
    # Cumulative falls up to time t-1
    node("CumYobs", t = 1:t.end, distr = "rconst",
         const = CumYobs[t-1] + Yobs[t-1]) +
    
    # Cumulative number of trainings
    node("CumTrain", t = 1:t.end, distr = "rconst", 
         const = CumTrain[t-1] + Train[t]) +
    
    # Death event (influenced by age, frailty, number of falls)
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       plogis(death_intercept + 0.05 * B3[0] + 0.5 * Frailty[0] + 0.3 * CumYobs[t]))) +
    
    # Recurrent event: True treatment effect (true_beta_arm) and frailty (0.4 * Frailty)
    node("Yobs", t = 1:t.end, distr = "rbern",
         prob = plogis(Y_intercept + 0.02 * B3[0] + 0.4 * Frailty[0] + true_beta_arm * A1[0]))
  
  # Simulate data
  IDAG <- set.DAG(D)
  Odat <- sim(DAG = IDAG, n = n_pts, wide = FALSE, rndseed = seed)
  
  Odat %<>%
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = B3[1],
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

# =================================================================
# Scenarios & Parallelisation (foreach)
# =================================================================

n_sim <- 1                    # Number of simulations per scenario

# Grid for scenarios
base_scenarios <- expand_grid(
  n_pts = 2778,
  true_beta_arm = -0.4,       # Ground truth treatment effect
  frailty_sd = c(1.0, 2.0),
  drop_beta_arm = c(0.0, -1.5),
  cens_frailty = c(0.0, 0.5)
) |> 
  mutate(scenario_id = row_number()) # Unique ID 1 to 8

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

# # load packages
# 
# library(simcausal)
# library(tidyverse)
# library(magrittr)
# library(foreach)
# library(doParallel)
# 
# # function to simulate data
# simulate_scenario <- function(
#     seed = 123,
#     n_pts = 2778,                   # number of patients
#     t.end = 52,                     # maximal time (one year in weeks)
#     frailty_sd = 1.0,               # Variance of Frailty
#     drop_intercept = -4.0,          # Basis-Dropout-Rate
#     drop_beta_arm = 0.0 ,           # Unbalanced Dropout (0 = balanced, >0 = Training group drop out less, 0> = Trianing group drop out more)
#     cens_intercept = -6.9,          # Baseline Censoring rate (exp(-6.9) ~ 0.001)
#     cens_age = 0.02,                # Influence of age on censoring
#     cens_frailty = 0.5,             # influence of frailty on censoring
#     train_intercept = -0.5,         # baseline train rate
#     death_intercept = -7,           # baseline death rate
#     Y_intercept = -3                # baseline fall rate
#     
#     ) {
#   # DAG initialise
#   D <- DAG.empty()
#   
#   # Baseline (t = 0)
#   D <- D +
#     node("Frailty", t=0, distr = "rnorm", mean=0, sd=frailty_sd) +    # unobserved frailty (A higher frailty score makes exercise less likely)
#     node("B3", t = 0, distr = "runif", min = 65, max = 95) +          # age
#     node("L1", t = 0, distr = "rnorm", mean = 50, sd = 10) +          # Baseline Fall-Risk
#     node("A1", t = 0, distr = "rbern", prob = 0.5) +                  # 2 arms, 1 = Training, 0=no training
#     node("Dropout", t = 0, distr = "rbern", prob = 0) +               # no dropouts at day 1
#     node("CumTrain", t = 0, distr = "rconst", const = 0) +            # start with 0 trainings
#     node("CumYobs", t = 0, distr = "rconst", const = 0) +             # Cumulated trainings at basline
#     node("Yobs", t = 0, distr = "rconst", const = 0) +                
#     node("Train", t=0, distr = "rconst", const =0 ) +                 
#     node("Death", t=0, distr = "rconst", const =0 )
#   
#   # t > 0
#   D <- D +
#     # Treatment status (constant)
#     node("A1", t = 1:t.end, distr = "rbern",
#          prob = ifelse(A1[t-1] == 1, 1, 0)) +            
#     
#     # informative dropout (dependent on baseline age, frailty, and treatment arm)
#     node("Dropout", t = 1:t.end, distr = "rbern",
#          prob = ifelse(Dropout[t-1] == 1, 1,
#                        plogis(drop_intercept + 0.03 * B3[0] + 0.6 * Frailty[0] + drop_beta_arm * A1[0]))) + 
#     
#     # Training participation (0 for control arm or post-dropout; driven by age, frailty, and past training)
#     node("Train", t = 1:t.end, distr = "rbern",
#          prob = ifelse(A1[0] == 0 | Dropout[t] == 1, 0, 
#                        plogis(train_intercept - 0.02 * B3[0] - 0.50 * Frailty[0] + 0.05 * CumTrain[t-1]))) +
#     # Cumulative falls up to time t-1
#     node("CumYobs", t = 1:t.end, distr = "rconst",
#          const = CumYobs[t-1] + Yobs[t-1]) +
#     
#     # Cumulative number of trainings
#     node("CumTrain", t = 1:t.end, distr = "rconst", 
#          const = CumTrain[t-1] + Train[t]) +
#     
#     # Death event (influenced by age, frailty, number of falls)
#     node("Death", t = 1:t.end, distr = "rbern",
#          prob = ifelse(Death[t-1] == 1, 1,
#                        plogis(death_intercept + 0.05 * B3[0] + 0.5 * Frailty[0] + 0.3 * CumYobs[t]))) +
#     
#     # Recurrent event (risk reduced by cumulative training)
#     node("Yobs", t = 1:t.end, distr = "rbern",
#          prob = plogis(Y_intercept + 0.02 * B3[0] + 0.4 * Frailty[0] - 0.05 * CumTrain[t]))
#   
#   # simulate data
#   IDAG <- set.DAG(D)
#   Odat <- sim(DAG = IDAG, n = n_pts, wide = F, rndseed = seed)  # Set wide = F for counting process formatted data
#   
#     
#   Odat %<>%
#     group_by(ID) %>%
#     mutate(
#       arm = A1[1],
#       age = B3[1],
#       base_risk = L1[1],
#       frailty = Frailty[1],
#       
#       # censoring time per patient (influenced by age, frailty)
#       censor_rate = exp(cens_intercept + cens_age * age + cens_frailty * frailty),
#       
#       exp_followtime = rexp(1, rate = censor_rate[1]) + 1,            # real follow up time (~exp)
#       entry_time = 0,                                                 # entry time (all starts at time 0?)
#       max_followtime = t.end - entry_time,                            
#       real_followtime = round(pmin(exp_followtime, max_followtime))
#     ) %>%
#     filter(t <= real_followtime) %>%                                  # filter(t <= real_followtime[ID]) ensures each patient's history is retained only up to their actual follow-up time.
#     mutate(
#       ever_dropped = max(Dropout),             
#       ever_died = max(Death),                  
#       
#       temp_drop_time = ifelse(Dropout == 1, t, NA),
#       dropout_time = ifelse(!all(is.na(temp_drop_time)), min(temp_drop_time, na.rm = TRUE), NA),
#       
#       temp_death_time = ifelse(Death == 1, t, NA),
#       death_time = ifelse(!all(is.na(temp_death_time)), min(temp_death_time, na.rm = TRUE), NA),
#       
#       end_of_obs = pmin(dropout_time, death_time, na.rm = TRUE)
#     ) %>%
#     filter(is.na(dropout_time) | t <= dropout_time) %>%
#     mutate(
#       # 0 = cen, 1 = fall, 2 = Dropout, 3 = death
#       status = case_when(
#         !is.na(death_time) & t == death_time ~ 3,
#         !is.na(dropout_time) & t == dropout_time ~ 2,
#         TRUE ~ Yobs
#       )
#     ) %>%
#     ungroup() %>%
#     select(-temp_drop_time, -temp_death_time, -end_of_obs)
#   
#   return(Odat)
# }
# 
# 
# # =================================================================
# # Scenarios & Parallelisation (foreach)
# # =================================================================
# 
# n_sim <- 1                    # number of simulation for each scenario
# 
# # grid for scenarios
# base_scenarios <- expand_grid(
#   n_pts = 2778,
#   frailty_sd = c(1.0, 2.0),
#   drop_beta_arm = c(0.0, -1.5),
#   cens_frailty = c(0.0, 0.5)
# ) |> 
#   mutate(scenario_id = row_number()) # ID for each scenario
# 
# # simulate number of sims for each senario
# scenarios <- expand_grid(
#   base_scenarios,
#   sim_run = 1:n_sim
# ) |> 
#   mutate(seed = 123000 + row_number())
# 
# if(!dir.exists("./results")) {
#   dir.create("./results")
# } else {
#   # delete all old files
#   file.remove(list.files("./results", pattern = "\\.rds$", full.names = TRUE))
# }
# # Setup parallel worker
# cat("Set up a parallel backend\n")
# no_cores <- parallel::detectCores() - 1
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# 
# cat("Start simulation with foreach...\n")
# 
# # foreach iterate over tibble 
# results_log <- foreach(i = 1:nrow(scenarios), 
#                    .combine = bind_rows,
#                    .packages = c("simcausal", "tidyverse", "magrittr")) %dopar% {
#                      
#                      p <- scenarios[i, ]      # extract parameter of current row
#                      
#                      # generate data
#                      sim_data <- simulate_scenario(
#                        seed = p$seed,
#                        n_pts = p$n_pts,
#                        frailty_sd = p$frailty_sd,
#                        drop_beta_arm = p$drop_beta_arm,
#                        cens_frailty = p$cens_frailty
#                      )
#                      
#                      # save
#                      filename <- sprintf("./results/Odat_scen%02d_run%03d.rds", 
#                                          p$scenario_id, p$sim_run)
#                      saveRDS(sim_data, file = filename)
#                      
#                      p$file_path <- filename
#                      return(p)
#                    }
# 
# # logbook with scenarios
# write_csv(results_log, "./results/00_simulation_log.csv")
# 
# stopCluster(cl)
# cat("All data records were successfully generated in parallel\n")
# 
