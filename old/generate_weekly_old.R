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
    drop_intercept = -5,             # Basis-Dropout-Rate
    drop_beta_arm = 0.0,             # Unbalanced Dropout (0 = balanced, >0 = Training group drop out less, 0> = Training group drop out more)
    cens_intercept = -6.9,           # Baseline Censoring rate (exp(-6.9) ~ 0.001) (administrative)
    #cens_age = 0.02,                 # Influence of age on censoring
    train_intercept = -0.5,          # Baseline train rate
    Y_intercept = -3,                # Baseline fall rate
    mortality_level = "low",         # "low" (~0.2%) or "high" (~2%)
    # Parameter for sceanrios
    eff_frailty_y = 0.4,       # Effect frailty on falls
    eff_frailty_cens = 0.6,    # Effect frailty on Dropout & death (informative cens)
    eff_frailty_train = -0.8,  # Effect frailty on training (<0 : train less >0: train more)
    eff_prev_falls = 0.0,      # Effect earlier falls on risk (Event-dependence)
    eff_time = 0.0             # time effect (for time-inhomogen effect)
    
) {
  
 
  # mortality parameter
  # These intercepts yield approx. 0.2% and 2% cumulative mortality at t=52
  death_intercept  <- ifelse(mortality_level == "high", -7.8, -10.1)
  
  
  # DAG initialise
  D <- DAG.empty()
  
  # Baseline (t = 0)
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +    # unobserved frailty (higher frailty -> less training, more falls, more dropout)
    node("age", t = 0,distr = "runif", min = 65, max = 80) +        # age 
    node("sex", t = 0, distr = "rbern", prob = 0.5) +                       # Binary variable for sex, with equal probability for male or female
    
    # Baseline Fall-Risk (historic falls before entry), (not used in the moment)
    node("L1", t = 0, distr = "rpois", 
         lambda = exp(Y_intercept + 0.02 * age[0] + 0.2 * sex[0] + eff_frailty_y * Frailty[0] + 1.5)) +
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
                       plogis(drop_intercept + 0.03 * age[0] + 0.1 * sex[0] + 
                                eff_frailty_cens * Frailty[0] + drop_beta_arm * A1[0]))) +
    # Training participation (driven by age, frailty, and past training)
    node("Train", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 0 | Dropout[t] == 1 | Death[t-1] == 1, 0, 
                       plogis(train_intercept - 0.02 * age[0] + eff_frailty_train * Frailty[0]))) + 
    
    # Cumulative falls up to time t-1
    node("CumYobs", t = 1:t.end, distr = "rconst",
         const = CumYobs[t-1] + Yobs[t-1]) +
    
    # Cumulative number of trainings
    node("CumTrain", t = 1:t.end, distr = "rconst", 
         const = CumTrain[t-1] + Train[t]) +
    
    # Death event (influenced by age, sex, frailty)
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       plogis(death_intercept + 0.05 * age[0] + 0.1 * sex[0] + eff_frailty_cens * Frailty[0]))) +
    
    # Recurrent Event (Falls) (influenced according to scenario)
    # Im moment nur boebachtete Stürze in Studie, nicht vor Beginn!!
    node("Yobs", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1 | Dropout[t] == 1, 0, # No falls if dead or dropped out
                       plogis(Y_intercept + 
                                0.02 * age[0] + 
                                0.2 * sex[0] +
                                true_beta_arm * A1[0] + 
                                eff_frailty_y * Frailty[0] + 
                                eff_prev_falls * CumYobs[t-1] + 
                                eff_time * t)))
  # Simulate data
  IDAG <- set.DAG(D)
  Odat <- sim(DAG = IDAG, n = n_pts, wide = FALSE, rndseed = seed)
  
  Odat %<>%
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = age[1],
      sex = sex[1],
      base_risk = L1[1],
      frailty = Frailty[1],
      
      # Administrative censoring 
      censor_rate = exp(cens_intercept + 0.02 * age), exp_followtime = rexp(1, rate = censor_rate[1]) + 1,
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
      
      cutoff_time = pmin(dropout_time, death_time, real_followtime, na.rm = TRUE)
    ) %>%
    filter(t <= cutoff_time) %>%              
    mutate(
      # 0 = cen, 1 = fall, 2 = Dropout, 3 = death
      status = case_when(
        !is.na(death_time) & t == death_time ~ 3,
        !is.na(dropout_time) & t == dropout_time ~ 2,
        TRUE ~ Yobs
      )
    ) %>%
    ungroup() %>%
    select(-temp_drop_time, -temp_death_time, -cutoff_time)
  
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
                          seed = 123,
                          Y_intercept = -3,
                          mortality_level = "low",
                          eff_frailty_y = 0.4, eff_frailty_cens = 0.6, 
                          eff_prev_falls = 0.0, eff_time = 0.0) {
  
  death_intercept  <- ifelse(mortality_level == "high", -7.8, -10.1)
  
  D <- DAG.empty()
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +
    node("age", t = 0, distr = "runif", min = 65, max = 80) +
    node("sex", t = 0, distr = "rbern", prob = 0.5) +
    node("A1", t = 0, distr = "rconst", const = arm_value) +
    node("CumYobs", t = 0, distr = "rconst", const = 0) +
    node("Yobs", t = 0, distr = "rconst", const = 0) +
    node("Death", t = 0, distr = "rconst", const = 0)
  
  D <- D +
    node("A1", t = 1:t.end, distr = "rconst", const = arm_value) + 
    node("CumYobs", t = 1:t.end, distr = "rconst", const = CumYobs[t-1] + Yobs[t-1]) +
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1, 
                       plogis(death_intercept + 0.05 * age[0] + 0.1 * sex[0] + eff_frailty_cens * Frailty[0]))) +
    node("Yobs", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 0, 
                       plogis(Y_intercept + 0.02 * age[0] + 0.2 * sex[0] + 
                                true_beta_arm * A1[0] + eff_frailty_y * Frailty[0] + 
                                eff_prev_falls * CumYobs[t-1] + eff_time * t)))
  sim(DAG = set.DAG(D), n = n_pts_true, wide = FALSE, rndseed = seed) 
}
############################################################################################
#------- Scenarios & Parallelisation (foreach) ------------------------------------------##
###########################################################################################

n_sim <- 50                    # Number of simulations per scenario (Set to 500 or 1000 later)

# 1. EXPLIZITE DEFINITION DER SZENARIEN (Kausale Struktur)
scenarios_def <- tribble(
  ~scenario_name,        ~eff_frailty_y, ~eff_frailty_cens, ~eff_frailty_train, ~eff_prev_falls, ~mortality_level,
  "1_Null_Scenario",     0.0,            0.0,               0.0,                0.0,             "low",
  "2_Indep_Censoring",   1.2,            0.0,               -0.8,               0.0,             "low",
  "3_Informative_Strong",1.2,            1.2,               -0.8,               0.0,             "low",
  "4_Informative_Weak",  1.2,            1.2,               -0.1,               0.0,             "low",
  "5_Event_Driven",      1.2,            1.2,               -0.8,               0.3,             "low",
  "6_High_Mortality",    1.2,            1.2,               -0.8,               0.0,             "high"
)

#
base_scenarios <- expand_grid(
  n_pts = 2778,
  true_beta_arm = c(0.0, -0.4),          # H0 and H1
  drop_beta_arm = c(0.0, -0.5)           # Balanced vs Unbalanced informative censoring
) %>% 
  cross_join(scenarios_def) %>%          # Verbindet alle Basis-Einstellungen mit den 6 kausalen Szenarien
  mutate(scenario_id = row_number())

# true values
unique_causal_scenarios <- base_scenarios %>% 
  distinct(true_beta_arm, scenario_name, mortality_level, eff_frailty_y, eff_frailty_cens, eff_prev_falls, .keep_all = TRUE)

cat("Set up a parallel backend for true values\n")
no_cores <- parallel::detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

true_estimands_df <- foreach(i = 1:nrow(unique_causal_scenarios), 
                             .combine = bind_rows,
                             .packages = c("simcausal", "dplyr")) %dopar% {
                               
                               p <- unique_causal_scenarios[i, ]
                               
                               # unique seed for each scenario
                               scenario_seed <- 888000 + i
                               
                               # Control group
                               dat_ctrl <- simulate_true(
                                 arm_value = 0, 
                                 true_beta_arm = p$true_beta_arm, 
                                 mortality_level = p$mortality_level, 
                                 seed = scenario_seed,
                                 eff_frailty_y = p$eff_frailty_y, 
                                 eff_frailty_cens = p$eff_frailty_cens, 
                                 eff_prev_falls = p$eff_prev_falls
                               )
                               mcf_ctrl <- sum(dat_ctrl$Yobs) / 50000
                               
                               # Intervention group
                               dat_int <- simulate_true(
                                 arm_value = 1, 
                                 true_beta_arm = p$true_beta_arm, 
                                 mortality_level = p$mortality_level, 
                                 seed = scenario_seed,
                                 eff_frailty_y = p$eff_frailty_y, 
                                 eff_frailty_cens = p$eff_frailty_cens, 
                                 eff_prev_falls = p$eff_prev_falls
                               )
                               mcf_int <- sum(dat_int$Yobs) / 50000
                               
                               # true causal effect
                               rate_ratio <- mcf_int / mcf_ctrl
                               
                               # results
                               data.frame(
                                 Scenario_Name = p$scenario_name,
                                 Mortality = p$mortality_level,
                                 True_Beta_Arm = p$true_beta_arm,
                                 MCF_Control = mcf_ctrl,
                                 MCF_Intervention = mcf_int,
                                 True_RateRatio = rate_ratio,
                                 True_LogRateRatio = log(rate_ratio)
                               )
                             }

stopCluster(cl)

# save true values
if(!dir.exists("./results")) dir.create("./results")
write_csv(true_estimands_df, "./results/00_true_estimands.csv")

cat("True values are simulated and saved.\n")


# Multiply by number of simulation runs
scenarios <- expand_grid(
  base_scenarios,
  sim_run = 1:n_sim
) %>% 
  mutate(seed = 123000 + row_number())

if(!dir.exists("./results")) {
  dir.create("./results")
} else {
  # Clean previous results if needed
  file.remove(list.files("./results", pattern = "\\.rds$", full.names = TRUE))
}

# Setup parallel worker again for the main loop
cat("Set up a parallel backend for main simulation\n")
cl <- makeCluster(no_cores)
registerDoParallel(cl)

cat("Start simulation with foreach...\n")

# Run simulations in parallel
results_log <- foreach(i = 1:nrow(scenarios), 
                       .combine = bind_rows,
                       .packages = c("simcausal", "tidyverse", "magrittr")) %dopar% {
                         
                         p <- scenarios[i, ]      # Current parameters
                         
                         # Generate data (Übergebe nun die neuen kausalen Parameter an die Hauptfunktion)
                         sim_data <- simulate_scenario(
                           seed = p$seed,
                           n_pts = p$n_pts,
                           true_beta_arm = p$true_beta_arm,
                           drop_beta_arm = p$drop_beta_arm,
                           mortality_level = p$mortality_level,
                           eff_frailty_y = p$eff_frailty_y,
                           eff_frailty_cens = p$eff_frailty_cens,
                           eff_frailty_train = p$eff_frailty_train,
                           eff_prev_falls = p$eff_prev_falls
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


# check for mortality and dropout
check_mortality <- Odat %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  summarise(
    Total_Patients = n(),
    Died = sum(status == 3),
    Dropped_Out = sum(status == 2),
    Death_Rate_Percent = round(mean(status == 3) * 100, 2),
    Dropout_Rate_Percent = round(mean(status == 2) * 100, 2)
  )

print(check_mortality)
# ############################################################################################
# #-------Scenarios & Parallelisation (foreach) ---------------------------------------------------------------##
# ###########################################################################################
# 
# n_sim <- 1                    # Number of simulations per scenario(Set to 500 or 1000 later)
# 
# # Grid for scenarios
# base_scenarios <- expand_grid(
#   n_pts = 2778,
#   true_beta_arm = c(0.0, -0.4),          # H0 and H1
#   #scenario_type = 1:5,                   # The 5 causal scenarios
#   scenario_type = 5,
#   mortality_level = c("low", "high"),    # 0.2% vs 2%
#   drop_beta_arm = c(0.0, -1.5)           # Balanced vs Unbalanced informative censoring
# ) |> 
#   mutate(scenario_id = row_number())
# 
# # Simulate true values
# # unique causal scenarios
# unique_causal_scenarios <- base_scenarios %>% 
#   distinct(true_beta_arm, scenario_type, mortality_level, .keep_all = TRUE)
# 
# # use parallel workers
# cat("Set up a parallel backend\n")
# no_cores <- parallel::detectCores() - 1
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# 
# # for each parallel loop for all true scenarios
# true_estimands_df <- foreach(i = 1:nrow(unique_causal_scenarios), 
#                              .combine = bind_rows,
#                              .packages = c("simcausal", "dplyr")) %dopar% {
#                                
#                                params <- unique_causal_scenarios[i, ]
#                                
#                                # unique seed for each scenario
#                                scenario_seed <- 888000 + i
#                                
#                                # controll group
#                                dat_ctrl <- simulate_true(
#                                  arm_value = 0, 
#                                  true_beta_arm = params$true_beta_arm, 
#                                  scenario_type = params$scenario_type, 
#                                  mortality_level = params$mortality_level, 
#                                  seed = scenario_seed
#                                )
#                                mcf_ctrl <- sum(dat_ctrl$Yobs) / 50000
#                                
#                                # Training group
#                                dat_int <- simulate_true(
#                                  arm_value = 1, 
#                                  true_beta_arm = params$true_beta_arm, 
#                                  scenario_type = params$scenario_type, 
#                                  mortality_level = params$mortality_level, 
#                                  seed = scenario_seed
#                                )
#                                mcf_int <- sum(dat_int$Yobs) / 50000
#                                
#                                # true causal effect
#                                rate_ratio <- mcf_int / mcf_ctrl
#                                
#                                # results
#                                data.frame(
#                                  Scenario_Type = params$scenario_type,
#                                  Mortality = params$mortality_level,
#                                  True_Beta_Arm = params$true_beta_arm,
#                                  #Frailty_SD = params$frailty_sd,
#                                  MCF_Control = mcf_ctrl,
#                                  MCF_Intervention = mcf_int,
#                                  True_RateRatio = rate_ratio,
#                                  True_LogRateRatio = log(rate_ratio)
#                                )
#                              }
# 
# # stop cluster
# stopCluster(cl)
# 
# # save true values in csv
# if(!dir.exists("./results")) dir.create("./results")
# write_csv(true_estimands_df, "./results/00_true_estimands.csv")
# 
# cat("True values are simulated and saved.\n")
# 
# 
# #test_Odat <- simulate_scenario(n_pts = 100) # to check error in DAG
# 
# 
# # Multiply by number of simulation runs
# scenarios <- expand_grid(
#   base_scenarios,
#   sim_run = 1:n_sim
# ) |> 
#   mutate(seed = 123000 + row_number())
# 
# if(!dir.exists("./results")) {
#   dir.create("./results")
# } else {
#   # Clean previous results
#   file.remove(list.files("./results", pattern = "\\.rds$", full.names = TRUE))
# }
# 
# # Setup parallel worker
# cat("Set up a parallel backend\n")
# no_cores <- parallel::detectCores() - 1
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# 
# cat("Start simulation with foreach...\n")
# 
# # Run simulations in parallel
# results_log <- foreach(i = 1:nrow(scenarios), 
#                        .combine = bind_rows,
#                        .packages = c("simcausal", "tidyverse", "magrittr")) %dopar% {
#                          
#                          p <- scenarios[i, ]      # Current parameters
#                          
#                          # Generate data
#                          sim_data <- simulate_scenario(
#                            seed = p$seed,
#                            n_pts = p$n_pts,
#                            true_beta_arm = p$true_beta_arm,
#                            #frailty_sd = p$frailty_sd,
#                            drop_beta_arm = p$drop_beta_arm,
#                            scenario_type = p$scenario_type,
#                            mortality_level = p$mortality_level
#                          )
#                          
#                          # Save individual dataset
#                          filename <- sprintf("./results/Odat_scen%02d_run%03d.rds", 
#                                              p$scenario_id, p$sim_run)
#                          saveRDS(sim_data, file = filename)
#                          
#                          # Append path to parameter row for logging
#                          p$file_path <- filename
#                          return(p)
#                        }
# 
# stopCluster(cl)
# 
# # Save logbook
# write_csv(results_log, "./results/00_simulation_log.csv")
# cat("All data records were successfully generated in parallel\n")
# 
# Odat <- readRDS("./results/Odat_scen16_run001.rds")
# 
# # check for mortality and dropout
# check_mortality <- Odat %>%
#   group_by(ID) %>%
#   slice_tail(n = 1) %>% 
#   ungroup() %>%
#   summarise(
#     Total_Patients = n(),
#     Died = sum(status == 3),
#     Dropped_Out = sum(status == 2),
#     Death_Rate_Percent = round(mean(status == 3) * 100, 2),
#     Dropout_Rate_Percent = round(mean(status == 2) * 100, 2)
#   )
# 
# print(check_mortality)
# 
