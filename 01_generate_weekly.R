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
    frailty_sd = 1.0,                # sd of Frailty
    drop_intercept = -5.3,           # Baseline dropout per week for average person (plogis(-5.3) = 0.5% -> ~23% per year)
    drop_beta_arm = 0.0,             # Unbalanced Dropout (0 = balanced, <0 = Training group drop out less, >0 = Training group drop out more)   
    cens_intercept = -6.9,           # Random censoring rate per week (exp(-6.9) ~ 0.001 -> ~5% per year)
    #cens_age = 0.02,                 # Influence of age on censoring
    train_intercept = 0.55,          # Baseline train rate for average person (plogis(0.55) = 63% of weeks)
    contact_intercept = 1.0,         # baseline probability that control group person answer phone (plogis(1) = 73%)
    Y_intercept = -3.5,              # Baseline fall rate per week for average person (exp(-3.5) = 0.03 -> ~1.6 falls per year)
    mortality_level = "low",         # "low" (~0.2%) or "high" (~2%)
    # Parameter for sceanrios
    eff_frailty_y = 0.4,       # Effect frailty on falls
    eff_frailty_cens = 0.6,    # Effect frailty on Dropout & death (informative/dependent cens)
    eff_frailty_train = -0.8,  # Effect frailty on training (<0 : train less >0: train more)
    eff_frailty_contact = -1.5,# Effect frailty on contact (<0: frail patients answer the phone less often) for control arm
    eff_prev_falls = 0.0,      # Effect earlier falls on risk (Event-dependence)
    report_prob = 1.0          # Prob that a fall is reported
    
) {
  
  
  # mortality parameter
  # weakly death rate for 0.2% or 2%.
  # for person with mean age and sex (in Death node specified)
  target_mortality <- ifelse(mortality_level == "high", 0.02, 0.002)
  death_intercept  <- qlogis(1 - (1 - target_mortality)^(1 / t.end))
  
  
  # DAG initialise
  D <- DAG.empty()
  
  # Baseline (t = 0)
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +    # unobserved frailty (higher frailty -> less training, more falls, more dropout)
    node("age", t = 0,distr = "runif", min = 65, max = 80) +                # age 
    node("sex", t = 0, distr = "rbern", prob = 0.5) +                       # Binary variable for sex, with equal probability for male or female
    # Baseline Fall-risk: number of falls in the year before the study (weekly rate * 52) (not used in the moment)
    node("L1", t = 0, distr = "rpois",
         lambda = 52 * exp(Y_intercept + 0.02 * (age[0] - 72.5) + 0.2 * (sex[0] - 0.5) + eff_frailty_y * Frailty[0])) +
    node("A1", t = 0, distr = "rbern", prob = 0.5) +                        # 2 arms, 1 = Training, 0 = no training
    node("Dropout", t = 0, distr = "rbern", prob = 0) +                     # no dropouts at day 0
    node("CumTrain", t = 0, distr = "rconst", const = 0) +                  # start with 0 trainings
    node("Contact", t = 0, distr = "rconst", const = 0) +                   # no call at baseline
    node("CumContact", t = 0, distr = "rconst", const = 0) +                # start with 0 answered calls
    node("CumYfull", t = 0, distr = "rconst", const = 0) +                  # Cumulated true falls at baseline
    node("Yfull", t = 0, distr = "rconst", const = 0) +                     # true falls (also after dropout)
    node("Yobs", t = 0, distr = "rconst", const = 0) +                      # falls observed in study (until dropout)
    node("Yrep", t = 0, distr = "rconst", const = 0) +                      # reported falls
    node("Train", t = 0, distr = "rconst", const = 0) +                  
    node("Death", t = 0, distr = "rconst", const = 0)
  
  # t > 0
  D <- D +
    # Treatment status (constant)
    node("A1", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[t-1] == 1, 1, 0)) + 
    
    # Informative dropout (dependent on intercept, age, sex, frailty (same in both arms), and treatment arm)
    node("Dropout", t = 1:t.end, distr = "rbern",
         prob = ifelse(Dropout[t-1] == 1, 1,
                       plogis(drop_intercept + 0.03 * (age[0] - 72.5) + 0.1 * (sex[0] - 0.5) +   # centered for average age and sex
                                eff_frailty_cens * Frailty[0] +              
                                drop_beta_arm * A1[0]))) +
    # Training participation (driven by age, frailty, sex?, motivation instead of frailty?)
    node("Train", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 0 | Dropout[t] == 1 | Death[t-1] == 1, 0, 
                       plogis(train_intercept - 0.02 * (age[0] - 72.5) + eff_frailty_train * Frailty[0]))) + # centered for average age
    
    # Cumulative true falls up to time t-1
    node("CumYfull", t = 1:t.end, distr = "rconst",
         const = CumYfull[t-1] + Yfull[t-1]) +
    
    # Cumulative number of trainings
    node("CumTrain", t = 1:t.end, distr = "rconst",
         const = CumTrain[t-1] + Train[t]) +
    
    # control arm: patient answer the phone (frailty person answer phone less)
    node("Contact", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 1 | Dropout[t] == 1 | Death[t-1] == 1, 0,
                       plogis(contact_intercept + eff_frailty_contact * Frailty[0]))) +
    
    # cumulative answered calls
    node("CumContact", t = 1:t.end, distr = "rconst",
         const = CumContact[t-1] + Contact[t]) +
    
    # Death event (influenced by age, sex, frailty)
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       plogis(death_intercept + 0.05 * (age[0] - 72.5) + 0.1 * (sex[0] - 0.5) +   # centered for average age and sex
                                eff_frailty_cens * Frailty[0]))) +
    
    # Recurrent Event (Falls) (influenced according to scenario)
    # Yfull = true falls, also after dropout (complete data for oracle M0) only death stops falls
    node("Yfull", t = 1:t.end, distr = "rpois",
         lambda = ifelse(Death[t] == 1, 0,                    # No falls if dead
                         exp(Y_intercept +
                               0.02 * (age[0] - 72.5) +       # centered for average age and sex
                               0.2 * (sex[0] - 0.5) +
                               true_beta_arm * A1[0] +
                               eff_frailty_y * Frailty[0] +
                               eff_prev_falls * CumYfull[t-1]     # later: Falls get more likely with earlier falls
                         ))) +


    # Yobs: falls observed in the study (Yfull until dropout, nothing after)
    node("Yobs", t = 1:t.end, distr = "rconst",
         const = Yfull[t] * (1 - Dropout[t])) +

    # reported falls, fall (Yobs is reported with prob report_prob), in study only Yrep known (possible for later)
    # each of the Yobs falls is reported independently -> Binomial(Yobs, report_prob)
    node("Yrep", t = 1:t.end, distr = "rbinom",
         size = Yobs[t], prob = report_prob)
  # Simulate data
  IDAG <- set.DAG(D)
  Odat <- sim(DAG = IDAG, n = n_pts, wide = FALSE, rndseed = seed)

  # seed for censoring
  set.seed(seed)

  # check for NA
  if (anyNA(Odat$Dropout) | anyNA(Odat$Yobs)) {
    stop("NA in Dropout oder Yobs: check parameter")
  }
  
  Odat %<>%
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = age[1],
      sex = sex[1],
      base_risk = L1[1],
      frailty = Frailty[1],
      
      # Random, non-informative censoring (independent of covariates)
      censor_rate = exp(cens_intercept),                     # ~5% per year     
      exp_followtime = rexp(1, rate = censor_rate[1]) + 1,   # minimum of one week observed
      entry_time = 0,
      max_followtime = t.end - entry_time, 
      real_followtime = round(pmin(exp_followtime, max_followtime)),

      # dropout and death time to compute full data and dropout data (in analyse)
      temp_drop_time = ifelse(Dropout == 1, t, NA),  # if droped out = time, instead NA
      dropout_time = ifelse(!all(is.na(temp_drop_time)), min(temp_drop_time, na.rm = TRUE), NA),

      temp_death_time = ifelse(Death == 1, t, NA),
      death_time = ifelse(!all(is.na(temp_death_time)), min(temp_death_time, na.rm = TRUE), NA),

      ever_died = as.integer(!is.na(death_time)),  

      # end of observation in the study (dropout, death or random censoring)
      cutoff_time = pmin(dropout_time, death_time, real_followtime, na.rm = TRUE)
    ) %>%
    # complete data only ends with death (or t.end)
    filter(is.na(death_time) | t <= death_time) %>%
    mutate(
      # TRUE = row is observed in the study (observed data with filter(observed))
      observed = t <= cutoff_time,

      # observed data: 0 = no event, 1 = fall, 2 = Dropout, 3 = death, NA = not observed
      status = case_when(
        !observed ~ NA_real_,
        !is.na(death_time) & t == death_time ~ 3,
        !is.na(dropout_time) & t == dropout_time ~ 2,
        TRUE ~ as.numeric(Yrep > 0)                # 1 = at least one reported fall (number in Yrep)
      ),
      # dropout only possible if not censored
      ever_dropped = as.integer(any(status == 2, na.rm = TRUE)),

      # complete data (for oracle M0): 0 = no event, 1 = fall, 3 = death
      status_full = case_when(
        !is.na(death_time) & t == death_time ~ 3,
        TRUE ~ as.numeric(Yfull > 0)               # number of falls in Yfull
      )
    ) %>%
    ungroup() %>%
    select(-temp_drop_time, -temp_death_time, -cutoff_time)    # drop helper variables
  
  return(Odat)
}

###########################################################################################
#-------Simulate true values ------------------------------------------------------------##
###### Monte Carlo Simulation of true causal effect ########################################
###########################################################################################
#For true values, simulate a large number of patients without dropouts or censoring for each arm

simulate_true <- function(n_pts_true = 50000, 
                          arm_value = 0, 
                          t.end = 52,
                          true_beta_arm = -0.4, 
                          frailty_sd = 1.0, 
                          seed = 123,
                          Y_intercept = -3.5,
                          mortality_level = "low",
                          eff_frailty_y = 0.4,  
                          eff_frailty_cens = 0.6, 
                          eff_prev_falls = 0.0,
                          with_death = 1) {          # 1 = with dead, 0 = without dead
  
  # mortality percentage
  target_mortality <- ifelse(mortality_level == "high", 0.02, 0.002)
  death_intercept  <- qlogis(1 - (1 - target_mortality)^(1 / t.end))
  
  D <- DAG.empty()
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +
    node("age", t = 0, distr = "runif", min = 65, max = 80) +
    node("sex", t = 0, distr = "rbern", prob = 0.5) +
    node("A1", t = 0, distr = "rconst", const = arm_value) +
    node("CumYfull", t = 0, distr = "rconst", const = 0) +    
    node("Yobs", t = 0, distr = "rconst", const = 0) +
    node("Death", t = 0, distr = "rconst", const = 0)
  
  D <- D +
    node("A1", t = 1:t.end, distr = "rconst", const = arm_value) + 
    node("CumYfull", t = 1:t.end, distr = "rconst", const = CumYfull[t-1] + Yobs[t-1]) +   
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       with_death *                                   
                         plogis(death_intercept + 0.05 * (age[0] - 72.5) + 0.1 * (sex[0] - 0.5) +  
                                  eff_frailty_cens * Frailty[0]))) +
    node("Yobs", t = 1:t.end, distr = "rpois",
         lambda = ifelse(Death[t] == 1, 0,
                         exp(Y_intercept + 0.02 * (age[0] - 72.5) + 0.2 * (sex[0] - 0.5) +
                               true_beta_arm * A1[0] + eff_frailty_y * Frailty[0] +
                               eff_prev_falls * CumYfull[t-1])))
  sim(DAG = set.DAG(D), n = n_pts_true, wide = FALSE, rndseed = seed) 
}
############################################################################################
#------- Scenarios & Parallelisation (foreach) ------------------------------------------##
###########################################################################################

n_sim <- 10                    # Number of simulations per scenario (Set to 500 or 1000 later)

# Scenarios
#  - Frailty effect 0 or informative 
#  - Frailty same impact on death and dropout in both arms
#  - eff_prev_falls: event probability get higher with previous falls (possible for later)
scenarios_def <- tribble(
  ~scenario_name,    ~eff_frailty_y, ~eff_frailty_cens, ~eff_frailty_train, ~eff_frailty_contact, ~eff_prev_falls,
  "1_Null",          0.0,            0.0,               0.0,                0.0,                 0.0,
  "2_Informative",   1.2,            1.2,              -1.5,               -1.5,                 0.0
)


base_scenarios <- expand_grid(
  n_pts = 2778,
  true_beta_arm = c(0.0, -0.4),          # H0 and H1
  drop_beta_arm = c(0.0, -1.0),          # 0 = same Dropout-Rate, -1.0 = control group drops out more often
  mortality_level = c("low", "high"),    # negligible and non-negligible mortality
  report_prob = c(1.0)                   # 1 = all falls reported (maybe change later)
) %>%
  cross_join(scenarios_def) %>%          
  mutate(scenario_id = row_number())

# true values
unique_causal_scenarios <- base_scenarios %>%
  distinct(scenario_name, true_beta_arm, mortality_level, .keep_all = TRUE)

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
                               
                               # without death (for LWWY), Control group
                               dat_ctrl_nd <- simulate_true(
                                 arm_value = 0,
                                 true_beta_arm = p$true_beta_arm,
                                 mortality_level = p$mortality_level,
                                 seed = scenario_seed,
                                 eff_frailty_y = p$eff_frailty_y,
                                 eff_frailty_cens = p$eff_frailty_cens,
                                 eff_prev_falls = p$eff_prev_falls,
                                 with_death = 0
                               )
                               mcf_ctrl_nd <- sum(dat_ctrl_nd$Yobs) / 50000
                               
                               # without death (for LWWY), Intervention group
                               dat_int_nd <- simulate_true(
                                 arm_value = 1,
                                 true_beta_arm = p$true_beta_arm,
                                 mortality_level = p$mortality_level,
                                 seed = scenario_seed,
                                 eff_frailty_y = p$eff_frailty_y,
                                 eff_frailty_cens = p$eff_frailty_cens,
                                 eff_prev_falls = p$eff_prev_falls,
                                 with_death = 0
                               )
                               mcf_int_nd <- sum(dat_int_nd$Yobs) / 50000
                               
                               # true causal effect
                               rate_ratio <- mcf_int / mcf_ctrl
                               
                               # results
                               data.frame(
                                 Scenario_Name = p$scenario_name,
                                 Mortality = p$mortality_level,
                                 True_Beta_Arm = p$true_beta_arm,
                                 MCF_Control = mcf_ctrl,
                                 MCF_Intervention = mcf_int,
                                 MCF_Control_NoDeath = mcf_ctrl_nd,          # without death (for LWYY)
                                 MCF_Intervention_NoDeath = mcf_int_nd,
                                 True_RateRatio = rate_ratio,
                                 True_LogRateRatio = log(rate_ratio),
                                 True_LogRateRatio_NoDeath = log(mcf_int_nd / mcf_ctrl_nd)   
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
                         
                         # Generate data 
                         sim_data <- simulate_scenario(
                           seed = p$seed,
                           n_pts = p$n_pts,
                           true_beta_arm = p$true_beta_arm,
                           drop_beta_arm = p$drop_beta_arm,
                           mortality_level = p$mortality_level,
                           eff_frailty_y = p$eff_frailty_y,
                           eff_frailty_cens = p$eff_frailty_cens,
                           eff_frailty_train = p$eff_frailty_train,
                           eff_frailty_contact = p$eff_frailty_contact,     
                           eff_prev_falls = p$eff_prev_falls,
                           report_prob = p$report_prob                      
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

Odat <- readRDS(results_log$file_path[1])
Odat <- Odat_scen01_run001
check_mortality <- Odat %>%
  group_by(ID) %>%
  mutate(falls_full = sum(Yfull)) %>%       # all falls (complete data)
  filter(observed) %>%                      # only observed data
  mutate(falls_true = sum(Yobs), falls_reported = sum(Yrep, na.rm = TRUE)) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  group_by(arm) %>%
  summarise(
    Prob_Reported = sum(falls_reported) / sum(falls_true),
    Prob_Observed = sum(falls_true) / sum(falls_full),   
    Falls_per_Pt_Full = mean(falls_full),                 
    Total_Patients = n(),
    Died = sum(status == 3, na.rm = TRUE),
    Dropped_Out = sum(status == 2, na.rm = TRUE),
    Death_Rate_Percent = round(mean(status == 3, na.rm = TRUE) * 100, 2),
    Dropout_Rate_Percent = round(mean(status == 2, na.rm = TRUE) * 100, 2)
  )

print(check_mortality)
