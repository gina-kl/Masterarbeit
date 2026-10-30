
##########################################################################################
#------- Data generation ---------------------------------------------------------------##
##########################################################################################
source("cluster_setup.R")   # run_on_cluster, n_cores_use, results_dir, part
source("functions.R")       # analysis functions (tie_offset, save_rds_atomic)
source("functions_dgp.R")   # data generating process: DAG, dropout, calibration, true values
source("GlobalVars.R")

# Load packages, Rscript on the cluster starts without any packages
library(simcausal)
# single packages instead of tidyverse (tidyverse cannot be installed on the cluster)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(magrittr)
library(foreach)
library(doParallel)
library(survival)   # true values: LWYY
library(mets)       # true values: Ghosh-Lin


#------------------- Scenarios -----------------------------------------------------------##
# scenario definitions, effects, blocks, scenario_id and causal_id: scenarios.R
source("scenarios.R")

# intercepts per scenario (calibration, functions_dgp.R), targets from GlobalVars.R
base_scenarios <- calibrate_scenarios(base_scenarios)

message("Scenario blocks: ", paste(scenario_blocks, collapse = ", "),
        " | Number of Scenarios (base_scenarios): ", nrow(base_scenarios))   

# for run on cluster: which part should run
run_true <- part %in% c("all", "true") # run true values
run_data <- part %in% c("all", "data") # run simulated data
if (!run_true && !run_data) stop("01_generate: argument must be 'true' or 'data', not: ", part)


#------------------- True values ---------------------------------------------------------##


# The true values do not depend on the dropout mechanism (Yfull continues after dropout), they are computed once per causal scenario
# Monte Carlo in blocks of n_chunk_true patients per arm (run parallel on cluster)
#   - one parallel task = one causal scenario x one block -> memory per worker depends only on n_chunk_true
#   - per block only small summaries are kept (falls per arm and quarter, deaths, model coefficients)
#   - mean numbers of falls: pooled over all blocks, LWYY / GL: mean of the block estimates
#     (bias of a block estimate is O(1 / n_chunk_true), negligible)
#   - Monte Carlo SE of each true value = sd over the blocks / sqrt(n_blocks)
#   - seed depends only on scenario and block -> a larger n_true keeps the first blocks and adds new ones



if (run_true) {   # cluster
  # unique causal scenarios
  # causal_id from scenarios.R (complete grid) -> same id and seeds whichever blocks are run
  unique_causal_scenarios <- base_scenarios %>%
    distinct(causal_id, .keep_all = TRUE) %>%
    dplyr::select(causal_id, effect, all_of(causal_vars), Y_intercept, death_intercept, train_intercept)   # intercepts follow from causal_vars
  
  if (n_true %% n_chunk_true != 0) stop("n_true must be a multiple of n_chunk_true")
  n_blocks <- n_true / n_chunk_true
  
  true_tasks <- expand_grid(causal_id = unique_causal_scenarios$causal_id,
                            block     = seq_len(n_blocks))
  # one task = true_block() (functions_dgp.R) for one causal scenario and one block

  cat("Set up a parallel backend for true values\n")
  
  n_workers_true <- if (run_on_cluster) n_cores_use else min(n_cores_use, n_cores_true)
  
  cat("Workers for true values:", n_workers_true, "| tasks:", nrow(true_tasks), "\n")
  cl <- makeCluster(n_workers_true)   # each worker needs memory for one block (see n_chunk_true)
  registerDoParallel(cl)
  
  
  true_blocks <- foreach(j = seq_len(nrow(true_tasks)),
                         .combine = bind_rows,
                         .packages = c("simcausal", "dplyr", "tidyr", "survival", "mets"),
                         .export = c("COEF", "sim_dag", "simulate_true", "prep_data_true", "tie_offset", "true_block",
                                     "model_mcf", "unique_causal_scenarios", "true_tasks", "n_chunk_true")) %dopar% {
                                       i <- true_tasks$causal_id[j]  # scenario
                                       k <- true_tasks$block[j]     # block
                                       # unique per scenario and block, separate from the data seeds (1e6 * causal_id + run, dropout + 5e8)
                                       # row by causal_id (not by position, the ids are not consecutive within a block)
                                       true_block(unique_causal_scenarios[unique_causal_scenarios$causal_id == i, ], k,
                                                  n_chunk_true, seed = 1e9 + 1e5 * i + k)
                                     }
  
  stopCluster(cl)
  
  # cluster: results_dir instead "./results" 
  saveRDS(true_blocks, file.path(results_dir, "00_true_blocks.rds"))    # block results: check MC error, add blocks later
  
  # Monte Carlo SE over the blocks (SE=s/sqrt(n))
  mcse <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ 
  
  true_estimands_df <- true_blocks %>%
    group_by(causal_id) %>%
    summarise(
      N_True = n() * n_chunk_true,                                   # patients per arm
      # Monte Carlo SE from the variation between the blocks (before the counts are summed up)
      MCSE_LogRateRatio         = mcse(log(y_int / y_ctrl)),
      MCSE_LogRateRatio_NoDeath = mcse(log(y_int_nd / y_ctrl_nd)),
      MCSE_Beta_LWYY            = mcse(beta_lwyy),
      MCSE_Beta_GL              = mcse(beta_gl),
      MCSE_MCF_LWYY_Control      = mcse(mcf_lwyy_ctrl),
      MCSE_MCF_LWYY_Intervention = mcse(mcf_lwyy_int),
      MCSE_LogRR_Q1             = mcse(log(y_int_Q1 / y_ctrl_Q1)),
      MCSE_LogRR_Q2             = mcse(log(y_int_Q2 / y_ctrl_Q2)),
      MCSE_LogRR_Q3             = mcse(log(y_int_Q3 / y_ctrl_Q3)),
      MCSE_LogRR_Q4             = mcse(log(y_int_Q4 / y_ctrl_Q4)),
      True_Beta_LWYY            = mean(beta_lwyy),                   # target of the adjusted LWYY
      True_Beta_GL              = mean(beta_gl),                     # target of the adjusted Ghosh-Lin
      MCF_LWYY_Control          = mean(mcf_lwyy_ctrl),               # target of the model-based LWYY falls per arm
      MCF_LWYY_Intervention     = mean(mcf_lwyy_int),                #   (LWYY on complete data, as True_Beta_LWYY)
      across(c(starts_with("y_"), starts_with("d_")), sum),          # pooled counts over all blocks
      .groups = "drop"
    ) %>%
    mutate(
      MCF_Control               = y_ctrl    / N_True,                # mean number of falls until t.end per arm
      MCF_Intervention          = y_int     / N_True,
      MCF_Control_NoDeath       = y_ctrl_nd / N_True,
      MCF_Intervention_NoDeath  = y_int_nd  / N_True,
      True_LogRateRatio         = log(MCF_Intervention / MCF_Control),                 # ratio of mean falls (nonparametric GL, marginal)
      True_LogRateRatio_NoDeath = log(MCF_Intervention_NoDeath / MCF_Control_NoDeath), # hypothetical, no death
      LogRR_Q1                  = log(y_int_Q1 / y_ctrl_Q1),         # time-varying effect per quarter
      LogRR_Q2                  = log(y_int_Q2 / y_ctrl_Q2),
      LogRR_Q3                  = log(y_int_Q3 / y_ctrl_Q3),
      LogRR_Q4                  = log(y_int_Q4 / y_ctrl_Q4),
      Mort_Control              = d_ctrl / N_True,                   # achieved mortality per arm
      Mort_Intervention         = d_int  / N_True
    ) %>%
    left_join(unique_causal_scenarios %>% dplyr::select(causal_id, effect, all_of(causal_vars)),   
              by = "causal_id") %>%
    # same truth for all scenarios with the same causal parameters (e.g. 2, 3, 4 only differ in dropout)
    inner_join(base_scenarios %>%
                 distinct(across(all_of(causal_vars)), scenario_name),
               by = causal_vars, relationship = "one-to-many") %>%
    mutate(Scenario_Name = scenario_name, Mortality = mortality_level, True_Beta_Arm = true_beta_arm,
           Effect = effect, Max_Frailty_Red = max_frailty_red) %>%                                     
    
    # strict H0 (no direct effect and no effect via training): both arms have the same distribution
    # -> all effects are exactly 0, levels per arm are the same (mean of both arms = less MC error)
    # the MCSE columns still show the Monte Carlo noise of the simulated values
    mutate(
      h0 = True_Beta_Arm == 0 & train_frailty_eff == 0,
      across(c(True_LogRateRatio, True_LogRateRatio_NoDeath, True_Beta_LWYY, True_Beta_GL, starts_with("LogRR_Q")),
             ~ ifelse(h0, 0, .x)),
      MCF_Control              = ifelse(h0, (MCF_Control + MCF_Intervention) / 2, MCF_Control),
      MCF_Intervention         = ifelse(h0, MCF_Control, MCF_Intervention),
      MCF_Control_NoDeath      = ifelse(h0, (MCF_Control_NoDeath + MCF_Intervention_NoDeath) / 2, MCF_Control_NoDeath),
      MCF_Intervention_NoDeath = ifelse(h0, MCF_Control_NoDeath, MCF_Intervention_NoDeath),
      MCF_LWYY_Control         = ifelse(h0, (MCF_LWYY_Control + MCF_LWYY_Intervention) / 2, MCF_LWYY_Control),
      MCF_LWYY_Intervention    = ifelse(h0, MCF_LWYY_Control, MCF_LWYY_Intervention),
      True_RateRatio           = exp(True_LogRateRatio)
    ) %>%
    dplyr::select(Scenario_Name, Mortality, True_Beta_Arm, Effect, Max_Frailty_Red,  
                  MCF_Control, MCF_Intervention, MCF_Control_NoDeath, MCF_Intervention_NoDeath,
                  MCF_LWYY_Control, MCF_LWYY_Intervention,
                  True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath,
                  True_Beta_LWYY, True_Beta_GL, Mort_Control, Mort_Intervention,
                  LogRR_Q1, LogRR_Q2, LogRR_Q3, LogRR_Q4,
                  N_True, starts_with("MCSE_"))
  
  write_csv(true_estimands_df, file.path(results_dir, "00_true_estimands.csv"))

  # scenario list for 03
  write_csv(base_scenarios, file.path(results_dir, "00_base_scenarios.csv"))

  cat("True values are simulated and saved.\n")
  
}   # ---- end of "true" -------------------------------------------------




#------------------- Simulate all scenarios ----------------------------------------------##



if (run_data) {   # cluster ---- start "data" -------------------------------------
  
  # local all scenarios, on cluster only scenario of array_task
  # array task k = k-th scenario of the chosen blocks (scenario_id itself is not consecutive)
  if (run_on_cluster && !(array_id %in% seq_len(nrow(base_scenarios)))) {
    stop("Array-Task ", array_id, " don't exist: only 1 - ", nrow(base_scenarios), " scenarios in block(s) ",
         paste(scenario_blocks, collapse = ", "))
  }
  scenario_ids <- if (run_on_cluster) base_scenarios$scenario_id[array_id] else base_scenarios$scenario_id
  
  # Multiply by number of simulation runs
  scenarios <- expand_grid(
    base_scenarios %>% filter(scenario_id %in% scenario_ids),   # cluster
    sim_run = 1:n_sim
  ) %>%
    mutate(seed      = 1e6 * causal_id + sim_run,   # same complete data for all scenarios with the same causal parameters
           seed_drop = seed + 5e8,                  # own stream for dropout (same uniforms for all dropout mechanisms)
           file_path = file.path(results_dir,                    
                                 sprintf("Odat_scen%02d_run%04d.rds", scenario_id, sim_run)))
  
  # all data sets of the scenario(s) are always generated again, existing files are overwritten
  # (cluster: after a restart the whole scenario is computed again, about 12 min)
  if (!run_on_cluster) {
    #local: remove all old data
    file.remove(list.files(results_dir, pattern = "^Odat_.*\\.rds$", full.names = TRUE))
  }

  # Setup parallel worker again for the main loop
  cat("Set up a parallel backend for main simulation\n")
  cl <- makeCluster(n_cores_use)  
  registerDoParallel(cl)

  cat("Start simulation with foreach:", nrow(scenarios), "data sets on", n_cores_use, "cores\n")
  t0 <- Sys.time()

  # Run simulations in parallel
  results_log <- foreach(i = 1:nrow(scenarios),
                         .combine = bind_rows,
                         .packages = c("simcausal", "dplyr", "tidyr", "tibble", "magrittr"),   # instead of tidyverse
                         .export = c("COEF", "sim_dag", "add_dropout", "simulate_scenario", "save_rds_atomic")) %dopar% {

                           p <- scenarios[i, ]      # Current parameters

                           # Generate data
                           sim_data <- simulate_scenario(
                             seed = p$seed,
                             seed_drop = p$seed_drop,
                             n_pts = p$n_pts,
                             true_beta_arm = p$true_beta_arm,
                             drop_beta_arm = p$drop_beta_arm,
                             Y_intercept = p$Y_intercept,
                             drop_intercept = p$drop_intercept,
                             death_intercept = p$death_intercept,
                             eff_frailty_y = p$eff_frailty_y,
                             eff_frailty_death = p$eff_frailty_death,
                             eff_frailty_drop = p$eff_frailty_drop,
                             eff_fall_drop = p$eff_fall_drop,
                             max_frailty_red = p$max_frailty_red,
                             eff_frailty_train = p$eff_frailty_train,
                             train_intercept = p$train_intercept,
                             train_frailty_eff = p$train_frailty_eff
                           )

                           # Save individual dataset (path set above: results_dir/Odat_scenXX_runYYYY.rds)
                           # atomic -> no broken file if the job is killed while saving
                           save_rds_atomic(sim_data, p$file_path)

                           return(p)
                         }

  stopCluster(cl)
  cat("Data sets generated in", format(round(Sys.time() - t0, 1)), "\n")

  # Save logbook (one row per data set)
  # local: one log for everything, cluster: one log per scenario (no simultaneous writing)
  log_file <- if (run_on_cluster) {
    file.path(results_dir, sprintf("00_simulation_log_scen%02d.csv", scenario_ids))   
  } else {
    file.path(results_dir, "00_simulation_log.csv")
  }
  write_csv(results_log, log_file)
  cat(nrow(results_log), "data sets in log:", log_file, "\n")
  
}  