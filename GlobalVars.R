############################################################################
### GlobalVars #############################################################
############################################################################

if (!exists("run_on_cluster")) run_on_cluster <- FALSE


# --- Variables for DGP -------------------------------
# Simulations
n_sim <- if(run_on_cluster) 1000 else 2                 
n_true  <- if(run_on_cluster) 4e6 else 20000  # patients per arm for Monte Carlo for true values
# 4e6 (400 blocks): MCSE of the true values ~0.001 (1e6: ~0.002 = 50-80 % of the MCSE of the bias in the simulation)
# if job_true takes too long: set back to 1e6 (100 blocks, about 4x faster, MCSE ~0.002), then also --time in job_true.sh

# For small test on cluster: (set by "bash cluster_jobs/submit_all.sh test"),
# without the environment variables N_SIM / N_TRUE nothing changes
if (Sys.getenv("N_SIM")  != "") n_sim  <- as.integer(Sys.getenv("N_SIM"))
if (Sys.getenv("N_TRUE") != "") n_true <- as.numeric(Sys.getenv("N_TRUE"))
if (Sys.getenv("N_SIM") != "" || Sys.getenv("N_TRUE") != "") message("TEST override: n_sim = ", n_sim, ", n_true = ", n_true)
n_chunk_true <- 10000          # patients per arm and block (n_true must be multiple of n_chunk_true)
n_cores_true <- 4              # workers for the true values (each needs memory for one block, ~0.5-1 GB per 10000), for local
# causal parameters: everything that changes Yfull or Death
causal_vars <- c("true_beta_arm", "max_frailty_red", "mortality_level", "eff_frailty_y", "eff_frailty_death",   
                 "eff_frailty_train", "train_frailty_eff")

# scenario blocks (defined in scenarios.R)
#   "main": 6 core scenarios without mortality (2_Informative H0 + low effect, 3 with low effect)
#   "supp": all other scenarios (1_NoFrailty, 4_ObsDropout, 5_NoIndirect, high effect, mortality "high")
scenario_blocks <- c("main")                                        # c("main"), c("supp") or c("main", "supp")
results_subdir  <- paste(sort(unique(scenario_blocks)), collapse = "_")   # cluster: results in <workspace>/results/<results_subdir>


# target values for intercepts
target_falls   <- 0.30                                          # percentage of units with at least one fall per year (control arm)
target_dropout <- 0.20                                          # perc with dropout per year (control arm; training arm: drop_beta_arm and frailty reduction)
target_death   <- c(none = 0, high = 0.02)                      # mortality per year (none = no death)
target_train   <- 0.63                                          # share of persons training in week 1 (training arm)



# --- Variables for analyse -------------------------------

# Bootstrap (own script 06_bootstrap_cluster.R, only the nonparametric MCF)
n_boot <- 200  # should be 200 or more
# test run on the cluster (bash cluster_jobs/submit_all.sh test sets N_BOOT=3), without N_BOOT nothing changes
if (Sys.getenv("N_BOOT") != "") {
  n_boot <- as.integer(Sys.getenv("N_BOOT"))
  message("TEST override: n_boot = ", n_boot)
}
boot_max_run <- 20 # bootstrap only for runs with sim_run <= boot_max_run, Inf = all runs, 0 = no bootstrap
boot_scenarios <- c("2_Informative", "4_ObsDropout") # bootstrap only in these scenarios (scenario_name), NULL = all scenarios
boot_seed_offset <- 2e8 # bootstrap seed = seed of run + offset 