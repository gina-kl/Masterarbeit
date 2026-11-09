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

# random numbers: one seed for the whole study, L'Ecuyer-CMRG streams (Section 7 in functions.R)
#   one stream per data set (causal_id, sim_run) and per block of the true values (causal_id, block),
#   substreams within a data stream: 0 = complete data, 1 = dropout, 2 = proxy, 2 + b = bootstrap sample b
master_seed <- 20261104

# scenario blocks (defined in scenarios.R)
#   "main": 6 core scenarios without mortality: 2_Informative (H0, low effect) and 4_ObsDropout (low effect, positive control), each with drop_beta_arm 0&-1
#   "supp": all other scenarios (1_NoFrailty, 3_NoInfDropout (negative control), 5_NoIndirect, high effect, mortality "high")
#   "mort": the main scenarios with very high mortality (15 %)
#   "sens": sensitivity analyses of 2_Informative with beta_C = -1: eff_frailty_y = 0.8 (H0, low effect), eff_frailty_train = 0 (low effect)

scenario_blocks <- c("main")                 # c("main"), c("supp"), c("mort") or combinations
if (Sys.getenv("SCENARIO_BLOCKS") != "") scenario_blocks <- strsplit(Sys.getenv("SCENARIO_BLOCKS"), ",")[[1]]

# length of the analysis intervals in weeks (to_intervals(), functions.R), must divide 52
# 1 = weekly reports, 13 = every 3 months (as in study)
# 4 = monthly: no quarterly MCF rows (quarters are no whole intervals)
# cluster: set with the environment variable, GlobalVars.R stays at 1:
#   INTERVAL_WEEKS=13 bash cluster_jobs/submit_all.sh analysis
interval_weeks <- 1
if (Sys.getenv("INTERVAL_WEEKS") != "") interval_weeks <- as.integer(Sys.getenv("INTERVAL_WEEKS"))

# folders (cluster_setup.R)
#   data_subdir:    data sets, logs, true values -> only the blocks (01)
#   results_subdir: analysis, bootstrap, aggregated results -> blocks + interval_weeks (03, 06)
#   -> an interval analysis uses the data and true values of the weekly run with the same blocks
data_subdir    <- paste(sort(unique(scenario_blocks)), collapse = "_")
results_subdir <- paste0(data_subdir, if (interval_weeks > 1) paste0("_int", interval_weeks))


# target values for intercepts
target_falls   <- 0.35                                          # percentage of units with at least one fall per year (control arm)
target_dropout <- 0.20                                          # perc with dropout per year (control arm; training arm: drop_beta_arm and frailty reduction)
target_death   <- c(none = 0, high = 0.02, very_high = 0.15)    # mortality per year (none = no death)
target_train   <- 0.63                                          # share of persons training in week 1 (training arm)

# proxy strength
# synthetic proxy P(t) = FrailtyT(t) + noise, measured in both arms
proxy_sd <- c(p025 = 0.25,p4 = 4, p1=1, p05=0.5, p2=2)   # no P models: proxy_sd <- NULL
proxy_scenarios <- c("2_Informative", "2_Informative_Y08")  # NOT <-NULL (then computet for all sceanrios)
proxy_icc <- 0   # share of the proxy error that is person-specific (0 = new error at every measurement, 1 = same error at each measurement, 0.5= half new, half the same)


# --- Variables for analyse -------------------------------

# Bootstrap (own script 06_bootstrap_cluster.R, only the nonparametric MCF)
n_boot <- 200  # should be 200 or more
# test run on the cluster (bash cluster_jobs/submit_all.sh test sets N_BOOT=3), without N_BOOT nothing changes
if (Sys.getenv("N_BOOT") != "") {
  n_boot <- as.integer(Sys.getenv("N_BOOT"))
  message("TEST override: n_boot = ", n_boot)
}
boot_max_run <- 200 # MCSE of the coverage about 1.5 percentage points; bootstrap only for runs with sim_run <= boot_max_run, Inf = all runs, 0 = no bootstrap
boot_scenarios <- c("2_Informative", "4_ObsDropout") # bootstrap only in these scenarios (scenario_name), NULL = all scenarios
# bootstrap sample b uses substream 2 + b of the data stream of the run (functions.R, Section 7)