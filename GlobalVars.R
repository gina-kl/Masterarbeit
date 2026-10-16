############################################################################
### GlobalVars #############################################################
############################################################################
# --- Variables for DGP -------------------------------
# Simulations
n_sim <- 10                    # Number of simulations per scenario (Set to 500 or 1000 later)
n_true  <- 20000               # patients per arm for Monte Carlo for true values (cluster: 1e6?)
n_chunk_true <- 10000          # patients per arm and block (n_true must be multiple of n_chunk_true)
n_cores_true <- 4              # workers for the true values (each needs memory for one block, ~0.5-1 GB per 10000)
# causal parameters: everything that changes Yfull or Death 
causal_vars <- c("true_beta_arm", "mortality_level", "eff_frailty_y", "eff_frailty_death",
                 "eff_frailty_train", "train_frailty_eff")


# target values for intercepts
target_falls   <- 0.30                                          # percentage of units with at least one fall per year (control arm)
target_dropout <- 0.20                                          # perc with dropout per year (control arm)
target_death   <- c(none = 0, low = 0.002, high = 0.02)         # mortality per year (none = no death)
target_train   <- 0.63                                          # share of persons training in week 1 (training arm)



# --- Variables for analyse -------------------------------

# Bootstrap
do_bootstrap <- TRUE 
n_boot <- 10
boot_max_run <- 20 # bootstrap only for runs with sim_run <= boot_max_run, Inf = all runs
boot_scenarios <- c("2_Informative", "4_ObsDropout") # bootstrap only in these scenarios (scenario_name), NULL = all scenarios
boot_seed_offset <- 5e6 # bootstrap seed = seed of run + offset