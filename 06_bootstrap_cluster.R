
#-- Bootstrap of the nonparametric mean number of falls (MCF) -----------------------------
#   only the nonparametric MCF (log ratio of the mean falls until week 52, all weight versions)
#   runs with sim_run <= boot_max_run in the scenarios boot_scenarios (GlobalVars.R)
#   per run: n_boot samples (persons drawn with replacement per arm, the weights are estimated again
#            in every sample), samples in parallel, one result file per run in <results>/boot
#   cluster: one array task per scenario (as job_data / job_analyse), needs only the data (job_data);
#            runs with a bootstrap file newer than their data set are skipped (restart after a TIMEOUT);
#            after new data (job_data) all runs are computed again
#   the results are joined to the analysis in 03 ("aggregate")

library(ipw)
#single packages instead of tidyverse (tidyverse cannot be installed on the cluster)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(magrittr)
library(splines)
library(foreach)
library(doParallel)

source("cluster_setup.R") # run_on_cluster, n_cores_use, results_dir, part
source("functions.R")     # boot_sample, boot_estimate_mcf, summarise_boot, compute_weights, weighted_mcf, ...
source("GlobalVars.R")    # n_boot, boot_max_run, boot_scenarios, master_seed
source("scenarios.R")     # base_scenarios of the chosen blocks (array task -> scenario_id)

if (!part %in% c("all", "boot")) stop("06_bootstrap: argument must be 'boot', not: ", part)

boot_dir <- file.path(results_dir, "boot")
dir.create(boot_dir, showWarnings = FALSE, recursive = TRUE)

# functions needed on the workers
boot_funs <- c("boot_sample", "boot_estimate_mcf", "compute_weights", "trunc_w", "lag_ipw",
               "weighted_mcf", "mcf_estimate", "prep_odat")

# runs with bootstrap (local 00_simulation_log.csv, cluster 00_simulation_log_scenXX.csv)
log_files <- list.files(data_dir, pattern = "^00_simulation_log.*\\.csv$", full.names = TRUE)
if (length(log_files) == 0) stop("no Simulationslog in ", data_dir, " - first run 01_generate")
sim_log <- bind_rows(lapply(log_files, read_csv, show_col_types = FALSE)) %>%
  mutate(file_path = file.path(data_dir, basename(file_path))) %>%
  distinct(file_path, .keep_all = TRUE) %>%
  filter(sim_run <= boot_max_run,
         is.null(boot_scenarios) | scenario_name %in% boot_scenarios) %>%
  arrange(scenario_id, sim_run)

if (run_on_cluster) {
  # array task k = k-th scenario of the chosen blocks (as in 01 and 03)
  if (!(array_id %in% seq_len(nrow(base_scenarios)))) stop("Array-Task ", array_id, " don't exist: only 1 - ", nrow(base_scenarios), " scenarios")
  scen_id <- base_scenarios$scenario_id[array_id]
  sim_log <- sim_log %>% filter(scenario_id == scen_id)
  # restart after a TIMEOUT: runs whose bootstrap file is newer than their data set are skipped (see the loop)
} else {
  file.remove(list.files(boot_dir, full.names = TRUE))   # local: remove old bootstrap results
}

if (nrow(sim_log) == 0) {
  cat("No bootstrap in this scenario (boot_scenarios / boot_max_run in GlobalVars.R)\n")
} else {

  cl <- makeCluster(n_cores_use)
  registerDoParallel(cl)
  cat("Bootstrap:", nrow(sim_log), "runs x", n_boot, "samples on", n_cores_use, "cores\n")
  t_start <- Sys.time()

  for (i in seq_len(nrow(sim_log))) {
    p        <- sim_log[i, ]
    out_file <- file.path(boot_dir, sprintf("boot_scen%02d_run%04d.rds", p$scenario_id, p$sim_run))
    # done in an earlier attempt (file newer than the data set; older files belong to older data and are computed again)
    if (run_on_cluster && file.exists(out_file) && file.mtime(out_file) > file.mtime(p$file_path)) next
    t_run    <- Sys.time()
    Odat_all <- prep_odat(readRDS(p$file_path), interval_weeks)
    s_run    <- rng_stream(p$stream_id)   # stream of this data set (functions.R, Section 7), once per run

    est <- foreach(b = seq_len(n_boot),
                   .packages = c("ipw", "splines", "dplyr", "tidyr"),
                   .export = c(boot_funs, "interval_weeks", "sub_stream", "use_stream", "s_run")) %dopar% {
      use_stream(sub_stream(s_run, 2 + b))   # bootstrap sample b: substream 2 + b (0 data, 1 dropout, 2 proxy)
      tryCatch(boot_estimate_mcf(boot_sample(Odat_all), interval_weeks), error = function(e) NULL)   # failed sample -> NULL
    }
    est <- do.call(rbind, est)                                   # one row per sample, failed samples (NULL) are dropped
    if (is.null(est)) {
      warning("all bootstrap samples failed: ", basename(p$file_path))
      next
    }

    save_rds_atomic(summarise_boot(est) %>% mutate(scenario_id = p$scenario_id, sim_run = p$sim_run), out_file)
    cat(format(Sys.time(), "%H:%M:%S"), basename(out_file), ":", nrow(est), "samples in",
        round(difftime(Sys.time(), t_run, units = "mins"), 2), "Minuten\n")
  }

  stopCluster(cl)
  cat("Laufzeit gesamt:", round(difftime(Sys.time(), t_start, units = "mins"), 2), "Minuten\n")
}
