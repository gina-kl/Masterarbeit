################################################################################
# cluster_setup.R – lokal / bwUniCluster
# load with source("cluster_setup.R")
################################################################################

# ---- cluster or local -----------------------------------------------------------
# automatic detection: inside a SLURM job SLURM_JOB_ID is set.
# override possible: set run_on_cluster before source("cluster_setup.R") (e.g. time_check.R)

if (!exists("run_on_cluster")) run_on_cluster <- Sys.getenv("SLURM_JOB_ID") != ""  

if (run_on_cluster) {
  #   Rscript 01_generate_weekly_cluster.R true        -> only true values (single job, no array)
  #   Rscript 01_generate_weekly_cluster.R data        -> only simulated data (one scenario per array task)
  #   Rscript 03_multiple_analyse_cluster.R analyse    -> analysis (one scenario per array task)
  #   Rscript 03_multiple_analyse_cluster.R aggregate  -> sum up all scenarios (single job, no array)
  #   Rscript 06_bootstrap_cluster.R boot              -> bootstrap of the MCF (one scenario per array task)
  valid_parts <- c("true", "data", "analyse", "aggregate", "boot")   # analyse / aggregate for 03, boot for 06
  args <- commandArgs(trailingOnly = TRUE)
  part <- if (length(args) >= 1) args[1] else stop("Missing argument: one of ", paste(valid_parts, collapse = ", "))  
  if (!part %in% valid_parts) stop("Argument should be one of ", paste(valid_parts, collapse = ", "), ", not: ", part)   

  # warning if --cpus-per-task is missing (otherwise silently only 1 core)
  if (Sys.getenv("SLURM_CPUS_PER_TASK") == "") warning("SLURM_CPUS_PER_TASK not set - only 1 core is used (set #SBATCH --cpus-per-task)")
  n_cores_use <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

  # array_task_number = position of the scenario in the chosen blocks (only for part = "data", "analyse" and "boot")
  array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
  if (part %in% c("data", "analyse", "boot") && is.na(array_id)) stop("part = '", part, "' must run as job array (#SBATCH --array=...)")
  # true / aggregate must not run as array: all tasks would compute the same and write the same files
  if (part %in% c("true", "aggregate") && !is.na(array_id)) stop("part = '", part, "' must not run as job array (remove #SBATCH --array)")

  # results in workspace
  sim_dir <- Sys.getenv("SIM_DIR")
  if (sim_dir == "") stop("SIM_DIR missing, set in Job-script")
  if (!startsWith(sim_dir, "/")) stop("SIM_DIR must be an absolute path, not: ", sim_dir)   # NEU
  source("GlobalVars.R")                                        # NEU: results_subdir (scenario blocks)
  results_dir <- file.path(sim_dir, "results", results_subdir)  # NEU: own folder per block selection

} else {
  part        <- "all"                               # local: true_values and all scenarios
  n_cores_use <- max(1, parallel::detectCores() - 1)
  array_id    <- NA
  results_dir <- "./results"
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

message("Modus: ", if (run_on_cluster) "CLUSTER" else "LOKAL",
        " | part: ", part,
        " | kernels: ", n_cores_use,
        if (!is.na(array_id)) paste0(" | Array-Task: ", array_id) else "",
        " | results: ", results_dir)
