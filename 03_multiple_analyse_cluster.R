
#-- Multiple Simulation: Analyse over all sceanrios and all runs -----------------------
#   Requires: n_sim in 01.generate weekly > 1 
#   all Odat_scen*_run*.rds data in results_dir (lokal ./results/, Cluster: Workspace)
#   local n_sim = 1 (GlobalVars.R) -> empirical_SE, mcse_bias, SE_ratio are NA (sd of one value);
#   for a local test with these columns set n_sim >= 2

library(ipw)
library(tidyverse)
library(magrittr)
library(survival)
library(splines)
library(mets)
library(foreach)
library(doParallel)

source("cluster_setup.R") # run_on_cluster, n_cores_use, results_dir, part
source("functions.R")   # expand_counting, prep_gl_data, extract_model_results,
# lag_ipw, weighted_mcf, model_mcf
source("GlobalVars.R")  # do_bootstrap, n_boot, boot_max_run, boot_scenarios, boot_seed_offset

# -- all of analyse one function -------------------------------------------------------------------

analyse_one_run <- function(file_path, sim_log, true_estimands, do_bootstrap = FALSE, n_boot = 200, boot_max_run = Inf,
                            boot_scenarios = NULL, boot_seed_offset = 5e6) {
  
  Odat_all <- readRDS(file_path)
  
  current_params <- sim_log %>% filter(file_path == !!file_path) # !! change file_name to variable
  
  # rue values per method (one row per causal scenario)
  truth <- true_estimands %>%
    filter(Scenario_Name == current_params$scenario_name[1],
           Mortality     == current_params$mortality_level[1],
           True_Beta_Arm == current_params$true_beta_arm[1])
  
  if (nrow(truth) != 1) {
    return(tibble(file_path = file_path, error = "no exact match in true_estimands",
                  scenario_id = current_params$scenario_id[1], sim_run = current_params$sim_run[1]))   # NEU: scenario_id, sim_run (completeness check)
  }
  
  # -- Data Preparation -------------------------------------------------
  
  Odat_all <- Odat_all %>%
    mutate(log_base_falls = log1p(base_risk),                     
           tstop          = t,
           fall_rate_lag  = cum_falls_lag / pmax(t - 1, 1)) %>%  # observed fall rate up to week t-1
    arrange(ID, t)
  
  Odat <- Odat_all %>% filter(observed)                          # observed data in the study
  
  # -- Weights --------------------------------------------------------------------
  # stabilised: numerator = baseline covariates + time (also in the outcome model)
  # ipwtm fitted on weeks without death; weights shifted by one week with lag_ipw()
  # unstabilised weights (no numerator) for the nonparametric mean number of falls
  Odat_w <- tryCatch(compute_weights(Odat), error = function(e) conditionMessage(e))
  
  if (!is.data.frame(Odat_w)) {
    return(tibble(file_path = file_path, error = paste("ipwtm ist fehlgeschlagen:", Odat_w),
                  scenario_id = current_params$scenario_id[1], sim_run = current_params$sim_run[1]))   # NEU: scenario_id, sim_run
  }
  
  # weight diagnostics  
  max_weight                  <- max(Odat_w$w_sym, na.rm = TRUE)
  max_weight_symtrain         <- max(Odat_w$w_symtrain, na.rm = TRUE)
  ess_sym <- sum(Odat_w$w_sym)^2 / sum(Odat_w$w_sym^2) / nrow(Odat_w)   # relative effective sample size
  
  # -- Models ------------------------------------
  # prep_counting(): prepare data for counting data
  
  fits <- tryCatch(fit_models(Odat_all, Odat_w) , error = function(e) conditionMessage(e))
  
  if (!is.list(fits)) {
    return(tibble(file_path = file_path, error = paste("Models failed:", fits),                       
                  scenario_id = current_params$scenario_id[1], sim_run = current_params$sim_run[1]))   
  }
  
  # true value that belongs to each model
  model_results <- bind_rows(lapply(names(fits), function(m) {
    true_beta <- if (startsWith(m, "LWYY")) truth$True_Beta_LWYY else truth$True_Beta_GL
    extract_model_results(fits[[m]], m, true_beta)
  }))
  
  # -- Mean number of falls until week 52 (nonparametric and model estimate) --------------------------------
  
  # mean falls per arm (model-based = same model as the LWYY / GL row,
  # log(mu_1 / mu_0) is exactly the arm coefficient -> only the levels per arm are new)
  # true values: LWYY -> world without death, GL and nonparametric MCF -> real world
  mcf_est <- mcf_estimate(Odat_all, Odat_w, fits)
  
  mcf_results <- tibble(Model                = sub("^MCF-", "", names(mcf_est)),
                        MCF_Control_hat      = sapply(mcf_est, `[`, 1),
                        MCF_Intervention_hat = sapply(mcf_est, `[`, 2)) %>%
    mutate(no_death              = startsWith(Model, "LWYY"),
           true_mcf_control      = if_else(no_death, truth$MCF_Control_NoDeath,      truth$MCF_Control),
           true_mcf_intervention = if_else(no_death, truth$MCF_Intervention_NoDeath, truth$MCF_Intervention)) %>%
    dplyr::select(-no_death)
  
  # nonparametric MCF: own rows, estimate = log ratio of mean falls (real world)
  mcf_np <- mcf_results %>%
    filter(startsWith(Model, "MCF")) %>%
    mutate(True_Beta  = truth$True_LogRateRatio,
           Estimate   = log(MCF_Intervention_hat / MCF_Control_hat),
           Rate_Ratio = exp(Estimate),
           Bias       = Estimate - True_Beta)
  
  # nonparametric MCF per quarter of year (for time varying effect)
  # only under H1: with beta = 0 there is no effect that could vary over time (true LogRR_Q* = 0)
  mcf_q <- if (current_params$true_beta_arm[1] != 0) bind_rows(
    "0: Full data (oracle M0)" = weighted_mcf_q(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full"),
    "1: Naive (Unweighted)"    = weighted_mcf_q(Odat_w %>% mutate(w1 = 1), "w1"),
    "3: IPCW sym"              = weighted_mcf_q(Odat_w, "wm_sym"),
    "4: IPCW sym_train"        = weighted_mcf_q(Odat_w, "wm_symtrain"),
    "5: IPCW oracle"           = weighted_mcf_q(Odat_w, "wm_oracle"),
    .id = "version") %>%
    pivot_wider(names_from = arm, values_from = MCF, names_prefix = "arm") %>%
    mutate(Model                = paste0("MCF-", q, " ", version),          # "MCF-Q1 3: IPCW sym" -> Method "MCF-Q1"
           True_Beta            = vapply(as.character(q), function(x) truth[[paste0("LogRR_", x)]], numeric(1)),
           Estimate             = log(arm1 / arm0),
           Rate_Ratio           = exp(Estimate),
           Bias                 = Estimate - True_Beta,
           MCF_Control_hat      = arm0,
           MCF_Intervention_hat = arm1) %>%
    dplyr::select(Model, True_Beta, Estimate, Rate_Ratio, Bias, MCF_Control_hat, MCF_Intervention_hat)
  
  
  # -- Extract ------------------------------------------
  # LWYY / GL rows get their model-based mean falls in the same row
  results <- bind_rows(model_results %>% left_join(mcf_results, by = "Model"),
                       mcf_np, mcf_q)
  
  
  # --- Bootstrap (if set in GlobalVars.R)----------
  # only for runs with sim_run <= boot_max_run in the scenarios boot_scenarios (NULL = all)
  # columns always exist (NA without bootstrap) -> aggregation works in both cases
  
  boot_res <- NULL
  boot_this_run <- do_bootstrap && current_params$sim_run[1] <= boot_max_run &&
    (is.null(boot_scenarios) || current_params$scenario_name[1] %in% boot_scenarios)
  if (boot_this_run) {
    set.seed(current_params$seed[1] + boot_seed_offset)
    boot_res <- tryCatch(run_bootstrap(Odat_all, n_boot), error = function(e) NULL)   # failed -> NA, run is kept
  }
  if (is.null(boot_res)) {
    boot_res <- tibble(Model = results$Model, Boot_SE = NA_real_, Boot_CI_low = NA_real_,
                       Boot_CI_high = NA_real_, n_boot_ok = NA_real_)
  }
  
  
  
  results %>%
    left_join(boot_res, by = "Model") %>%   # join with bootstrap results
    mutate(
      Boot_Covered_wald = Estimate - qnorm(0.975) * Boot_SE <= True_Beta &
        True_Beta <= Estimate + qnorm(0.975) * Boot_SE,     # Wald CI with bootstrap SE
      Boot_Covered_perc = Boot_CI_low <= True_Beta & True_Beta <= Boot_CI_high,  # percentile CI
      Boot_Reject       = 2 * pnorm(-abs(Estimate / Boot_SE)) < 0.05,         # Wald test with bootstrap SE
      Boot_Reject_perc  = Boot_CI_high < 0 | Boot_CI_low > 0,                 # 0 outside the percentile CI (two-sided)
      
      Method = sub(" .*", "", Model),                       # LWYY / GL / MCF 
      file_path         = file_path,
      scenario_id       = current_params$scenario_id[1],
      scenario_name     = current_params$scenario_name[1],
      sim_run           = current_params$sim_run[1],
      true_beta_arm_dgp = current_params$true_beta_arm[1],
      drop_beta_arm     = current_params$drop_beta_arm[1],
      mortality_level   = current_params$mortality_level[1],
      max_weight                  = max_weight,
      max_weight_symtrain         = max_weight_symtrain,
      ess_sym                     = ess_sym,
      error = NA_character_
    )
}


# Cluster: "analyse" = one scenario per array-task, "aggregate" = sum up all
run_analyse   <- part %in% c("all", "analyse")
run_aggregate <- part %in% c("all", "aggregate")

if (!run_analyse && !run_aggregate) stop("03_multiple_analyse: argument must be 'analyse' or 'aggregate', not: ", part)

# cluster: results are saved in chunks
# if one job is cancelled (time), only remaining chunks are computed
# local: one chunk with all runs
chunk_size <- 96   # multiple of 16,24,..( no half-empty last round per chunk)
chunk_dir  <- file.path(results_dir, "analysis_chunks")
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
chunk_file <- function(ch) file.path(chunk_dir, paste0("res_", ch, ".rds"))

# functions needed on the workers (one list for .export and for the settings check)
ana_funs <- c("analyse_one_run", "extract_model_results", "prep_counting", "lag_ipw", "weighted_mcf", "compute_weights",
              "model_mcf", "fit_models", "mcf_estimate", "boot_sample", "boot_estimate", "run_bootstrap",
              "weighted_mcf_q", "tie_offset")


# --- parallalise ----------------------------------------------------------------

if (run_analyse) {   # cluster ------------- "analyse" ----------------------------
  
  # local 00_simulation_log.csv, cluster 00_simulation_log_scenXX.csv
  log_files <- list.files(results_dir, pattern = "^00_simulation_log.*\\.csv$", full.names = TRUE)
  if (length(log_files) == 0) stop("no Simulationslog in ", results_dir, " - first run 01_generate")
  sim_log <- bind_rows(lapply(log_files, read_csv, show_col_types = FALSE)) %>%
    mutate(file_path = file.path(results_dir, basename(file_path))) %>%
    distinct(file_path, .keep_all = TRUE) %>%
    arrange(scenario_id, sim_run)
  true_estimands <- read_csv(file.path(results_dir, "00_true_estimands.csv"), show_col_types = FALSE)

  if (run_on_cluster) {
    # only scenario of this array task
    sim_log <- sim_log %>% filter(scenario_id == array_id)
    if (nrow(sim_log) == 0) stop("no data for scenario ", array_id, " in Log")   # NEU: typo

    # existing chunks are only skipped if they were computed with the same code, settings,
    # true values and data (one settings file per scenario -> array tasks never write the same file)
    data_sig_file <- file.path(results_dir, sprintf("00_data_settings_scen%02d.rds", array_id))
    sig_ana <- list(
      code     = lapply(setNames(ana_funs, ana_funs), function(f) deparse(get(f))),   # deparse ignores comments
      settings = list(do_bootstrap = do_bootstrap, n_boot = n_boot, boot_max_run = boot_max_run,
                      boot_scenarios = boot_scenarios, boot_seed_offset = boot_seed_offset,
                      chunk_size = chunk_size),
      truth    = true_estimands %>%
        semi_join(sim_log %>% distinct(scenario_name, mortality_level, true_beta_arm),
                  by = c(Scenario_Name = "scenario_name", Mortality = "mortality_level", True_Beta_Arm = "true_beta_arm")) %>%
        as.list(),
      data     = if (file.exists(data_sig_file)) readRDS(data_sig_file) else NULL
    )
    check_settings(sig_ana, file.path(chunk_dir, sprintf("settings_scen%02d.rds", array_id)),
                   what = "Analysis chunks")
  } else {
    # time check (only local)
    t1 <- Sys.time()
    test_result <- analyse_one_run(sim_log$file_path[1], sim_log, true_estimands)
    print(Sys.time() - t1)
    print(test_result)
    
    # Optional zum Testen: erstmal nur eine Teilmenge laufen lassen, z.B.:
    #sim_log <- sim_log %>% filter(sim_run <= 2)
    
    # lokal remove old chunks
    file.remove(list.files(chunk_dir, full.names = TRUE))
  }
  
  # cluster
  sim_log <- sim_log %>%
    mutate(chunk = if (run_on_cluster) {
      sprintf("scen%02d_chunk%03d", scenario_id, (sim_run - 1) %/% chunk_size + 1)   # NEU: chunk_size (chunk_size_cluster did not exist)
    } else {
      "local_all"                                   # local: all in ohne chunk
    })
  chunks <- unique(sim_log$chunk)
  todo   <- chunks[!file.exists(chunk_file(chunks))]
  cat(length(chunks) - length(todo), "von", length(chunks), "package finished,",
      length(todo), "to compute (", nrow(sim_log), "runs)\n")
  
  if (length(todo) > 0) {
    
    cl <- makeCluster(n_cores_use)   
    registerDoParallel(cl)
    
    t_start <- Sys.time()
    
    for (ch in todo) {               
      idx <- which(sim_log$chunk == ch)
      t_ch <- Sys.time()
      
      res_chunk <- foreach(
        i = idx,
        .combine = bind_rows,
        .packages = c("ipw", "survival", "splines", "mets", "dplyr", "tidyr", "readr"),
        .export = c(ana_funs,                                                             # NEU: list from above
                    "do_bootstrap", "n_boot", "boot_max_run", "boot_scenarios", "boot_seed_offset")
      ) %dopar% {
        # every unexpected error (e.g. readRDS, MCF, extraction) only marks this run as failed,
        # before, foreach stopped and the whole chunk was lost (and failed again after every restart)
        tryCatch(
          # named arguments: no mix-up of the order
          analyse_one_run(sim_log$file_path[i], sim_log, true_estimands,
                          do_bootstrap = do_bootstrap, n_boot = n_boot, boot_max_run = boot_max_run,
                          boot_scenarios = boot_scenarios, boot_seed_offset = boot_seed_offset),
          error = function(e) tibble(file_path   = sim_log$file_path[i],
                                     error       = paste("unexpected error:", conditionMessage(e)),
                                     scenario_id = sim_log$scenario_id[i],
                                     sim_run     = sim_log$sim_run[i]))
      }

      save_rds_atomic(res_chunk, chunk_file(ch))   # NEU: atomic -> no broken chunk file after a time limit
      cat(format(Sys.time(), "%H:%M:%S"), "Paket", ch, "fertig:", length(idx), "Läufe in",
          round(difftime(Sys.time(), t_ch, units = "mins"), 2), "Minuten\n")
    }
    
    stopCluster(cl)
    
    t_end <- Sys.time()
    cat("Laufzeit gesamt:", round(difftime(t_end, t_start, units = "mins"), 2), "Minuten\n")
    
  }   # Ende if (length(todo) > 0)
  
}   # --------- end of "analyse" ------------------------------------------------


if (run_aggregate) {   # ---------------- "aggregate" ------------------------
  
   chunk_files <- list.files(chunk_dir, pattern = "^res_.*\\.rds$", full.names = TRUE)
  if (length(chunk_files) == 0) stop("Keine Analyse-Ergebnisse in ", chunk_dir)
  all_results <- bind_rows(lapply(chunk_files, readRDS))
  
  # check if complete
  # failed runs now also have scenario_id -> they count as analysed and are no extra NA group
  runs_per_scen <- all_results %>% distinct(scenario_id, file_path) %>% count(scenario_id, name = "n_runs")
  cat("analysed runs per scenario (expected per", n_sim, "):\n")
  print(runs_per_scen, n = 100)
  # expected scenarios from 00_base_scenarios.csv (written by 01, part "true") instead of the fixed 28
  base_file <- file.path(results_dir, "00_base_scenarios.csv")
  if (file.exists(base_file)) {
    expected_ids <- read_csv(base_file, show_col_types = FALSE)$scenario_id
    missing_ids  <- setdiff(expected_ids, runs_per_scen$scenario_id)
    if (length(missing_ids) > 0) {
      warning("Scenarios without results: ", paste(missing_ids, collapse = ", "), " - start these analyse tasks again")
    }
  } else {
    warning(base_file, " missing - number of scenarios not checked")
  }
  if (any(runs_per_scen$n_runs < n_sim)) {
    warning("Not all runs analysed in scenario(s) ",
            paste(runs_per_scen$scenario_id[runs_per_scen$n_runs < n_sim], collapse = ", "),
            " - start these analyse tasks again")
  }
  
  # check failed runs
  failed_runs <- all_results %>% filter(!is.na(error))
  if (nrow(failed_runs) > 0) {
    cat(nrow(failed_runs), "Läufe sind fehlgeschlagen:\n")
    print(failed_runs %>% dplyr::select(file_path, error) %>% distinct())
  }
  
  all_results <- all_results %>% filter(is.na(error))
  
  write_csv(all_results, file.path(results_dir, "00_all_results_raw.csv"))   # cluster: results_dir
  
  # -- Results (performance measures) ---------------------------------------------
  # bias against the true value of each method, coverage with MCSE,
  # model SE vs empirical SE, nonparametric mean number of falls
  
  aggregated_results <- all_results %>%
    group_by(scenario_id, scenario_name, drop_beta_arm, mortality_level,
             true_beta_arm_dgp, Method, Model) %>%
    summarise(
      n_runs        = n(),
      true_beta     = first(True_Beta),                                  # true value of this method
      mean_estimate = mean(Estimate, na.rm = TRUE),                      # mean log rate / mean ratio
      mean_bias     = mean(Bias, na.rm = TRUE),
      mcse_bias     = sd(Estimate, na.rm = TRUE) / sqrt(n_runs),         # Monte Carlo SE of the bias
      empirical_SE  = sd(Estimate, na.rm = TRUE),
      mean_model_SE = mean(SE, na.rm = TRUE),                            # robust sandwich SE (NA for MCF)
      SE_ratio      = mean_model_SE / empirical_SE,                      # ~1: SE fits
      RMSE          = sqrt(mean(Bias^2, na.rm = TRUE)),
      coverage      = mean(Covered, na.rm = TRUE),                       # nominal 0.95 (NA for MCF)
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_runs),
      power = mean(Reject, na.rm = TRUE),                       # beta = -0.4 -> power
      # bootstrap (only runs with bootstrap; NaN / 0 if do_bootstrap = FALSE)
      n_boot_runs        = sum(!is.na(Boot_SE)),
      mean_n_boot_ok     = mean(n_boot_ok, na.rm = TRUE),                # successful bootstrap samples per run
      mean_boot_SE       = mean(Boot_SE, na.rm = TRUE),
      boot_SE_ratio      = mean_boot_SE / empirical_SE,                  # ~1: bootstrap SE fits
      boot_coverage_wald = mean(Boot_Covered_wald, na.rm = TRUE),        # Wald CI with bootstrap SE
      boot_coverage_perc = mean(Boot_Covered_perc, na.rm = TRUE),        # percentile CI
      boot_power         = mean(Boot_Reject, na.rm = TRUE),
      boot_power_perc    = mean(Boot_Reject_perc, na.rm = TRUE),
      # mean number of falls until week 52 (MCF: nonparametric, LWYY / GL: model-based)
      true_mcf_control      = first(true_mcf_control),
      true_mcf_intervention = first(true_mcf_intervention),
      mean_mcf_control      = mean(MCF_Control_hat, na.rm = TRUE),
      mean_mcf_intervention = mean(MCF_Intervention_hat, na.rm = TRUE),
      bias_mcf_control      = mean_mcf_control - true_mcf_control,
      bias_mcf_intervention = mean_mcf_intervention - true_mcf_intervention,
      # weights
      mean_max_weight                  = mean(max_weight, na.rm = TRUE),
      mean_max_weight_symtrain         = mean(max_weight_symtrain, na.rm = TRUE),   
      mean_ess_sym                     = mean(ess_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario_id, Model)
  
  print(aggregated_results, n = 100)
  
  write_csv(aggregated_results, file.path(results_dir, "00_aggregated_results.csv"))   # cluster:  results_dir

  # cluster: one combined log 00_simulation_log.csv (as written locally by 01) with paths "./results/<file>"
  # -> 02_analyse.R / 04_overview.R work locally after copying the cluster results to ./results
  if (run_on_cluster) {
    scen_logs <- list.files(results_dir, pattern = "^00_simulation_log_scen.*\\.csv$", full.names = TRUE)
    bind_rows(lapply(scen_logs, read_csv, show_col_types = FALSE)) %>%
      mutate(file_path = file.path(".", "results", basename(file_path))) %>%
      arrange(scenario_id, sim_run) %>%
      write_csv(file.path(results_dir, "00_simulation_log.csv"))
  }
  
  cat("\nFertig. Ergebnisse gespeichert in:\n")
  cat(" -", file.path(results_dir, "00_all_results_raw.csv"), "(jede einzelne Schätzung)\n")
  cat(" -", file.path(results_dir, "00_aggregated_results.csv"), "(Bias/RMSE/Coverage pro Szenario)\n")
  
}   # -----------------end of  "aggregate" ----------------------------------------------


# # Test for one bootstrap
# 
# Odat_test <- readRDS(sim_log$file_path[1]) %>%
#   mutate(log_base_falls = log1p(base_risk), tstop = t,
#          fall_rate_lag = cum_falls_lag / pmax(t - 1, 1)) %>%
#   arrange(ID, t)
# 
# set.seed(1)
# bs <- boot_sample(Odat_test)
# 
# # Kontrollen der Stichprobe
# n_distinct(bs$ID) == n_distinct(Odat_test$ID)        # gleiche Personenzahl -> TRUE
# table(distinct(bs, ID, arm)$arm)                     # Armgrößen wie im Original
# table(distinct(Odat_test, ID, arm)$arm)
# system.time(est <- boot_estimate(bs))               # Laufzeit EINER Stichprobe
# est                                                  # 25 benannte Werte, möglichst ohne NA
# 
# 
