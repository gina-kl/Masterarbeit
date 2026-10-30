
#-- Multiple Simulation: Analyse over all sceanrios and all runs -----------------------
#   Requires: n_sim in 01.generate weekly > 1
#   all Odat_scen*_run*.rds data in results_dir (lokal ./results/, Cluster: Workspace)
#   local n_sim = 1 (GlobalVars.R) -> empirical_SE, mcse_bias, SE_ratio are NA (sd of one value);
#   for a local test with these columns set n_sim >= 2
#   Bootstrap: own script 06_bootstrap_cluster.R, its results are joined here in "aggregate"

library(ipw)
#single packages instead of tidyverse (tidyverse cannot be installed on the cluster)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(magrittr)
library(survival)
library(splines)
library(mets)
library(foreach)
library(doParallel)

source("cluster_setup.R") # run_on_cluster, n_cores_use, results_dir, part
source("functions.R")     # analysis functions: prep_odat, compute_weights, fit_models, mcf_estimate, mcf_mets, ...
source("GlobalVars.R")    # n_sim
source("scenarios.R")     # base_scenarios of the chosen blocks (array task -> scenario_id)

# array task k = k-th scenario of the chosen blocks (as in 01)
scen_id <- NA
if (run_on_cluster && !is.na(array_id)) {
  if (!(array_id %in% seq_len(nrow(base_scenarios)))) stop("Array-Task ", array_id, " don't exist: only 1 - ", nrow(base_scenarios), " scenarios")
  scen_id <- base_scenarios$scenario_id[array_id]
}

# -- analyse function -------------------------------------------------------------------

analyse_one_run <- function(file_path, sim_log, true_estimands) {

  Odat_all <- readRDS(file_path)

  current_params <- sim_log %>% filter(file_path == !!file_path) # !! change file_name to variable

  # true values per method (one row per causal scenario)
  truth <- true_estimands %>%
    filter(Scenario_Name == current_params$scenario_name[1],
           Mortality     == current_params$mortality_level[1],
           True_Beta_Arm == current_params$true_beta_arm[1])

  if (nrow(truth) != 1) {
    return(tibble(file_path = file_path, error = "no exact match in true_estimands",
                  scenario_id = current_params$scenario_id[1], sim_run = current_params$sim_run[1]))   # scenario_id, sim_run (completeness check)
  }

  # -- Data Preparation -------------------------------------------------

  Odat_all <- prep_odat(Odat_all)                                # log_base_falls, tstop, log_falls_lag (functions.R)

  # synthetic proxy P(t) = FrailtyT(t) + sd * eps (GlobalVars.R: proxy_sd, proxy_seed, proxy_scenarios)
  # one eps per person and week; seed of the complete data + offset -> same noise in all scenarios
  # with the same causal parameters; the same eps for all three strengths -> the levels are paired
  # only without mortality (there GL = LWYY -> only LWYY fits with P), no bootstrap for the P models (06)
  use_proxy <- length(proxy_sd) > 0 &&
    (is.null(proxy_scenarios) || current_params$scenario_name[1] %in% proxy_scenarios) &&
    current_params$mortality_level[1] == "none"
  if (use_proxy) {
    set.seed(current_params$seed[1] + proxy_seed)
    eps <- rnorm(nrow(Odat_all))
    Odat_all <- Odat_all %>%
      mutate(P_p2  = frailty_t + proxy_sd[["p2"]]  * eps,
             P_p1  = frailty_t + proxy_sd[["p1"]]  * eps,
             P_p05 = frailty_t + proxy_sd[["p05"]] * eps)
  }

  Odat     <- Odat_all %>% filter(observed) # observed data in the study

  # realised strength of all proxies: correlation with the true frailty
  #   L1: one row per person, with the baseline frailty
  #   falls history, training share, P: over the observed person-weeks, with the current frailty FrailtyT(t)
  #   training share only in the training arm and from week 2 on (before that it is 0 for everyone)
  per_person <- Odat_all %>% filter(t == 1)
  obs_train  <- Odat %>% filter(arm == 1, t > 1)
  safe_cor   <- function(x, y) if (sd(x) > 0 && sd(y) > 0) cor(x, y) else NA_real_   # NA if one side is constant
  proxy_cor <- tibble(
    cor_L1        = safe_cor(per_person$log_base_falls, per_person$frailty),
    cor_falls_lag = safe_cor(Odat$log_falls_lag, Odat$frailty_t),
    cor_train     = safe_cor(obs_train$train_rate_lag, obs_train$frailty_t),
    cor_P_p2      = if (use_proxy) safe_cor(Odat$P_p2,  Odat$frailty_t) else NA_real_,
    cor_P_p1      = if (use_proxy) safe_cor(Odat$P_p1,  Odat$frailty_t) else NA_real_,
    cor_P_p05     = if (use_proxy) safe_cor(Odat$P_p05, Odat$frailty_t) else NA_real_
  )

  # realised dropout per arm (control arm calibrated to 20 %, training arm depends on beta_C and the frailty reduction)
  per_pt    <- Odat_all %>% distinct(ID, arm, ever_dropped)      # one row per person
  drop_ctrl <- mean(per_pt$ever_dropped[per_pt$arm == 0])
  drop_int  <- mean(per_pt$ever_dropped[per_pt$arm == 1])

  # -- Weights --------------------------------------------------------------------
  # stabilised: numerator = baseline covariates + time (also in the outcome model)
  # ipwtm fitted on weeks without death; weights shifted by one week with lag_ipw()
  # unstabilised weights (no numerator) for the nonparametric mean number of falls
  # *_tr: truncated at the 99% percentile per arm (sensitivity analysis)
  Odat_w <- tryCatch(compute_weights(Odat, use_proxy = use_proxy), error = function(e) conditionMessage(e))

  if (!is.data.frame(Odat_w)) {
    return(tibble(file_path = file_path, error = paste("ipwtm ist fehlgeschlagen:", Odat_w),
                  scenario_id = current_params$scenario_id[1], sim_run = current_params$sim_run[1]))
  }

  # weight diagnostics (before truncation): maximum, 99% percentile and relative effective sample size (should be ~1)
  # of the stabilised weights w_* (models) and the unstabilised weights wm_* (MCF), e.g. column max_wm_sym
  w_cols <- c("w_sym", "w_symtrain", "w_oracle", "wm_sym", "wm_symtrain", "wm_oracle",
              if (use_proxy) c("w_P_p2", "w_P_p1", "w_P_p05", "wm_P_p2", "wm_P_p1", "wm_P_p05"))
  w_diag <- Odat_w %>%
    summarise(across(all_of(w_cols),
                     list(max = max,
                          q99 = ~ quantile(.x, 0.99, names = FALSE),
                          ess = ~ sum(.x)^2 / sum(.x^2) / n()),
                     .names = "{.fn}_{.col}"))

  # -- Models ------------------------------------
  # prep_counting(): prepare data for counting data

  fits <- tryCatch(fit_models(Odat_all, Odat_w, use_proxy), error = function(e) conditionMessage(e))

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
  # true values: LWYY -> limit of the LWYY fit on complete data (MCF_LWYY_*, as True_Beta_LWYY),
  #              GL and nonparametric MCF -> real world with death (MCF_*)
  mcf_est <- mcf_estimate(Odat_all, Odat_w, fits, use_proxy)

  mcf_results <- tibble(Model                = sub("^MCF-", "", names(mcf_est)),
                        MCF_Control_hat      = sapply(mcf_est, `[`, 1),
                        MCF_Intervention_hat = sapply(mcf_est, `[`, 2)) %>%
    mutate(lwyy                  = startsWith(Model, "LWYY"),
           true_mcf_control      = if_else(lwyy, truth$MCF_LWYY_Control,      truth$MCF_Control),
           true_mcf_intervention = if_else(lwyy, truth$MCF_LWYY_Intervention, truth$MCF_Intervention)) %>%
    dplyr::select(-lwyy)

  # nonparametric MCF: own rows, estimate = log ratio of mean falls (real world)
  # robust SE of the log ratio from mets (mcf_mets, functions.R) -> CI, test and coverage as for LWYY / GL
  # mean falls of mets (MCF_*_mets) only as check against weighted_mcf (MCF_*_hat)
  mcf_np_est <- mcf_results %>% filter(startsWith(Model, "MCF"))
  mets <- tryCatch(mcf_mets(Odat_all, Odat_w, use_proxy),            # mets failed -> NA, run is kept
                   error = function(e) tibble(Model = mcf_np_est$Model, MCF_Control_mets = NA_real_,
                                              MCF_Intervention_mets = NA_real_, SE_log_ratio = NA_real_))
  mcf_np_est <- mcf_np_est %>% left_join(mets, by = "Model")
  mcf_np <- wald_results(mcf_np_est$Model,
                         est       = log(mcf_np_est$MCF_Intervention_hat / mcf_np_est$MCF_Control_hat),
                         se        = mcf_np_est$SE_log_ratio,
                         true_beta = truth$True_LogRateRatio) %>%
    left_join(mcf_np_est %>% dplyr::select(-SE_log_ratio), by = "Model")

  # nonparametric MCF per quarter of year (for time varying effect)
  # only under H1: with beta = 0 there is no effect that could vary over time (true LogRR_Q* = 0)
  q_breaks <- c(0, 13, 26, 39, 52)                                 # weeks
  mcf_q <- if (current_params$true_beta_arm[1] != 0) {
    q_list <- list(
      "0: Full data (oracle M0)" = weighted_mcf(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full", breaks = q_breaks),
      "1: Naive (Unweighted)"    = weighted_mcf(Odat_w %>% mutate(w1 = 1), "w1", breaks = q_breaks),
      "3: IPCW sym"              = weighted_mcf(Odat_w, "wm_sym", breaks = q_breaks),
      "4: IPCW sym_train"        = weighted_mcf(Odat_w, "wm_symtrain", breaks = q_breaks),
      "5: IPCW oracle"           = weighted_mcf(Odat_w, "wm_oracle", breaks = q_breaks))
    if (use_proxy) {
      q_list[["8: IPCW sym+P sd=2"]]    <- weighted_mcf(Odat_w, "wm_P_p2",  breaks = q_breaks)
      q_list[["9: IPCW sym+P sd=1"]]    <- weighted_mcf(Odat_w, "wm_P_p1",  breaks = q_breaks)
      q_list[["10: IPCW sym+P sd=0.5"]] <- weighted_mcf(Odat_w, "wm_P_p05", breaks = q_breaks)
    }
    bind_rows(q_list, .id = "version") %>%
      mutate(q = paste0("Q", period)) %>%                              # period 1-4 -> Q1-Q4
      dplyr::select(-period) %>%
      pivot_wider(names_from = arm, values_from = MCF, names_prefix = "arm") %>%
      mutate(Model                = paste0("MCF-", q, " ", version),          # "MCF-Q1 3: IPCW sym" -> Method "MCF-Q1"
             True_Beta            = vapply(q, function(x) truth[[paste0("LogRR_", x)]], numeric(1)),
             Estimate             = log(arm1 / arm0),
             Rate_Ratio           = exp(Estimate),
             Bias                 = Estimate - True_Beta,
             MCF_Control_hat      = arm0,
             MCF_Intervention_hat = arm1) %>%
      dplyr::select(Model, True_Beta, Estimate, Rate_Ratio, Bias, MCF_Control_hat, MCF_Intervention_hat)
  }


  # -- Extract ------------------------------------------
  # LWYY / GL rows get their model-based mean falls in the same row
  results <- bind_rows(model_results %>% left_join(mcf_results, by = "Model"),
                       mcf_np, mcf_q)

  results %>%
    cross_join(w_diag) %>%                                  # weight diagnostics (one row) in every row
    cross_join(proxy_cor) %>%                               # proxy strength (one row) in every row
    mutate(
      Method = sub(" .*", "", Model),                       # LWYY / GL / MCF
      file_path         = file_path,
      scenario_id       = current_params$scenario_id[1],
      scenario_name     = current_params$scenario_name[1],
      sim_run           = current_params$sim_run[1],
      true_beta_arm_dgp = current_params$true_beta_arm[1],
      effect            = current_params$effect[1],          # none / low / high
      drop_beta_arm     = current_params$drop_beta_arm[1],
      mortality_level   = current_params$mortality_level[1],
      drop_ctrl         = drop_ctrl,                         # realised dropout per arm
      drop_int          = drop_int,
      error = NA_character_
    )
}


# Cluster: "analyse" = one scenario per array-task, "aggregate" = sum up all
run_analyse   <- part %in% c("all", "analyse")
run_aggregate <- part %in% c("all", "aggregate")

if (!run_analyse && !run_aggregate) stop("03_multiple_analyse: argument must be 'analyse' or 'aggregate', not: ", part)

# results: one file per scenario (cluster) or one file for all (local)
# every job computes its scenario completely and overwrites the old file
# (a scenario takes about 20-30 min on the cluster -> after a time limit the scenario is simply computed again)
res_dir <- file.path(results_dir, "analysis")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# functions needed on the workers
ana_funs <- c("analyse_one_run", "prep_odat", "extract_model_results", "wald_results", "prep_counting", "lag_ipw",
              "weighted_mcf", "compute_weights", "trunc_w", "model_mcf", "fit_models", "mcf_estimate",
              "mets_mcf", "mcf_mets", "tie_offset")


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
    sim_log <- sim_log %>% filter(scenario_id == scen_id)
    if (nrow(sim_log) == 0) stop("no data for scenario ", scen_id, " in Log")
    res_file <- file.path(res_dir, sprintf("res_scen%02d.rds", scen_id))   # an old file is overwritten
  } else {
    # time check (only local)
    t1 <- Sys.time()
    test_result <- analyse_one_run(sim_log$file_path[1], sim_log, true_estimands)
    print(Sys.time() - t1)
    print(test_result)

    # Optional zum Testen: erstmal nur eine Teilmenge laufen lassen, z.B.:
    #sim_log <- sim_log %>% filter(sim_run <= 2)

    # lokal remove old results
    file.remove(list.files(res_dir, full.names = TRUE))
    res_file <- file.path(res_dir, "res_local.rds")
  }

  cl <- makeCluster(n_cores_use)
  registerDoParallel(cl)

  t_start <- Sys.time()
  cat("Start analysis:", nrow(sim_log), "runs on", n_cores_use, "cores\n")

  res <- foreach(
    i = seq_len(nrow(sim_log)),
    .combine = bind_rows,
    .packages = c("ipw", "survival", "splines", "mets", "dplyr", "tidyr", "readr"),
    .export = c(ana_funs, "proxy_sd", "proxy_seed", "proxy_scenarios")
  ) %dopar% {
    # every unexpected error (e.g. readRDS, MCF, extraction) only marks this run as failed,
    # otherwise foreach would stop and all runs of the scenario would be lost
    tryCatch(
      analyse_one_run(sim_log$file_path[i], sim_log, true_estimands),
      error = function(e) tibble(file_path   = sim_log$file_path[i],
                                 error       = paste("unexpected error:", conditionMessage(e)),
                                 scenario_id = sim_log$scenario_id[i],
                                 sim_run     = sim_log$sim_run[i]))
  }

  stopCluster(cl)
  save_rds_atomic(res, res_file)   # atomic -> no broken result file after a time limit

  cat("Laufzeit gesamt:", round(difftime(Sys.time(), t_start, units = "mins"), 2), "Minuten\n")

}   # --------- end of "analyse" ------------------------------------------------


if (run_aggregate) {   # ---------------- "aggregate" ------------------------

  # only the result files of the scenarios in the chosen blocks (old files of other designs are not read)
  res_files <- if (run_on_cluster) {
    file.path(res_dir, sprintf("res_scen%02d.rds", base_scenarios$scenario_id))
  } else {
    file.path(res_dir, "res_local.rds")
  }
  res_files <- res_files[file.exists(res_files)]
  if (length(res_files) == 0) stop("Keine Analyse-Ergebnisse in ", res_dir)
  all_results <- bind_rows(lapply(res_files, readRDS))

  # check if complete
  # failed runs now also have scenario_id -> they count as analysed and are no extra NA group
  runs_per_scen <- all_results %>% distinct(scenario_id, file_path) %>% count(scenario_id, name = "n_runs")
  cat("analysed runs per scenario (expected per", n_sim, "):\n")
  print(runs_per_scen, n = 100)
  missing_ids <- setdiff(base_scenarios$scenario_id, runs_per_scen$scenario_id)   # scenarios of the chosen blocks
  if (length(missing_ids) > 0) {
    warning("Scenarios without results: ", paste(missing_ids, collapse = ", "), " - start these analyse tasks again")
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

  # check: nonparametric MCF of mets = weighted_mcf (largest difference over all runs and both arms)
  d_mets <- abs(c(all_results$MCF_Control_mets - all_results$MCF_Control_hat,
                  all_results$MCF_Intervention_mets - all_results$MCF_Intervention_hat))
  cat("largest difference mets - weighted_mcf:", if (all(is.na(d_mets))) NA else max(d_mets, na.rm = TRUE), "\n")

  # -- Bootstrap (06_bootstrap_cluster.R): one file per run, joined by scenario, run and model ------------
  # runs / models without bootstrap -> NA
  boot_files <- list.files(file.path(results_dir, "boot"), pattern = "^boot_.*\\.rds$", full.names = TRUE)
  boot_res   <- bind_rows(lapply(boot_files, readRDS))
  if (nrow(boot_res) == 0) {
    boot_res <- tibble(scenario_id = numeric(), sim_run = numeric(), Model = character(), Boot_SE = numeric(),
                       Boot_CI_low = numeric(), Boot_CI_high = numeric(), n_boot_ok = numeric())
  }
  cat(n_distinct(boot_res[c("scenario_id", "sim_run")]), "runs with bootstrap\n")

  all_results <- all_results %>%
    left_join(boot_res, by = c("scenario_id", "sim_run", "Model")) %>%
    mutate(
      Boot_Covered_wald = Estimate - qnorm(0.975) * Boot_SE <= True_Beta &
        True_Beta <= Estimate + qnorm(0.975) * Boot_SE,     # Wald CI with bootstrap SE
      Boot_Covered_perc = Boot_CI_low <= True_Beta & True_Beta <= Boot_CI_high,  # percentile CI
      Boot_Reject       = 2 * pnorm(-abs(Estimate / Boot_SE)) < 0.05,         # Wald test with bootstrap SE
      Boot_Reject_perc  = Boot_CI_high < 0 | Boot_CI_low > 0                  # 0 outside the percentile CI (two-sided)
    )

  write_csv(all_results, file.path(results_dir, "00_all_results_raw.csv"))   # cluster: results_dir

  # -- Results (performance measures) ---------------------------------------------
  # bias against the true value of each method, coverage with MCSE,
  # model SE vs empirical SE, nonparametric mean number of falls

  aggregated_results <- all_results %>%
    group_by(scenario_id, scenario_name, effect, drop_beta_arm, mortality_level,
             true_beta_arm_dgp, Method, Model) %>%
    summarise(
      n_runs        = n(),
      true_beta     = first(True_Beta),                                  # true value of this method
      mean_estimate = mean(Estimate, na.rm = TRUE),                      # mean log rate / mean ratio
      mean_bias     = mean(Bias, na.rm = TRUE),
      mcse_bias     = sd(Estimate, na.rm = TRUE) / sqrt(n_runs),         # Monte Carlo SE of the bias
      empirical_SE  = sd(Estimate, na.rm = TRUE),
      mean_model_SE = mean(SE, na.rm = TRUE),                            # robust sandwich SE (MCF: robust SE from mets)
      SE_ratio      = mean_model_SE / empirical_SE,                      # ~1: SE fits
      RMSE          = sqrt(mean(Bias^2, na.rm = TRUE)),
      coverage      = mean(Covered, na.rm = TRUE),                       # nominal 0.95 (NA for MCF per quarter)
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_runs),
      power = mean(Reject, na.rm = TRUE),                       # beta = -0.4 -> power
      # bootstrap (only runs with bootstrap; NaN / 0 without bootstrap)
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
      # realised dropout per arm (mean over the runs)
      mean_drop_ctrl = mean(drop_ctrl, na.rm = TRUE),
      mean_drop_int  = mean(drop_int, na.rm = TRUE),
      # realised proxy strength: mean over the runs, e.g. mean_cor_L1, mean_cor_P_p05 (NaN where not computed)
      across(starts_with("cor_"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      # weights (before truncation): mean over the runs of max, 99% percentile and relative ESS, e.g. mean_max_wm_sym
      across(matches("^(max|q99|ess)_w"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
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
