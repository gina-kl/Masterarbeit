###############################################################################
### Mehrfach-Simulation: Analyse über alle Szenarien und Läufe             ###
### Voraussetzung: n_sim in der Datengenerierung wurde erhöht (z.B. 1000)  ###
### und alle Odat_scen*_run*.rds Dateien liegen in ./results/             ###
###############################################################################

library(ipw)
library(tidyverse)
library(magrittr)
library(survival)
library(splines)
library(mets)
library(foreach)
library(doParallel)

source("functions.R")   # expand_counting, prep_gl_data, extract_model_results,
# lag_ipw, weighted_mcf, model_mcf

# -----------------------------------------------------------------------
# all of analyse one function
# -----------------------------------------------------------------------

analyse_one_run <- function(file_path, sim_log, true_estimands) {
  
  Odat_all <- readRDS(file_path)
  
  current_params <- sim_log %>% filter(file_path == !!file_path)
  
  # rue values per method (one row per causal scenario)
  truth <- true_estimands %>%
    filter(Scenario_Name == current_params$scenario_name[1],
           Mortality     == current_params$mortality_level[1],
           True_Beta_Arm == current_params$true_beta_arm[1])
  
  if (nrow(truth) != 1) {
    return(tibble(file_path = file_path, error = "Kein eindeutiger Match in true_estimands"))
  }
  
  # --- Data Preparation ---
  
  Odat_all <- Odat_all %>%
     mutate(log_base_falls = log1p(base_risk),                     
           tstop          = t,
           fall_rate_lag  = cum_falls_lag / pmax(t - 1, 1)) %>%  # observed fall rate up to week t-1
    arrange(ID, t)
  
  Odat <- Odat_all %>% filter(observed)                          # observed data in the study
  
  # --- Weights ---
  # stabilised: numerator = baseline covariates + time (also in the outcome model)
  # ipwtm fitted on weeks without death; weights shifted by one week with lag_ipw()
  # unstabilised weights (no numerator) for the nonparametric mean number of falls
  Odat_w <- tryCatch(compute_weights(Odat), error = function(e) conditionMessage(e))
  
  if (!is.data.frame(Odat_w)) {
    return(tibble(file_path = file_path, error = paste("ipwtm ist fehlgeschlagen:", Odat_w)))
  }
  
  # weight diagnostics  
  max_weight                  <- max(Odat_w$w_sym, na.rm = TRUE)
  max_weight_symtrain         <- max(Odat_w$w_symtrain, na.rm = TRUE)
  ess_sym <- sum(Odat_w$w_sym)^2 / sum(Odat_w$w_sym^2) / nrow(Odat_w)   # relative effective sample size
  
  # --- Models ---
  # prep_counting(): prepare data for counting data
  
  fits <- tryCatch({
    cp      <- prep_counting(Odat_w, "observed")
    cp_full <- prep_counting(Odat_all, "full")
    
    f_lwyy <- Surv(start, stop, status == 1) ~ arm + age + sex + log_base_falls + cluster(ID)
    f_gl   <- Event(start, stop, status) ~ arm + age + sex + log_base_falls + cluster(ID)
    
    # GL: cens.model = ~ strata(arm) -> internal censoring Kaplan-Meier per arm
    list(
      "LWYY 0: Full data (oracle M0)" = coxph(f_lwyy, data = cp_full, ties = "breslow"),
      "LWYY 1: Naive (Unweighted)"    = coxph(f_lwyy, data = cp, ties = "breslow"),
      "LWYY 3: IPCW sym"              = coxph(f_lwyy, data = cp, ties = "breslow", weights = w_sym),
      "LWYY 4: IPCW sym_train"        = coxph(f_lwyy, data = cp, ties = "breslow", weights = w_symtrain),
      "LWYY 5: IPCW oracle"           = coxph(f_lwyy, data = cp, ties = "breslow", weights = w_oracle),
      "GL 0: Full data (oracle M0)"   = recreg(f_gl, data = cp_full, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm)),
      "GL 1: Naive (Unweighted)"      = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm)),
      "GL 3: IPCW sym"                = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_sym),
      "GL 4: IPCW sym_train"          = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_symtrain),
      "GL 5: IPCW oracle"             = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_oracle)
    )
  }, error = function(e) conditionMessage(e))
  
  if (!is.list(fits)) {
    return(tibble(file_path = file_path, error = paste("Modelle sind fehlgeschlagen:", fits)))
  }
  
  # true value that belongs to each model
  model_results <- bind_rows(lapply(names(fits), function(m) {
    true_beta <- if (startsWith(m, "LWYY")) truth$True_Beta_LWYY else truth$True_Beta_GL
    extract_model_results(fits[[m]], m, true_beta)
  }))
  
  # --- Mean number of falls until week 52 (nonparametric Ghosh-Lin, weighted) ---
  
  mcf_list <- list(
    "MCF 0: Full data (oracle M0)" = weighted_mcf(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full"),
    "MCF 1: Naive (Unweighted)"    = weighted_mcf(Odat_w %>% mutate(w1 = 1), "w1"),
    "MCF 3: IPCW sym"              = weighted_mcf(Odat_w, "wm_sym"),
    "MCF 4: IPCW sym_train"        = weighted_mcf(Odat_w, "wm_symtrain"),
    "MCF 5: IPCW oracle"           = weighted_mcf(Odat_w, "wm_oracle")
  )
  # estimated mean number of falls per arm (control, intervention):
  # nonparametric (weighted_mcf) and model-based from all LWYY and GL fits (model_mcf, functions.R)
  base_cov <- Odat_all %>% filter(t == 1) %>% dplyr::select(ID, age, sex, log_base_falls)   # one row per person
  mcf_est <- c(
    lapply(mcf_list, function(mm) c(mm$MCF[mm$arm == 0], mm$MCF[mm$arm == 1])),
    setNames(lapply(fits, function(f) tryCatch(model_mcf(f, base_cov), error = function(e) c(NA, NA))),
             paste0("MCF-", names(fits)))                                  # "MCF-LWYY ..." / "MCF-GL ..."
  )

  # true values per row: MCF-LWYY -> world without death, all others -> real world
  mcf_results <- bind_rows(lapply(names(mcf_est), function(m) {
    mu       <- mcf_est[[m]]
    no_death <- startsWith(m, "MCF-LWYY")
    tibble(Model                 = m,
           True_Beta             = if (no_death) truth$True_LogRateRatio_NoDeath else truth$True_LogRateRatio,
           Estimate              = log(mu[2] / mu[1]),                    # log ratio of mean falls
           Bias                  = Estimate - True_Beta,
           MCF_Control_hat       = mu[1],
           MCF_Intervention_hat  = mu[2],
           true_mcf_control      = if (no_death) truth$MCF_Control_NoDeath      else truth$MCF_Control,
           true_mcf_intervention = if (no_death) truth$MCF_Intervention_NoDeath else truth$MCF_Intervention)
  }))

  # --- Extract ---
  bind_rows(model_results, mcf_results) %>%
    mutate(
      Method = sub(" .*", "", Model),                       # LWYY / GL / MCF / MCF-LWYY / MCF-GL  
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
# -----------------------------------------------------------------------
# parallalise
# -----------------------------------------------------------------------

sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

# time check
t1 <- Sys.time()
test_result <- analyse_one_run(sim_log$file_path[1], sim_log, true_estimands)
print(Sys.time() - t1)
print(test_result)

# Optional zum Testen: erstmal nur eine Teilmenge laufen lassen, z.B.:
# sim_log <- sim_log %>% filter(sim_run <= 5)

no_cores <- parallel::detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

t_start <- Sys.time()

all_results <- foreach(
  i = 1:nrow(sim_log),
  .combine = bind_rows,
  .packages = c("ipw", "survival", "splines", "mets", "dplyr", "tidyr", "readr"),
  .export = c("analyse_one_run", "extract_model_results", "prep_counting", "lag_ipw", "weighted_mcf", "compute_weights",
              "model_mcf")  
) %dopar% {
  analyse_one_run(sim_log$file_path[i], sim_log, true_estimands)
}

stopCluster(cl)

t_end <- Sys.time()
cat("Laufzeit gesamt:", round(difftime(t_end, t_start, units = "mins"), 2), "Minuten\n")

# check failed runs
failed_runs <- all_results %>% filter(!is.na(error))
if (nrow(failed_runs) > 0) {
  cat(nrow(failed_runs), "Läufe sind fehlgeschlagen:\n")
  print(failed_runs %>% dplyr::select(file_path, error) %>% distinct())
}

all_results <- all_results %>% filter(is.na(error))

write_csv(all_results, "./results/00_all_results_raw.csv")

# -----------------------------------------------------------------------
# Results (performance measures, Morris et al. 2019)
# -----------------------------------------------------------------------
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
    # mean number of falls until week 52 (only MCF rows)
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

write_csv(aggregated_results, "./results/00_aggregated_results.csv")

cat("\nFertig. Ergebnisse gespeichert in:\n")
cat(" - ./results/00_all_results_raw.csv (jede einzelne Schätzung)\n")
cat(" - ./results/00_aggregated_results.csv (Bias/RMSE/Coverage pro Szenario)\n")
