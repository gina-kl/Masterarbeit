###############################################################################
### Mehrfach-Simulation: Analyse über alle Szenarien und Läufe             ###
### Voraussetzung: n_sim in der Datengenerierung wurde erhöht (z.B. 100)   ###
### und alle Odat_scen*_run*.rds Dateien liegen in ./results/             ###
###############################################################################

library(ipw)
library(tidyverse)
library(magrittr)
library(survival)
library(splines)     # NEU: für ns() im Gewichtsmodell
library(foreach)
library(doParallel)

source("functions.R")

# -----------------------------------------------------------------------
# 1) Alles, was vorher in Section 1-5 des Einzel-Analyseskripts stand,
#    wird hier in eine Funktion gepackt, die EINEN Lauf verarbeitet.
# -----------------------------------------------------------------------

analyse_one_run <- function(file_path, sim_log, true_estimands) {
  
  Odat <- readRDS(file_path)
  
  current_params <- sim_log %>% filter(file_path == !!file_path)
  
  true_beta_dgp <- current_params$true_beta_arm[1]
  
  true_beta_marginal <- true_estimands %>%
    filter(
      Scenario_Name == current_params$scenario_name[1],
      Mortality     == current_params$mortality_level[1],
      True_Beta_Arm == current_params$true_beta_arm[1]
    ) %>%
    pull(True_LogRateRatio)
  
  if (length(true_beta_marginal) != 1) {
    return(tibble(file_path = file_path, error = "Kein eindeutiger Match in true_estimands"))
  }
  
  # NEU: wahrer Wert ohne Tod
  true_beta_nodeath <- true_estimands %>%
    filter(
      Scenario_Name == current_params$scenario_name[1],
      Mortality     == current_params$mortality_level[1],
      True_Beta_Arm == current_params$true_beta_arm[1]
    ) %>%
    pull(True_LogRateRatio_NoDeath)
  
  # --- Section 1: Data Preparation (Zensur-Cutoff ist bereits in der
  #     Datengenerierung korrigiert; dieser Filter ist nur eine zusaetzliche
  #     Absicherung, falls doch noch Reste vorhanden sind) ---
  Odat <- Odat %>%
    group_by(ID) %>%
    mutate(cutoff_t = pmin(dropout_time, death_time, real_followtime, na.rm = TRUE)) %>%
    filter(t <= cutoff_t) %>%
    ungroup() %>%
    dplyr::select(-cutoff_t)
  
  Odat_prepared <- Odat %>%
    group_by(ID) %>%
    mutate(
      ttl_events_pts = sum(Yobs, na.rm = TRUE),
      ttl_time_pts = max(t)
    ) %>%
    ungroup()
  
  # --- Section 3: Weights ---
  Odat_counting <- Odat_prepared %>%
    group_by(ID) %>%
    filter(t != 0) %>%
    mutate(
      tstart = t - 1,
      tstop = t,
      dropout_event = ifelse(!is.na(dropout_time) & dropout_time == tstop, 1, 0),
      cum_events = lag(cumsum(Yrep)),              # GEÄNDERT: berichtete Stürze
      cum_events = replace_na(cum_events, 0),
      CumTrain_lag = lag(CumTrain),
      CumTrain_lag = replace_na(CumTrain_lag, 0),
      CumContact_lag = lag(CumContact),                     # NEU: angenommene Anrufe bis t-1 (nur Kontrollarm)
      CumContact_lag = replace_na(CumContact_lag, 0),
      fall_rate_lag    = cum_events     / pmax(tstart, 1),  # NEU: Raten statt Anzahlen
      train_rate_lag   = CumTrain_lag   / pmax(tstart, 1),  # NEU
      contact_rate_lag = CumContact_lag / pmax(tstart, 1)   # NEU
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  Odat_arm0 <- subset(Odat_counting, arm == 0)
  Odat_arm1 <- subset(Odat_counting, arm == 1)
  
  # GEÄNDERT: drei Varianten von Gewichten (Erklärung siehe analyse.R, Section 3)
  #  falls = Sturzhistorie, oracle = wahre Frailty,
  #  proxy = Sturzhistorie + Proxy (Trainingsarm: Training, Kontrollarm: angenommene Anrufe)
  ipw_arm0_falls <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) * fall_rate_lag,
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm0
    ), error = function(e) NULL
  )
  
  ipw_arm0_proxy <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) * (fall_rate_lag + contact_rate_lag),
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm0
    ), error = function(e) NULL
  )
  
  ipw_arm0_oracle <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) + frailty,
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm0
    ), error = function(e) NULL
  )
  
  ipw_arm1_falls <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) * fall_rate_lag,
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm1
    ), error = function(e) NULL
  )
  
  ipw_arm1_proxy <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) * (fall_rate_lag + train_rate_lag),
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm1
    ), error = function(e) NULL
  )
  
  ipw_arm1_oracle <- tryCatch(
    ipwtm(
      exposure = dropout_event, family = "binomial", link = "logit",
      numerator = ~ age + sex + ns(tstop, df = 3),
      denominator = ~ age + sex + ns(tstop, df = 3) + frailty,
      id = ID, tstart = tstart, timevar = tstop, type = "cens", data = Odat_arm1
    ), error = function(e) NULL
  )
  
  if (is.null(ipw_arm0_falls) || is.null(ipw_arm0_proxy) || is.null(ipw_arm0_oracle) ||
      is.null(ipw_arm1_falls) || is.null(ipw_arm1_proxy) || is.null(ipw_arm1_oracle)) {
    return(tibble(file_path = file_path, error = "ipwtm ist fehlgeschlagen"))
  }
  
  Odat_arm0$w_falls       <- ipw_arm0_falls$ipw.weights
  Odat_arm0$w_proxy       <- ipw_arm0_proxy$ipw.weights
  Odat_arm0$w_oracle      <- ipw_arm0_oracle$ipw.weights
  
  Odat_arm1$w_falls       <- ipw_arm1_falls$ipw.weights
  Odat_arm1$w_proxy       <- ipw_arm1_proxy$ipw.weights
  Odat_arm1$w_oracle      <- ipw_arm1_oracle$ipw.weights
  
  Odat_counting_weights <- bind_rows(Odat_arm0, Odat_arm1) %>%
    arrange(ID, tstop) %>%
    as.data.frame() %>%
    mutate(status_cox = ifelse(status == 1, 1, 0))
  
  max_weight <- max(Odat_counting_weights$w_proxy, na.rm = TRUE)   # GEÄNDERT: Proxy-Gewichte
  
  # --- Section 4: Models ---
  # NEU: Woche des Dropouts/Todes (status 2/3) weglassen, dort ist kein Sturz möglich
  Odat_lwyy <- Odat_counting_weights %>% filter(status %in% c(0, 1))
  
  # GEÄNDERT: ties = "breslow" (viele gleichzeitige Stürze pro Woche) und weitere Modelle
  cox_unweighted <- tryCatch(
    coxph(Surv(tstart, tstop, status_cox) ~ arm + age + sex + cluster(ID),
          data = Odat_lwyy, ties = "breslow", robust = TRUE),
    error = function(e) NULL
  )
  
  cox_falls <- tryCatch(
    coxph(Surv(tstart, tstop, status_cox) ~ arm + age + sex + cluster(ID),
          data = Odat_lwyy, weights = w_falls, ties = "breslow", robust = TRUE),
    error = function(e) NULL
  )
  
  cox_weighted <- tryCatch(
    coxph(Surv(tstart, tstop, status_cox) ~ arm + age + sex + cluster(ID),
          data = Odat_lwyy, weights = w_proxy, ties = "breslow", robust = TRUE),
    error = function(e) NULL
  )
  
  
  cox_oracle <- tryCatch(
    coxph(Surv(tstart, tstop, status_cox) ~ arm + age + sex + cluster(ID),
          data = Odat_lwyy, weights = w_oracle, ties = "breslow", robust = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(cox_unweighted) || is.null(cox_falls) || is.null(cox_weighted) ||
      is.null(cox_oracle)) {
    return(tibble(file_path = file_path, error = "coxph ist fehlgeschlagen"))
  }
  
  # --- Section 5: Extract ---
  bind_rows(
    extract_model_results(cox_unweighted,  "1: Naive (Unweighted)",      true_beta_marginal, true_beta_nodeath),
    extract_model_results(cox_falls,       "2: LWYY + IPCW falls",       true_beta_marginal, true_beta_nodeath),
    extract_model_results(cox_weighted,    "3: LWYY + IPCW proxy",       true_beta_marginal, true_beta_nodeath),
    extract_model_results(cox_oracle,      "4: LWYY + IPCW oracle",      true_beta_marginal, true_beta_nodeath)
  ) %>%
    mutate(
      file_path = file_path,
      scenario_id = current_params$scenario_id[1],
      scenario_name = current_params$scenario_name[1],
      sim_run = current_params$sim_run[1],
      true_beta_arm_dgp = true_beta_dgp,
      true_beta_marginal = true_beta_marginal,
      drop_beta_arm = current_params$drop_beta_arm[1],
      mortality_level = current_params$mortality_level[1],
      report_prob = current_params$report_prob[1],     # NEU
      true_beta_nodeath = true_beta_nodeath,            # NEU
      max_weight = max_weight,
      error = NA_character_
    )
}
# -----------------------------------------------------------------------
# 2) Über alle Läufe loopen (parallelisiert)
# -----------------------------------------------------------------------

sim_log <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)

# time check
# GEÄNDERT: nach unten verschoben -- vorher wurde sim_log benutzt, bevor es eingelesen war
t1 <- Sys.time()
test_result <- analyse_one_run(sim_log$file_path[1], sim_log, true_estimands)
print(Sys.time() - t1)

# Optional zum Testen: erstmal nur eine Teilmenge laufen lassen, z.B.:
# sim_log <- sim_log %>% filter(sim_run <= 5)

no_cores <- parallel::detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

t_start <- Sys.time()

all_results <- foreach(
  i = 1:nrow(sim_log),
  .combine = bind_rows,
  .packages = c("ipw", "survival", "splines", "dplyr", "tidyr", "readr"),   # GEÄNDERT: splines
  .export = c("analyse_one_run", "extract_model_results")                  # NEU: Funktionen an die Worker übergeben
) %dopar% {
  analyse_one_run(sim_log$file_path[i], sim_log, true_estimands)
}

stopCluster(cl)

t_end <- Sys.time()
cat("Laufzeit gesamt:", round(difftime(t_end, t_start, units = "mins"), 2), "Minuten\n")

# Fehlgeschlagene Läufe prüfen
failed_runs <- all_results %>% filter(!is.na(error))
if (nrow(failed_runs) > 0) {
  cat(nrow(failed_runs), "Läufe sind fehlgeschlagen:\n")
  print(failed_runs %>% dplyr::select(file_path, error) %>% distinct())
}

all_results <- all_results %>% filter(is.na(error))

write_csv(all_results, "./results/00_all_results_raw.csv")

# -----------------------------------------------------------------------
# 3) Aggregation pro Szenario: Bias, empirische SE, RMSE, Coverage, Power
# -----------------------------------------------------------------------

aggregated_results <- all_results %>%
  group_by(scenario_id, scenario_name, drop_beta_arm, report_prob, mortality_level,          # GEÄNDERT: report_prob
           true_beta_arm_dgp, true_beta_marginal, true_beta_nodeath, Model) %>%          # GEÄNDERT: true_beta_nodeath
  summarise(
    n_runs = n(),
    mean_max_weight = mean(max_weight, na.rm = TRUE),
    mean_LogHR = mean(Log_HR, na.rm = TRUE),
    mean_bias = mean(Bias, na.rm = TRUE),
    mean_bias_nodeath = mean(Bias_NoDeath, na.rm = TRUE),               # NEU: Bias gegenüber wahrem Wert ohne Tod
    mcse_bias = sd(Log_HR, na.rm = TRUE) / sqrt(n()),                   # NEU: Monte-Carlo-Standardfehler des Bias
    empirical_SE = sd(Log_HR, na.rm = TRUE),
    mean_model_SE = mean(Robust_SE, na.rm = TRUE),
    RMSE = sqrt(mean(Bias^2, na.rm = TRUE)),
    coverage_95 = mean(
      (Log_HR - 1.96 * Robust_SE <= true_beta_marginal) &
        (Log_HR + 1.96 * Robust_SE >= true_beta_marginal),
      na.rm = TRUE
    ),
    coverage_95_nodeath = mean(                                          # NEU: gegenüber wahrem Wert ohne Tod
      (Log_HR - 1.96 * Robust_SE <= true_beta_nodeath) &
        (Log_HR + 1.96 * Robust_SE >= true_beta_nodeath),
      na.rm = TRUE
    ),
    power_or_type1 = mean(Significant == "Yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario_id, Model)

print(aggregated_results, n = 100)

write_csv(aggregated_results, "./results/00_aggregated_results.csv")

cat("\nFertig. Ergebnisse gespeichert in:\n")
cat(" - ./results/00_all_results_raw.csv (jede einzelne Schätzung)\n")
cat(" - ./results/00_aggregated_results.csv (Bias/RMSE/Coverage pro Szenario)\n")