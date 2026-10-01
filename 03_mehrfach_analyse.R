###############################################################################
### Mehrfach-Simulation: Analyse über alle Szenarien und Läufe             ###
### Voraussetzung: n_sim in der Datengenerierung wurde erhöht (z.B. 100)   ###
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

source("functions.R")

# -----------------------------------------------------------------------
# all of section 1-5 in one function
# -----------------------------------------------------------------------

analyse_one_run <- function(file_path, sim_log, true_estimands,
                            trunc_level = 0.01) {   # IPW truncation at 1st / 99th percentile
  
  Odat_all <- readRDS(file_path)   

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
  
  # true beta without death
  true_beta_nodeath <- true_estimands %>%
    filter(
      Scenario_Name == current_params$scenario_name[1],
      Mortality     == current_params$mortality_level[1],
      True_Beta_Arm == current_params$true_beta_arm[1]
    ) %>%
    pull(True_LogRateRatio_NoDeath)

  # true mean number of falls without death (for LWYY)
  true_mcf <- true_estimands %>%
    filter(
      Scenario_Name == current_params$scenario_name[1],
      Mortality     == current_params$mortality_level[1],
      True_Beta_Arm == current_params$true_beta_arm[1]
    )

  # method depends on mortality: LWYY only for low, Ghosh-Lin only for high mortality
  method <- ifelse(current_params$mortality_level[1] == "low", "LWYY", "GL")

  # ---Data Preparation ---
  # observed data in the study 
  Odat <- Odat_all %>% filter(observed)

  # complete data for oracle M0 (no dropout, no random censoring, only death / t.end)
  Odat_full <- Odat_all %>%
    filter(t != 0) %>%
    mutate(tstart = t - 1, tstop = t) %>%
    as.data.frame()

  Odat_prepared <- Odat %>%
    group_by(ID) %>%
    mutate(
      ttl_events_pts = sum(Yobs, na.rm = TRUE),
      ttl_time_pts = max(t)
    ) %>%
    ungroup()
  
  # --- Weights ---
  Odat_counting <- Odat_prepared %>%
    group_by(ID) %>%
    filter(t != 0) %>%
    mutate(
      tstart = t - 1,
      tstop = t,
      dropout_event = ifelse(!is.na(dropout_time) & dropout_time == tstop, 1, 0),
      cum_events = lag(cumsum(Yrep)),              # reported falls
      cum_events = replace_na(cum_events, 0),
      CumTrain_lag = lag(CumTrain),
      CumTrain_lag = replace_na(CumTrain_lag, 0),
      CumContact_lag = lag(CumContact),                     # control arm
      CumContact_lag = replace_na(CumContact_lag, 0),
      fall_rate_lag    = cum_events     / pmax(tstart, 1),   
      train_rate_lag   = CumTrain_lag   / pmax(tstart, 1),   
      contact_rate_lag = CumContact_lag / pmax(tstart, 1)   
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  Odat_arm0 <- subset(Odat_counting, arm == 0)
  Odat_arm1 <- subset(Odat_counting, arm == 1)
  
  # three different proxys
  #  oracle = with frailty
  #  falls = fall history, oracle = true Frailty,
  #  proxy = fall history + Proxy (treatment arm: training, control arm: control calls)
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
  
  Odat_arm0$w_falls       <- trunc_weights(ipw_arm0_falls$ipw.weights, trunc_level)
  Odat_arm0$w_proxy       <- trunc_weights(ipw_arm0_proxy$ipw.weights, trunc_level)
  Odat_arm0$w_oracle      <- trunc_weights(ipw_arm0_oracle$ipw.weights, trunc_level)
  
  Odat_arm1$w_falls       <- trunc_weights(ipw_arm1_falls$ipw.weights, trunc_level)
  Odat_arm1$w_proxy       <- trunc_weights(ipw_arm1_proxy$ipw.weights, trunc_level)
  Odat_arm1$w_oracle      <- trunc_weights(ipw_arm1_oracle$ipw.weights, trunc_level)
  
  Odat_counting_weights <- bind_rows(Odat_arm0, Odat_arm1) %>%
    arrange(ID, tstop) %>%
    as.data.frame()

  max_weight <- max(Odat_counting_weights$w_proxy, na.rm = TRUE)   # proxy weights (truncated)
  # untruncated proxy weights (to see how extreme they were before truncation)
  max_weight_untrunc <- max(ipw_arm0_proxy$ipw.weights, ipw_arm1_proxy$ipw.weights, na.rm = TRUE)

  # --- Models ---
  fits_ok <- tryCatch({

    if (method == "LWYY") {
      set.seed(current_params$seed[1])   # seed for split event times
      
      Odat_lwyy <- Odat_counting_weights %>% filter(status %in% c(0, 1)) %>%
        expand_falls("status", "Yrep")
      Odat_lwyy_full <- Odat_full %>% filter(status_full %in% c(0, 1)) %>%
        expand_falls("status_full", "Yfull")

      fit_full <- coxph(Surv(tstart, tstop, status_full) ~ arm + age + sex + cluster(ID),
                        data = Odat_lwyy_full, ties = "breslow", robust = TRUE)
      fit_naive <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                         data = Odat_lwyy, ties = "breslow", robust = TRUE)
      fit_falls <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                         data = Odat_lwyy, weights = w_falls, ties = "breslow", robust = TRUE)
      fit_weighted <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                            data = Odat_lwyy, weights = w_proxy, ties = "breslow", robust = TRUE)
      fit_oracle <- coxph(Surv(tstart, tstop, status) ~ arm + age + sex + cluster(ID),
                          data = Odat_lwyy, weights = w_oracle, ties = "breslow", robust = TRUE)
    }

    if (method == "GL") {
      set.seed(current_params$seed[1])
      Odat_gl <- prep_gl_data(expand_falls(Odat_counting_weights, "status", "Yrep"),
                                   status_col = "status")
      Odat_gl_full <- prep_gl_data(expand_falls(Odat_full, "status_full", "Yfull"),
                                   status_col = "status_full")

      fit_full <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                         data = Odat_gl_full, cause = 1, death.code = 3, cens.code = 0)
      fit_naive <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                          data = Odat_gl, cause = 1, death.code = 3, cens.code = 0)
      fit_falls <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                          data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                          weights = Odat_gl$w_falls)
      fit_weighted <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                             data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                             weights = Odat_gl$w_proxy)
      fit_oracle <- recreg(Event(tstart, tstop, status_gl) ~ arm + age + sex + cluster(ID),
                           data = Odat_gl, cause = 1, death.code = 3, cens.code = 0,
                           weights = Odat_gl$w_oracle)
    }
    TRUE
  }, error = function(e) conditionMessage(e))
  
  if (!isTRUE(fits_ok)) {
    return(tibble(file_path = file_path, error = paste(method, "ist fehlgeschlagen:", fits_ok)))
  }

  # --- Mean number of falls until week 52 (G-computation, only LWYY) ---
  # survfit needs a coxph model -> for Ghosh-Lin NA
  fits <- list("0: Full data (oracle M0)" = fit_full,
               "1: Naive (Unweighted)"    = fit_naive,
               "2: IPCW falls"            = fit_falls,
               "3: IPCW proxy"            = fit_weighted,
               "4: IPCW oracle"           = fit_oracle)

  if (method == "LWYY") {
    # one row per patient
    baseline_pts <- Odat_counting_weights %>% distinct(ID, age, sex)
    exposed   <- baseline_pts %>% mutate(arm = 1)     # Training group
    unexposed <- baseline_pts %>% mutate(arm = 0)     # Control group

    mean_falls <- bind_rows(lapply(names(fits), function(m) {
      pred_exposed   <- survfit(fits[[m]], newdata = exposed, se.fit = FALSE)  # se.fit False for faster time
      pred_unexposed <- survfit(fits[[m]], newdata = unexposed, se.fit = FALSE)

      data.frame(
        time           = pred_exposed$time,
        mean_exposed   = rowMeans(pred_exposed$cumhaz),
        mean_unexposed = rowMeans(pred_unexposed$cumhaz)
      ) %>%
        filter(time <= 52) %>%
        tail(1) %>%
        transmute(Model = paste(method, m), mean_exposed, mean_unexposed)
    }))
  } else {
    mean_falls <- tibble(Model = paste(method, names(fits)),
                         mean_exposed = NA_real_, mean_unexposed = NA_real_)
  }

  # --- Extract ---
  bind_rows(
    extract_model_results(fit_full,     paste(method, "0: Full data (oracle M0)"), true_beta_marginal, true_beta_nodeath),
    extract_model_results(fit_naive,    paste(method, "1: Naive (Unweighted)"),    true_beta_marginal, true_beta_nodeath),
    extract_model_results(fit_falls,    paste(method, "2: IPCW falls"),            true_beta_marginal, true_beta_nodeath),
    extract_model_results(fit_weighted, paste(method, "3: IPCW proxy"),            true_beta_marginal, true_beta_nodeath),
    extract_model_results(fit_oracle,   paste(method, "4: IPCW oracle"),           true_beta_marginal, true_beta_nodeath)
  ) %>%
    left_join(mean_falls, by = "Model") %>%
    mutate(
      Method = method,
      true_mcf_exposed_nodeath   = true_mcf$MCF_Intervention_NoDeath,
      true_mcf_unexposed_nodeath = true_mcf$MCF_Control_NoDeath,
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
      max_weight_untrunc = max_weight_untrunc,
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
  .export = c("analyse_one_run", "extract_model_results", "prep_gl_data", "expand_falls", "trunc_weights")      
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
# Bias, empirische SE, RMSE, Coverage, Power
# -----------------------------------------------------------------------

aggregated_results <- all_results %>%
  group_by(scenario_id, scenario_name, drop_beta_arm, report_prob, mortality_level,          
           true_beta_arm_dgp, true_beta_marginal, true_beta_nodeath, Method, Model) %>%  
  summarise(
    n_runs = n(),
    mean_max_weight = mean(max_weight, na.rm = TRUE),
    mean_max_weight_untrunc = mean(max_weight_untrunc, na.rm = TRUE),
    mean_LogHR = mean(Log_HR, na.rm = TRUE),
    mean_bias = mean(Bias, na.rm = TRUE),
    mean_bias_nodeath = mean(Bias_NoDeath, na.rm = TRUE),               
    mcse_bias = sd(Log_HR, na.rm = TRUE) / sqrt(n()),                   
    empirical_SE = sd(Log_HR, na.rm = TRUE),
    mean_model_SE = mean(Robust_SE, na.rm = TRUE),
    RMSE = sqrt(mean(Bias^2, na.rm = TRUE)),
    coverage_95 = mean(
      (Log_HR - 1.96 * Robust_SE <= true_beta_marginal) &
        (Log_HR + 1.96 * Robust_SE >= true_beta_marginal),
      na.rm = TRUE
    ),
    coverage_95_nodeath = mean(                                         
      (Log_HR - 1.96 * Robust_SE <= true_beta_nodeath) &
        (Log_HR + 1.96 * Robust_SE >= true_beta_nodeath),
      na.rm = TRUE
    ),
    power_or_type1 = mean(Significant == "Yes", na.rm = TRUE),
    # mean number of falls until week 52 (only LWYY, NA for GL)
    true_mcf_exposed_nodeath   = first(true_mcf_exposed_nodeath),
    true_mcf_unexposed_nodeath = first(true_mcf_unexposed_nodeath),
    mean_falls_exposed   = mean(mean_exposed, na.rm = TRUE),
    mean_falls_unexposed = mean(mean_unexposed, na.rm = TRUE),
    bias_falls_exposed   = mean_falls_exposed   - true_mcf_exposed_nodeath,
    bias_falls_unexposed = mean_falls_unexposed - true_mcf_unexposed_nodeath,
    .groups = "drop"
  ) %>%
  arrange(scenario_id, Model)

print(aggregated_results, n = 100)

write_csv(aggregated_results, "./results/00_aggregated_results.csv")

cat("\nFertig. Ergebnisse gespeichert in:\n")
cat(" - ./results/00_all_results_raw.csv (jede einzelne Schätzung)\n")
cat(" - ./results/00_aggregated_results.csv (Bias/RMSE/Coverage pro Szenario)\n")