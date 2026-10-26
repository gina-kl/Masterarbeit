
##########################################################################################
#------- Data generation ---------------------------------------------------------------##
##########################################################################################
run_on_cluster <- TRUE
source("cluster_setup.R")   # run_on_cluster, n_cores_use, results_dir, part
source("functions.R")
source("GlobalVars.R")

# Load packages, Rscript on the cluster starts without any packages
library(simcausal)
# single packages instead of tidyverse (tidyverse cannot be installed on the cluster)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(magrittr)
library(foreach)
library(doParallel)
library(survival)   # true values: LWYY
library(mets)       # true values: Ghosh-Lin

#-------Fixed coefficients----------------------------------------------------------------##

COEF <- list(
  age_cent = 72.5,     # centered age (mid age)
  sex_cent = 0.5,      # centered sex 
  y_age = 0.02,        # effect age on falls
  y_sex = 0.1,         # effect sex on falls
  cens_age = 0.03,     # effect age on dropout (older patient drop out more often)
  cens_sex = 0.1,      # effect sex on dropout 
  d_age = 0.05,        # effect age on death
  d_sex = 0.1,         # effect sex on death
  train_age = -0.02,   # effect age on training (older patient train less more often)
  train_sex = 0.3      # >0: men train more
)


# Use one DAG for simulated and true values
# arm_value = NA -> randomised
# arm_value = 0/1: fixed arm (large population for the true values)
# Dropout is not part of the DAG -> generated afterwards

sim_dag <- function(n, seed, arm_value = NA,
                    t.end = 52,                      # weeks (one year)
                    true_beta_arm = -0.4,            # direct (conditional) effect of training on falls
                    frailty_sd = 1.0,                # sd of frailty
                    Y_intercept = -3.5,              # calibrated (calib_Y)
                    death_intercept = -Inf,          # calibrated (calib_logit), -Inf = no death
                    train_intercept = 0.55,          # calibrated (calib_logit, week 1)
                    train_frailty_eff = 0.02,        # frailty reduction per training session
                    max_frailty_red = 0.5,           # maximal frailty reduction
                    eff_frailty_y = 0.4,             # effect frailty on falls
                    eff_frailty_death = 0.6,         # effect frailty on death
                    eff_frailty_train = -0.8         # effect frailty on training (<0: frail train less)
                   ) {
  # fixed coefficients as local cov
  age_cent <- COEF$age_cent
  sex_cent <- COEF$sex_cent
  y_age <- COEF$y_age
  y_sex <- COEF$y_sex
  d_age <- COEF$d_age
  d_sex <- COEF$d_sex
  train_age <- COEF$train_age
  train_sex <- COEF$train_sex
  
  D <- DAG.empty()
  
  
  #--- Baseline (t = 0) --------------------------------------------------------------------
  D <- D +
    node("Frailty", t = 0, distr = "rnorm", mean = 0, sd = frailty_sd) +   # unobserved frailty
    node("age", t = 0, distr = "runif", min = 65, max = 80) +
    node("sex", t = 0, distr = "rbern", prob = 0.5) +                     # 0 = female, 1 = male
    node("L1", t = 0, distr = "rpois",                                    # proxy: falls in the year before the study
         lambda = 52 * exp(Y_intercept + y_age * (age[0] - age_cent) + y_sex * (sex[0] - sex_cent) + 
                             eff_frailty_y * Frailty[0]))                 # baseline frailty value
  if (is.na(arm_value)) {
    D <- D + node("A1", t = 0, distr = "rbern", prob = 0.5)              # randomisation (trial)
  } else {
    D <- D + node("A1", t = 0, distr = "rconst", const = arm_value)      # fixed arm (true values)
  }
  D <- D +
    node("FrailtyT", t = 0, distr = "rconst", const = Frailty[0]) +      # time varying frailty
    node("Train",    t = 0, distr = "rconst", const = 0) +               # training 
    node("CumTrain", t = 0, distr = "rconst", const = 0) +               # number of trainings
    node("Yfull",    t = 0, distr = "rconst", const = 0) +               # number of falls
    node("Death",    t = 0, distr = "rconst", const = 0)                 

  #--- t>0 --------------------------------------------------------------------------------------------
  D <- D +
    
    # time-varying frailty, decreases with cumulated training (only training arm trains)
    node("FrailtyT", t = 1:t.end, distr = "rconst",
         const = Frailty[0] - ifelse(train_frailty_eff * CumTrain[t-1] > max_frailty_red,
                                     max_frailty_red, train_frailty_eff * CumTrain[t-1])) +
    
    # training participation (driven by age, sex, current frailty);
    # continues after dropout (needed for the complete data Yfull), stops at death
    node("Train", t = 1:t.end, distr = "rbern",
         prob = ifelse(A1[0] == 0 | Death[t-1] == 1, 0,
                       plogis(train_intercept + train_age * (age[0] - age_cent) +
                                train_sex * (sex[0] - sex_cent) + eff_frailty_train * FrailtyT[t]))) +
    
    node("CumTrain", t = 1:t.end, distr = "rconst", const = CumTrain[t-1] + Train[t]) +
    
    # falls in week t if alive at the start of week t (death happens at the end)
    node("Yfull", t = 1:t.end, distr = "rpois",
         lambda = ifelse(Death[t-1] == 1, 0,
                         exp(Y_intercept + y_age * (age[0] - age_cent) + y_sex * (sex[0] - sex_cent) +
                               true_beta_arm * A1[0] + eff_frailty_y * FrailtyT[t]))) +
    
    # death at the end of week t
    node("Death", t = 1:t.end, distr = "rbern",
         prob = ifelse(Death[t-1] == 1, 1,
                       plogis(death_intercept + d_age * (age[0] - age_cent) + d_sex * (sex[0] - sex_cent) +
                                eff_frailty_death * FrailtyT[t])))
   

  sim(DAG = set.DAG(D), n = n, wide = FALSE, rndseed = seed)
}


# dropout afterwards, from its own random stream
# same uniforms U[person, week] for all dropout mechanisms (only depend on seed_drop) 
# scenarios with the same causal parameters have identical complete data and coupled dropout
# Dropout[t] = 1 if Dropout[t-1] = 1; 0 if death at the end of week t (or before); else Bernoulli(p_t)
# Yobs[t]    = Yfull[t] * (1 - Dropout[t-1])  (falls of the dropout week are observed)
# eff_fall_drop: dropout depends on log(1 + falls up to t-1) (= log_falls_lag in the weight models). While a
# person is in the study, the observed falls are the complete falls -> computed from Yfull
add_dropout <- function(dat, seed_drop, t.end = 52,
                        drop_intercept, drop_beta_arm = 0,
                        eff_frailty_drop = 0, eff_L1_drop = 0, eff_fall_drop = 0) {
  set.seed(seed_drop)
  ids <- sort(unique(dat$ID))
  U   <- matrix(runif(length(ids) * t.end), nrow = length(ids))   # U[person, week]

  dat %>%
    arrange(ID, t) %>%
    group_by(ID) %>%
    mutate(arm_i = A1[1], age_i = age[1], sex_i = sex[1], L1_i = L1[1],   # baseline values (t = 0)
           log_falls_lag = log1p(cumsum(Yfull) - Yfull)) %>%            # log(1 + falls up to t-1)
    ungroup() %>%
    mutate(
      p_drop   = plogis(drop_intercept +
                          COEF$cens_age * (age_i - COEF$age_cent) + COEF$cens_sex * (sex_i - COEF$sex_cent) +
                          eff_frailty_drop * FrailtyT + eff_L1_drop * log1p(L1_i) +   # log1p computes log(1+x) also für x<<1
                          eff_fall_drop * log_falls_lag +                             
                          drop_beta_arm * arm_i),
      u        = ifelse(t >= 1, U[cbind(match(ID, ids), pmax(t, 1))], 1),
      drop_now = t >= 1 & Death == 0 & u < p_drop          # dropout at the end of week t, not in the death week
    ) %>%
    group_by(ID) %>%
    mutate(Dropout = as.integer(cumsum(drop_now) > 0),      # absorbing
           Yobs    = Yfull * (1 - lag(Dropout, default = 0L))) %>%
    ungroup() %>%
    dplyr::select(-arm_i, -age_i, -sex_i, -L1_i, -log_falls_lag, -p_drop, -u, -drop_now)   
}


# Function to simulate data
simulate_scenario <- function(
    seed = 123,                      # seed of the complete data (same for all scenarios with the same causal parameters)
    seed_drop = seed + 5e8,          # seed of the dropout uniforms (own stream)
    n_pts = 1800,                   # number of patients (~900 per arm)
    t.end = 52,                      # maximal time (one year in weeks)
    true_beta_arm = -0.4,            # true direct treatment effect on recurrent falls
    frailty_sd = 1.0,                # sd of Frailty
    drop_intercept = -5.3,           # Baseline dropout per week (calibrated per scenario with calib_logit)
    drop_beta_arm = 0.0,             # extra arm effect on dropout (<0: training group drops out less), with frailty effects dropout differs between arms also for 0
    train_intercept = 0.55,          # Baseline train rate (calibrated per scenario with calib_logit: 63 % train in week 1)
    Y_intercept = -3.5,              # Baseline fall rate per week (calibrated per scenario with calib_Y)
    death_intercept = -Inf,          # Intercept for mortality, calibrated per scenario with calib_logit (-Inf = no death)
    # Parameter for scenarios
    train_frailty_eff = 0.02,  # reduction of frailty per training session (0 = no effect)
    max_frailty_red   = 0.5,   # maximal reduction of frailty (stabilizes afterwards)
    eff_frailty_y = 0.4,       # Effect frailty on falls
    eff_frailty_death = 0.6,    # Effect frailty on death
    eff_frailty_drop = 0.6,    # Effect frailty on dropout (informative censoring)
    eff_frailty_train = -0.8,  # Effect frailty on training (<0 : train less >0: train more)
    eff_L1_drop = 0.0,         # Effect observed baseline falls on dropout
    eff_fall_drop = 0.0        # Effect observed fall rate up to t-1 on dropout (positive control, should be corrected by ipcw)

) {
  
  # complete data (DAG without dropout, defined in sim_dag())
  Odat <- sim_dag(n = n_pts, seed = seed, t.end = t.end,
                  true_beta_arm = true_beta_arm, frailty_sd = frailty_sd,
                  Y_intercept = Y_intercept, death_intercept = death_intercept,
                  train_intercept = train_intercept,
                  train_frailty_eff = train_frailty_eff, max_frailty_red = max_frailty_red,
                  eff_frailty_y = eff_frailty_y, eff_frailty_death = eff_frailty_death,
                  eff_frailty_train = eff_frailty_train)

  # dropout 
  Odat <- add_dropout(Odat, seed_drop = seed_drop, t.end = t.end,
                      drop_intercept = drop_intercept, drop_beta_arm = drop_beta_arm,
                      eff_frailty_drop = eff_frailty_drop, eff_L1_drop = eff_L1_drop,
                      eff_fall_drop = eff_fall_drop)  

  # check for NA
  if (anyNA(Odat$Dropout) | anyNA(Odat$Yobs)) {
    stop("NA in Dropout oder Yobs: check parameter")
  }
  
  Odat %<>%
    arrange(ID, t) %>%  # arrange for ID and t
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = age[1],
      sex = sex[1],
      base_risk = L1[1],             # proxy 1: falls in the year before the study
      frailty = Frailty[1],
      frailty_t = FrailtyT,          # used for oracle
      
      # dropout and death time to compute full data and dropout data (in analyse.R)
      temp_drop_time = ifelse(Dropout == 1, t, NA),  # if droped out = time, instead NA
      dropout_time = ifelse(!all(is.na(temp_drop_time)), min(temp_drop_time, na.rm = TRUE), NA),
      
      temp_death_time = ifelse(Death == 1, t, NA),
      death_time = ifelse(!all(is.na(temp_death_time)), min(temp_death_time, na.rm = TRUE), NA),
      
      ever_died = as.integer(!is.na(death_time)),
      
      # proxies, known at the start of week t
      cum_train_lag = lag(CumTrain, default = 0),      # trainings up to week t-1
      cum_falls_lag = lag(cumsum(Yobs), default = 0)   # observed falls up to week t-1
    ) %>%
    # complete data until death_time (falls of the death week count)
    filter(t >= 1, is.na(death_time) | t <= death_time) %>%
    mutate(
      # TRUE = alive and in study at the start of week t
      # (dropout at the END of week dropout_time -> that week is still observed)
      observed = is.na(dropout_time) | t <= dropout_time,
      
      # number of falls (several falls per week possible)
      n_falls      = ifelse(observed, Yobs, NA),   # observed falls
      n_falls_full = Yfull,                        # complete data (oracle M0)
      
      # separate indicators: death at the end of the week and falls in that week are possible
      drop       = as.integer(observed & !is.na(dropout_time) & t == dropout_time),
      death      = as.integer(observed & !is.na(death_time) & t == death_time),
      death_full = as.integer(!is.na(death_time) & t == death_time),
      
      train_rate_lag = ifelse(t > 1, cum_train_lag / (t - 1), 0),   # share of weeks trained (until t-1)
      
      # observed data: 0 = no event, 1 = fall, 2 = Dropout, 3 = death, NA = not observed
      status = case_when(
        !observed ~ NA_real_,  # not observed (after dropout)
        death == 1 ~ 3,        # death at the end of this week
        drop == 1 ~ 2,         # dropout at the end of this week
        TRUE ~ as.numeric(n_falls > 0)  # at least one fall
      ),
      # complete data (for oracle M0): 0 = no event, 1 = fall, 3 = death (same note)
      status_full = case_when(
        death_full == 1 ~ 3,
        TRUE ~ as.numeric(Yfull > 0)
      )
    ) %>%
    group_by(ID) %>%
    mutate(ever_dropped = as.integer(any(drop == 1))) %>%
    ungroup() %>%
    select(-temp_drop_time, -temp_death_time)    # drop helper variables
  
  return(Odat)
}





#-------Simulate true values ------------------------------------------------------------##
###### Monte Carlo Simulation of the true values ##########################################


# simulate a very large trial per scenario without dropout (complete data, only death stops falls). Then
#   - mean number of falls per arm  -> True_LogRateRatio
#   - fit the same models as in the analysis (LWYY, Ghosh-Lin, adjusted for age, sex,
#     log_base_falls on complete data -> True_Beta_LWYY, True_Beta_GL)
# The true value of a method is what it estimates with complete data and very large n.
# Both arms use the same seed (same frailty, age, sex) -> less Monte Carlo noise.

simulate_true <- function(n_pts_true = 50000,
                          seed = 123,
                          ...) {                 # all other DGP parameters -> sim_dag()
  bind_rows(
    sim_dag(n = n_pts_true, seed = seed, arm_value = 0, ...),                              # control
    sim_dag(n = n_pts_true, seed = seed, arm_value = 1, ...) %>%                           # training
      mutate(ID = ID + n_pts_true)  
  )
}


#------------------- Scenarios -----------------------------------------------------------##
# Calibration of intercepts
# The intercepts are chosen so that the rates in the whole population (with frailty) hit the targets
# coefficients from COEF -> automatically the same as in sim_dag()

# reference population with the same distribution as in the data generation (frailty_sd = 1)
set.seed(1)
cd <- tibble(age = runif(2e5, 65, 80), sex = rbinom(2e5, 1, 0.5), Z = rnorm(2e5),  # Z=frailty
             U = runif(2e5))                                                       # uniform for L1 (scenario 4)

# Y_intercept so that number of units with >= 1 fall in 52 weeks = target (for control arm, no death)
# P(N >= 1) = 1 - exp(-52 * weekly rate) (Poisson), averaged over all persons
calib_Y <- function(eff_frailty_y, target = 0.30, t.end = 52) {
  lp <- COEF$y_age * (cd$age - COEF$age_cent) + COEF$y_sex * (cd$sex - COEF$sex_cent) + eff_frailty_y * cd$Z
  uniroot(function(a) mean(1 - exp(-t.end * exp(a + lp))) - target, c(-15, 0))$root
}

# intercept for weekly logit model (dropout / death)
# so that share with event within 52 weeks = target
# P(event) = 1 - (1 - weekly prob)^52, averaged over all persons
calib_logit <- function(lp, target, t.end = 52) {
  if (target == 0) return(-Inf)                       # -Inf -> plogis = 0 -> no event
  uniroot(function(b) mean(1 - (1 - plogis(b + lp))^t.end) - target, c(-20, 5))$root
}

# training: calibrated only for week 1, because later weeks depend on the own training history
# (training -> frailty -> training). In week 1 CumTrain = 0 -> FrailtyT = Frailty, so
# calib_logit with t.end = 1 gives E[plogis(b + lp)] = target

# dropout depending on log(1 + falls up to t-1) (scenario 4): the weekly probability changes over time,
# so there is no closed formula. The weekly falls of the reference population are simulated
# (control arm: no training -> constant rate, no death as in the other calibrations), then
# P(no dropout until week 52) = prod_t (1 - plogis(b + lp + eff_fall_drop * log_falls_lag_t)), averaged
calib_drop_fall <- function(lp, eff_fall_drop, Y_intercept, eff_frailty_y, target, t.end = 52) {
  lambda <- exp(Y_intercept + COEF$y_age * (cd$age - COEF$age_cent) + COEF$y_sex * (cd$sex - COEF$sex_cent) +
                  eff_frailty_y * cd$Z)
  set.seed(2)                                                     # fixed -> same intercept in every array task
  Y <- matrix(rpois(nrow(cd) * t.end, rep(lambda, t.end)), nrow = nrow(cd))   # Y[person, week]
  cum_lag <- matrix(0, nrow(cd), t.end)                                     # falls up to t-1
  for (k in 2:t.end) cum_lag[, k] <- cum_lag[, k - 1] + Y[, k - 1]
  log_falls_lag <- log1p(cum_lag)                                           # = log_falls_lag in the analysis
  uniroot(function(b) {
    log_surv <- rowSums(log1p(-plogis(b + lp + eff_fall_drop * log_falls_lag)))  # lp recycled over the weeks
    mean(1 - exp(log_surv)) - target
  }, c(-20, 5))$root
}

# scenario definitions, effects, blocks, scenario_id and causal_id: scenarios.R
source("scenarios.R")

base_scenarios <- base_scenarios %>%
  rowwise() %>%                          # calibration for each scenario separately
  mutate(
    Y_intercept     = calib_Y(eff_frailty_y, target_falls),
    drop_intercept  = {
      lp_drop <- COEF$cens_age * (cd$age - COEF$age_cent) + COEF$cens_sex * (cd$sex - COEF$sex_cent) +
        eff_frailty_drop * cd$Z +
        # observed baseline falls (eff_L1_drop, currently 0 in all scenarios)
        eff_L1_drop * log1p(qpois(cd$U, 52 * exp(Y_intercept +COEF$y_age * (cd$age - COEF$age_cent) + COEF$y_sex * (cd$sex - COEF$sex_cent) +
                                                   eff_frailty_y * cd$Z)))
      if (eff_fall_drop == 0) calib_logit(lp_drop, target_dropout)
      else calib_drop_fall(lp_drop, eff_fall_drop, Y_intercept, eff_frailty_y, target_dropout)   # scenario 4
    },
    death_intercept = calib_logit(COEF$d_age * (cd$age - COEF$age_cent) + COEF$d_sex * (cd$sex - COEF$sex_cent) +
                                    eff_frailty_death * cd$Z, target_death[[mortality_level]]),
    train_intercept = calib_logit(COEF$train_age * (cd$age - COEF$age_cent) + COEF$train_sex * (cd$sex - COEF$sex_cent) +
                                    eff_frailty_train * cd$Z, target_train, t.end = 1)
  ) %>%
  ungroup()

message("Scenario blocks: ", paste(scenario_blocks, collapse = ", "),
        " | Number of Scenarios (base_scenarios): ", nrow(base_scenarios))   

# for run on cluster: which part should run
run_true <- part %in% c("all", "true") # run true values
run_data <- part %in% c("all", "data") # run simulated data
if (!run_true && !run_data) stop("01_generate: argument must be 'true' or 'data', not: ", part)


#------------------- True values ---------------------------------------------------------##


# The true values do not depend on the dropout mechanism (Yfull continues after dropout), they are computed once per causal scenario
# Monte Carlo in blocks of n_chunk_true patients per arm (run parallel on cluster)
#   - one parallel task = one causal scenario x one block -> memory per worker depends only on n_chunk_true
#   - per block only small summaries are kept (falls per arm and quarter, deaths, model coefficients)
#   - mean numbers of falls: pooled over all blocks, LWYY / GL: mean of the block estimates
#     (bias of a block estimate is O(1 / n_chunk_true), negligible)
#   - Monte Carlo SE of each true value = sd over the blocks / sqrt(n_blocks)
#   - seed depends only on scenario and block -> a larger n_true keeps the first blocks and adds new ones



if (run_true) {   # cluster
  # unique causal scenarios
  # causal_id from scenarios.R (complete grid) -> same id and seeds whichever blocks are run
  unique_causal_scenarios <- base_scenarios %>%
    distinct(causal_id, .keep_all = TRUE) %>%
    dplyr::select(causal_id, effect, all_of(causal_vars), Y_intercept, death_intercept, train_intercept)   # intercepts follow from causal_vars
  
  if (n_true %% n_chunk_true != 0) stop("n_true must be a multiple of n_chunk_true")
  n_blocks <- n_true / n_chunk_true
  
  true_tasks <- expand_grid(causal_id = unique_causal_scenarios$causal_id,
                            block     = seq_len(n_blocks))
  
  # one block of the true value simulation for scenario p
  # both arms with the same seed, once with and once without death
  true_block <- function(p, block, n_chunk, seed) {
    
    # only the death intercept differs between the calls, only the columns needed later are kept
    sim_both <- function(death_intercept) {
      simulate_true(
        n_pts_true        = n_chunk,
        seed              = seed,                   # same seed in both arms
        true_beta_arm     = p$true_beta_arm,
        Y_intercept       = p$Y_intercept,          # calibrated
        death_intercept   = death_intercept,        # -Inf = no death
        eff_frailty_y     = p$eff_frailty_y,
        eff_frailty_death  = p$eff_frailty_death,
        eff_frailty_train = p$eff_frailty_train,
        train_intercept   = p$train_intercept,      # calibrated
        train_frailty_eff = p$train_frailty_eff,
        max_frailty_red   = p$max_frailty_red       # depends on the effect (low / high)
      ) %>%
        mutate(arm = as.integer(ID > n_chunk)) %>%  # arm from ID (A1 is only stored at t = 0)
        dplyr::select(ID, t, arm, age, sex, L1, Yfull, Death)
    }
    
    # with death (scenario)
    dat <- sim_both(p$death_intercept)
    falls <- dat %>%                                # falls per arm and quarter
      filter(t >= 1) %>%
      mutate(q = cut(t, c(0, 13, 26, 39, 52), labels = paste0("Q", 1:4))) %>%
      group_by(arm, q) %>%
      summarise(y = sum(Yfull), .groups = "drop") %>%  
      tidyr::complete(arm, q, fill = list(y = 0))   # quarter without falls -> 0 (keeps Q1..Q4)
    deaths <- dat %>%                               # deaths per arm (death is absorbing)
      filter(t == max(t)) %>%
      group_by(arm) %>%
      summarise(d = sum(Death), .groups = "drop")
    ag <- prep_data_true(dat)  # counting process data for gl and lwyy
    rm(dat)
    
    # same models as in the analysis, fitted on the complete data of this block
    fit_lwyy <- coxph(Surv(start, stop, status == 1) ~ arm + age + sex + log_base_falls,
                      data = ag, ties = "breslow")
    fit_gl   <- recreg(Event(start, stop, status) ~ arm + age + sex + log_base_falls + cluster(ID),
                       data = ag, cause = 1, death.code = 3, cens.code = 0, control = list(tol = 1e-12))
    rm(ag)
    
    # without death (hypothetical)
    # mortality "none" (death_intercept = -Inf): the run above has no death already,
    # same seed + same parameters -> sim_both(-Inf) would give exactly the same data
    # -> reuse the falls of the run above instead of simulating again (saves time)
    if (is.infinite(p$death_intercept)) {
      falls_nd <- falls %>%
        group_by(arm) %>%
        summarise(y = sum(y), .groups = "drop")
    } else {
      dat_nd   <- sim_both(-Inf)
      falls_nd <- dat_nd %>%
        filter(t >= 1) %>%
        group_by(arm) %>%
        summarise(y = sum(Yfull), .groups = "drop")
      rm(dat_nd); gc()
    }
    
    y_q <- function(a) falls %>% filter(arm == a) %>% arrange(q) %>% pull(y)   # falls per quarter, Q1..Q4
    
    tibble(causal_id = p$causal_id,
           block     = block,
           y_ctrl    = sum(y_q(0)),
           y_int     = sum(y_q(1)),
           y_ctrl_nd = sum(falls_nd$y[falls_nd$arm == 0]),
           y_int_nd  = sum(falls_nd$y[falls_nd$arm == 1]),
           d_ctrl    = sum(deaths$d[deaths$arm == 0]),
           d_int     = sum(deaths$d[deaths$arm == 1]),
           !!!setNames(as.list(y_q(0)), paste0("y_ctrl_Q", 1:4)),
           !!!setNames(as.list(y_q(1)), paste0("y_int_Q", 1:4)),
           beta_lwyy = unname(coef(fit_lwyy)["arm"]),
           beta_gl   = unname(coef(fit_gl)["arm"]))
  }
  
  cat("Set up a parallel backend for true values\n")
  
  n_workers_true <- if (run_on_cluster) n_cores_use else min(n_cores_use, n_cores_true)
  
  cat("Workers for true values:", n_workers_true, "| tasks:", nrow(true_tasks), "\n")
  cl <- makeCluster(n_workers_true)   # each worker needs memory for one block (see n_chunk_true)
  registerDoParallel(cl)
  
  
  true_blocks <- foreach(j = seq_len(nrow(true_tasks)),
                         .combine = bind_rows,
                         .packages = c("simcausal", "dplyr", "tidyr", "survival", "mets"),
                         .export = c("COEF", "sim_dag", "simulate_true", "prep_data_true", "tie_offset", "true_block",
                                     "unique_causal_scenarios", "true_tasks", "n_chunk_true")) %dopar% {
                                       i <- true_tasks$causal_id[j]
                                       k <- true_tasks$block[j]
                                       # unique per scenario and block, separate from the data seeds (1e6 * causal_id + run, dropout + 5e8)
                                       # row by causal_id (not by position, the ids are not consecutive within a block)
                                       true_block(unique_causal_scenarios[unique_causal_scenarios$causal_id == i, ], k,
                                                  n_chunk_true, seed = 1e9 + 1e5 * i + k)
                                     }
  
  stopCluster(cl)
  
  # cluster: results_dir instead "./results" 
  saveRDS(true_blocks, file.path(results_dir, "00_true_blocks.rds"))    # block results: check MC error, add blocks later
  
  # Monte Carlo SE over the blocks (SE=s/sqrt(n))
  mcse <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ 
  
  true_estimands_df <- true_blocks %>%
    group_by(causal_id) %>%
    summarise(
      N_True = n() * n_chunk_true,                                   # patients per arm
      # Monte Carlo SE from the variation between the blocks (before the counts are summed up)
      MCSE_LogRateRatio         = mcse(log(y_int / y_ctrl)),
      MCSE_LogRateRatio_NoDeath = mcse(log(y_int_nd / y_ctrl_nd)),
      MCSE_Beta_LWYY            = mcse(beta_lwyy),
      MCSE_Beta_GL              = mcse(beta_gl),
      MCSE_LogRR_Q1             = mcse(log(y_int_Q1 / y_ctrl_Q1)),
      MCSE_LogRR_Q2             = mcse(log(y_int_Q2 / y_ctrl_Q2)),
      MCSE_LogRR_Q3             = mcse(log(y_int_Q3 / y_ctrl_Q3)),
      MCSE_LogRR_Q4             = mcse(log(y_int_Q4 / y_ctrl_Q4)),
      True_Beta_LWYY            = mean(beta_lwyy),                   # target of the adjusted LWYY
      True_Beta_GL              = mean(beta_gl),                     # target of the adjusted Ghosh-Lin
      across(c(starts_with("y_"), starts_with("d_")), sum),          # pooled counts over all blocks
      .groups = "drop"
    ) %>%
    mutate(
      MCF_Control               = y_ctrl    / N_True,                # mean number of falls until t.end per arm
      MCF_Intervention          = y_int     / N_True,
      MCF_Control_NoDeath       = y_ctrl_nd / N_True,
      MCF_Intervention_NoDeath  = y_int_nd  / N_True,
      True_LogRateRatio         = log(MCF_Intervention / MCF_Control),                 # ratio of mean falls (nonparametric GL, marginal)
      True_LogRateRatio_NoDeath = log(MCF_Intervention_NoDeath / MCF_Control_NoDeath), # hypothetical, no death
      LogRR_Q1                  = log(y_int_Q1 / y_ctrl_Q1),         # time-varying effect per quarter
      LogRR_Q2                  = log(y_int_Q2 / y_ctrl_Q2),
      LogRR_Q3                  = log(y_int_Q3 / y_ctrl_Q3),
      LogRR_Q4                  = log(y_int_Q4 / y_ctrl_Q4),
      Mort_Control              = d_ctrl / N_True,                   # achieved mortality per arm
      Mort_Intervention         = d_int  / N_True
    ) %>%
    left_join(unique_causal_scenarios %>% dplyr::select(causal_id, effect, all_of(causal_vars)),   
              by = "causal_id") %>%
    # same truth for all scenarios with the same causal parameters (e.g. 2, 3, 4 only differ in dropout)
    inner_join(base_scenarios %>%
                 distinct(across(all_of(causal_vars)), scenario_name),
               by = causal_vars, relationship = "one-to-many") %>%
    mutate(Scenario_Name = scenario_name, Mortality = mortality_level, True_Beta_Arm = true_beta_arm,
           Effect = effect, Max_Frailty_Red = max_frailty_red) %>%                                     
    
    # strict H0 (no direct effect and no effect via training): both arms have the same distribution
    # -> all effects are exactly 0, levels per arm are the same (mean of both arms = less MC error)
    # the MCSE columns still show the Monte Carlo noise of the simulated values
    mutate(
      h0 = True_Beta_Arm == 0 & train_frailty_eff == 0,
      across(c(True_LogRateRatio, True_LogRateRatio_NoDeath, True_Beta_LWYY, True_Beta_GL, starts_with("LogRR_Q")),
             ~ ifelse(h0, 0, .x)),
      MCF_Control              = ifelse(h0, (MCF_Control + MCF_Intervention) / 2, MCF_Control),
      MCF_Intervention         = ifelse(h0, MCF_Control, MCF_Intervention),
      MCF_Control_NoDeath      = ifelse(h0, (MCF_Control_NoDeath + MCF_Intervention_NoDeath) / 2, MCF_Control_NoDeath),
      MCF_Intervention_NoDeath = ifelse(h0, MCF_Control_NoDeath, MCF_Intervention_NoDeath),
      True_RateRatio           = exp(True_LogRateRatio)
    ) %>%
    dplyr::select(Scenario_Name, Mortality, True_Beta_Arm, Effect, Max_Frailty_Red,  
                  MCF_Control, MCF_Intervention, MCF_Control_NoDeath, MCF_Intervention_NoDeath,
                  True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath,
                  True_Beta_LWYY, True_Beta_GL, Mort_Control, Mort_Intervention,
                  LogRR_Q1, LogRR_Q2, LogRR_Q3, LogRR_Q4,
                  N_True, starts_with("MCSE_"))
  
  write_csv(true_estimands_df, file.path(results_dir, "00_true_estimands.csv"))

  # scenario list for 03
  write_csv(base_scenarios, file.path(results_dir, "00_base_scenarios.csv"))

  cat("True values are simulated and saved.\n")
  
}   # ---- end of "true" -------------------------------------------------




#------------------- Simulate all scenarios ----------------------------------------------##



if (run_data) {   # cluster ---- start "data" -------------------------------------
  
  # local all scenarios, on cluster only scenario of array_task
  # array task k = k-th scenario of the chosen blocks (scenario_id itself is not consecutive)
  if (run_on_cluster && !(array_id %in% seq_len(nrow(base_scenarios)))) {
    stop("Array-Task ", array_id, " don't exist: only 1 - ", nrow(base_scenarios), " scenarios in block(s) ",
         paste(scenario_blocks, collapse = ", "))
  }
  scenario_ids <- if (run_on_cluster) base_scenarios$scenario_id[array_id] else base_scenarios$scenario_id
  
  # Multiply by number of simulation runs
  scenarios <- expand_grid(
    base_scenarios %>% filter(scenario_id %in% scenario_ids),   # cluster
    sim_run = 1:n_sim
  ) %>%
    mutate(seed      = 1e6 * causal_id + sim_run,   # same complete data for all scenarios with the same causal parameters
           seed_drop = seed + 5e8,                  # own stream for dropout (same uniforms for all dropout mechanisms)
           file_path = file.path(results_dir,                    
                                 sprintf("Odat_scen%02d_run%04d.rds", scenario_id, sim_run)))
  
  scenarios_all <- scenarios   # complete list for the log (also data sets skipped after a restart)
  
  if (run_on_cluster) {
    # existing data sets are only skipped if they were generated with the same DGP code and parameters
    # (one settings file per scenario -> array tasks never write the same file)
    # deparse() ignores comments -> only real code changes count
    for (id in scenario_ids) {
      sig_data <- list(sim_dag           = deparse(sim_dag),
                       add_dropout       = deparse(add_dropout),          
                       simulate_scenario = deparse(simulate_scenario),
                       COEF              = COEF,
                       params            = base_scenarios %>% filter(scenario_id == id) %>% as.list())
      check_settings(sig_data, file.path(results_dir, sprintf("00_data_settings_scen%02d.rds", id)),
                     what = "Data sets")
    }

    done <- file.exists(scenarios$file_path)
    if (any(done)) cat(sum(done), "-", nrow(scenarios), "data already exist. Are skipped\n")
    scenarios <- scenarios[!done, ]
  } else {
    #local: remove all old data
    file.remove(list.files(results_dir, pattern = "^Odat_.*\\.rds$", full.names = TRUE))
  }
  
  if (nrow(scenarios) > 0) {
    
    # Setup parallel worker again for the main loop
    cat("Set up a parallel backend for main simulation\n")
    cl <- makeCluster(n_cores_use)   
    registerDoParallel(cl)
    
    cat("Start simulation with foreach:", nrow(scenarios), "data sets on", n_cores_use, "cores\n")
    t0 <- Sys.time()
    
    # Run simulations in parallel
    results_log <- foreach(i = 1:nrow(scenarios),
                           .combine = bind_rows,
                           .packages = c("simcausal", "dplyr", "tidyr", "tibble", "magrittr"),   # instead of tidyverse
                           .export = c("COEF", "sim_dag", "add_dropout", "simulate_scenario", "save_rds_atomic")) %dopar% {   

                             p <- scenarios[i, ]      # Current parameters

                             # Generate data
                             sim_data <- simulate_scenario(
                               seed = p$seed,
                               seed_drop = p$seed_drop,   
                               n_pts = p$n_pts,
                               true_beta_arm = p$true_beta_arm,
                               drop_beta_arm = p$drop_beta_arm,
                               Y_intercept = p$Y_intercept,
                               drop_intercept = p$drop_intercept,
                               death_intercept = p$death_intercept,
                               eff_frailty_y = p$eff_frailty_y,
                               eff_frailty_death = p$eff_frailty_death,
                               eff_frailty_drop = p$eff_frailty_drop,
                               eff_L1_drop = p$eff_L1_drop,
                               eff_fall_drop = p$eff_fall_drop,           
                               max_frailty_red = p$max_frailty_red,       
                               eff_frailty_train = p$eff_frailty_train,
                               train_intercept = p$train_intercept,
                               train_frailty_eff = p$train_frailty_eff
                             )
                             
                             # Save individual dataset (path set above: results_dir/Odat_scenXX_runYYYY.rds)
                             # atomic -> no broken file if the job is killed while saving
                             save_rds_atomic(sim_data, p$file_path)
                             
                             return(p)
                           }
    
    stopCluster(cl)
    cat("Data sets generated in", format(round(Sys.time() - t0, 1)), "\n")
    
  }   # end if (nrow(scenarios) > 0)
  
  # Save logbook 
  # local: one log for everything, cluster: one log per scenario (no simultaneous writing)
  results_log <- scenarios_all %>% filter(file.exists(file_path))
  missing <- nrow(scenarios_all) - nrow(results_log)
  if (missing > 0) warning(missing, " data sets are missing - restart the job")
  log_file <- if (run_on_cluster) {
    file.path(results_dir, sprintf("00_simulation_log_scen%02d.csv", scenario_ids))   
  } else {
    file.path(results_dir, "00_simulation_log.csv")
  }
  write_csv(results_log, log_file)
  cat(nrow(results_log), "data sets in log:", log_file, "\n")
  
}  