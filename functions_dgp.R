###############################################################################
### Functions for the data generating process (used in 01 and time_check.R)
### - Section 1: Fixed coefficients and DAG (complete data)
### - Section 2: Dropout and simulation of one trial
### - Section 3: Calibration of the intercepts
### - Section 4: True values (Monte Carlo)
###############################################################################

#------------------- Section 1: Fixed coefficients and DAG ---------------------------------##

COEF <- list(
  age_cent = 72.5,     # centered age (mid age)
  sex_cent = 0.5,      # centered sex 
  y_age = 0.05,        # effect age on falls (falls increase with age (RKI))
  y_sex = -0.1,        # effect sex on falls (woman fall more (RKI))
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

# rng: state of the random number stream (rng_stream(), functions.R), complete data = substream 0
sim_dag <- function(n, rng, arm_value = NA,
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
   

  DAG <- set.DAG(D)                    # lock first: set.DAG may draw random numbers for its own test simulation
  use_stream(rng)                      # start of the stream of this data set (functions.R)
  sim(DAG = DAG, n = n, wide = FALSE)  # no rndseed: simcausal must not reset the seed itself
}


#------------------- Section 2: Dropout and simulation of one trial ------------------------##

# dropout afterwards, from its own random stream (substream 1 of the data stream, rng_drop)
# same uniforms U[person, week] for all dropout mechanisms (only depend on rng_drop)
# scenarios with the same causal parameters have identical complete data and coupled dropout
# Dropout[t] = 1 if Dropout[t-1] = 1; 0 if death at the end of week t (or before), else Bernoulli(p_t)
# Yobs[t]    = Yfull[t] * (1 - Dropout[t-1])  (falls of the dropout week are observed)
# Scenario 4: eff_fall_drop: dropout depends on log(1 + falls up to t-1) (= log_falls_lag in the weight models). While a
# person is in the study, the observed falls are the complete falls -> computed from Yfull
add_dropout <- function(dat, rng_drop, t.end = 52,
                        drop_intercept, drop_beta_arm = 0,
                        eff_frailty_drop = 0, eff_fall_drop = 0) {
  use_stream(rng_drop)
  ids <- sort(unique(dat$ID))
  U   <- matrix(runif(length(ids) * t.end), nrow = length(ids))   # U[person, week]

  dat %>%
    arrange(ID, t) %>%
    group_by(ID) %>%
    mutate(arm_i = A1[1], age_i = age[1], sex_i = sex[1],                 # baseline values (t = 0)
           log_falls_lag = log1p(cumsum(Yfull) - Yfull)) %>%            # log(1 + falls up to t-1) cumsum(Yfull) = all falls up to week t, Yfull = falls in week t
    # log1p computes log(1+x) also für x<<1
    ungroup() %>%
    mutate(
      p_drop   = plogis(drop_intercept +
                          COEF$cens_age * (age_i - COEF$age_cent) + COEF$cens_sex * (sex_i - COEF$sex_cent) +
                          eff_frailty_drop * FrailtyT +
                          eff_fall_drop * log_falls_lag +                             
                          drop_beta_arm * arm_i),
      u        = ifelse(t >= 1, U[cbind(match(ID, ids), pmax(t, 1))], 1),  # random number for each person i and week t
      drop_now = t >= 1 & Death == 0 & u < p_drop          # dropout at the end of week t, not in the death week, if random number smaller than probability. U unform distributed -> P(U<p)=p
    ) %>%
    group_by(ID) %>%
    mutate(Dropout = as.integer(cumsum(drop_now) > 0),      # absorbing, only one dropout possible
           Yobs    = Yfull * (1 - lag(Dropout, default = 0L))) %>%
    ungroup() %>%
    dplyr::select(-arm_i, -age_i, -sex_i, -log_falls_lag, -p_drop, -u, -drop_now)
}


# Function to simulate data
simulate_scenario <- function(
    rng,                             # stream of the complete data (rng_stream(data_stream_id(causal_id, sim_run)),
                                     #   same for all scenarios with the same causal parameters)
    rng_drop = sub_stream(rng, 1),   # dropout uniforms: substream 1 of the same stream
    n_pts = 1800,                    # number of patients (~900 per arm)
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
    eff_fall_drop = 0.0       # Scenario 4:Effect observed fall rate up to t-1 on dropout (positive control, should be corrected by ipcw)

) {
  
  # complete data (DAG without dropout, defined in sim_dag())
  Odat <- sim_dag(n = n_pts, rng = rng, t.end = t.end,
                  true_beta_arm = true_beta_arm, frailty_sd = frailty_sd,
                  Y_intercept = Y_intercept, death_intercept = death_intercept,
                  train_intercept = train_intercept,
                  train_frailty_eff = train_frailty_eff, max_frailty_red = max_frailty_red,
                  eff_frailty_y = eff_frailty_y, eff_frailty_death = eff_frailty_death,
                  eff_frailty_train = eff_frailty_train)

  # dropout 
  Odat <- add_dropout(Odat, rng_drop = rng_drop, t.end = t.end,
                      drop_intercept = drop_intercept, drop_beta_arm = drop_beta_arm,
                      eff_frailty_drop = eff_frailty_drop,
                      eff_fall_drop = eff_fall_drop)   

  # check for NA
  if (anyNA(Odat$Dropout) | anyNA(Odat$Yobs)) {
    stop("NA in Dropout oder Yobs: check parameter")
  }
  
  Odat %<>%  # save mutate changes
    arrange(ID, t) %>%  # arrange for ID and t
    group_by(ID) %>%
    mutate(
      arm = A1[1],
      age = age[1],
      sex = sex[1],
      base_risk = L1[1],             # proxy 1: falls in the year before the study
      frailty = Frailty[1],          # baseline frailty
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



# Summary of one simulated data set, one row with the columns ctrl_* and int_* (per arm)
# saved in the simulation log -> 04_overview.R describes all runs without the data sets
# shares as proportions; fq1 ... fq4: quartiles of the baseline frailty in this data set
summarise_odat <- function(Odat) {
  per <- Odat %>%
    group_by(ID) %>%
    summarise(arm = first(arm), age = first(age), sex = first(sex), base_risk = first(base_risk),
              frailty = first(frailty), dropout_time = first(dropout_time), death_time = first(death_time),
              falls_full = sum(n_falls_full), falls_obs = sum(n_falls, na.rm = TRUE),
              dropped = first(ever_dropped), died = first(ever_died),
              train_share = mean(Train), train_w1 = first(Train), .groups = "drop") %>%
    mutate(fq = ntile(frailty, 4))

  stats_arm <- per %>%
    group_by(arm) %>%
    summarise(n = n(), age_mean = mean(age), male = mean(sex),
              L1_0 = mean(base_risk == 0), L1_1 = mean(base_risk == 1), L1_2plus = mean(base_risk >= 2),
              fall1 = mean(falls_full >= 1), fall2 = mean(falls_full >= 2), fall5 = mean(falls_full >= 5),
              falls_mean = mean(falls_full), falls_obs_mean = mean(falls_obs), falls_max = max(falls_full),
              top10 = sum(sort(falls_full, decreasing = TRUE)[seq_len(ceiling(0.1 * n()))]) / sum(falls_full),   # share of all falls of the top 10 % persons
              dropout = mean(dropped), death = mean(died), train_w1 = mean(train_w1), train_share = mean(train_share),
              drop_w13 = mean(!is.na(dropout_time) & dropout_time <= 13), drop_w26 = mean(!is.na(dropout_time) & dropout_time <= 26),
              drop_w39 = mean(!is.na(dropout_time) & dropout_time <= 39), drop_w52 = mean(!is.na(dropout_time) & dropout_time <= 52),
              death_w13 = mean(!is.na(death_time) & death_time <= 13), death_w26 = mean(!is.na(death_time) & death_time <= 26),
              death_w39 = mean(!is.na(death_time) & death_time <= 39), death_w52 = mean(!is.na(death_time) & death_time <= 52),
              .groups = "drop")

  stats_fq <- per %>%
    group_by(arm, fq) %>%
    summarise(fall1 = mean(falls_full >= 1), falls_mean = mean(falls_full), dropout = mean(dropped),
              death = mean(died), train_share = mean(train_share), .groups = "drop") %>%
    pivot_wider(names_from = fq, values_from = c(fall1, falls_mean, dropout, death, train_share),
                names_glue = "{.value}_fq{fq}")

  stats_arm %>%
    left_join(stats_fq, by = "arm") %>%
    pivot_longer(-arm, names_to = "measure", values_to = "value") %>%
    mutate(measure = paste0(ifelse(arm == 0, "ctrl_", "int_"), measure)) %>%
    dplyr::select(measure, value) %>%
    pivot_wider(names_from = measure, values_from = value)
}


#------------------- Section 3: Calibration of intercepts ----------------------------------##
# The intercepts are chosen so that the rates in the whole population (with frailty) hit the targets
# coefficients from COEF -> automatically the same as in sim_dag()
# cd: reference population with the same distribution as in the data generation (see calibrate_scenarios())

# Y_intercept so that number of units with >= 1 fall in 52 weeks = target (for control arm, no death)
# P(N >= 1) = 1 - exp(-52 * weekly rate) (Poisson), averaged over all persons
calib_Y <- function(cd, eff_frailty_y, target = 0.30, t.end = 52) {
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

# dropout intercept (scenarios without eff_fall_drop): death as competing event, a dropout at the end of
# week t is only possible if the person is alive at the end of week t (add_dropout())
# P(dropout until week 52) = sum_t [(1-p)(1-q)]^(t-1) (1-q) p   (p, q: weekly dropout / death probability)
calib_drop <- function(lp, target, q = 0, t.end = 52) {
  if (all(q == 0)) return(calib_logit(lp, target, t.end))       # without death
  if (target == 0) return(-Inf)
  uniroot(function(b) {
    p <- plogis(b + lp)
    r <- (1 - p) * (1 - q)                                        # alive and still in the study after one week
    mean((1 - q) * p * (1 - r^t.end) / (1 - r)) - target
  }, c(-20, 5))$root
}

# training: calibrated only for week 1, because later weeks depend on the own training history
# (training -> frailty -> training). In week 1 CumTrain = 0 -> FrailtyT = Frailty, so
# calib_logit with t.end = 1 gives E[plogis(b + lp)] = target

# dropout depending on log(1 + falls up to t-1) (scenario 4): the weekly probability changes over time,
# so there is no closed formula. The weekly falls of the reference population are simulated
# (control arm: no training -> constant rate, no death as in the other calibrations), then
# P(no dropout until week 52) = prod_t (1 - plogis(b + lp + eff_fall_drop * log_falls_lag_t)), averaged
calib_drop_fall <- function(cd, lp, eff_fall_drop, Y_intercept, eff_frailty_y, target, q= 0,  t.end = 52) {
  lambda <- exp(Y_intercept + COEF$y_age * (cd$age - COEF$age_cent) + COEF$y_sex * (cd$sex - COEF$sex_cent) +
                  eff_frailty_y * cd$Z)                 # weekly fall rate
  set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")   # fixed -> same intercept in every array task
  # generator given explicitly: after rng_stream() the session uses L'Ecuyer-CMRG, the intercepts must not change
  Y <- matrix(rpois(nrow(cd) * t.end, rep(lambda, t.end)), nrow = nrow(cd))   # Y[person, week]
  cum_lag <- matrix(0, nrow(cd), t.end)                                     # falls up to t-1
  for (k in 2:t.end) cum_lag[, k] <- cum_lag[, k - 1] + Y[, k - 1]
  log_falls_lag <- log1p(cum_lag)                                           # = log_falls_lag in the analysis
  uniroot(function(b) {
    p <- plogis(b + lp + eff_fall_drop * log_falls_lag)                     # weekly dropout probability [person, week]
    if (all(q == 0)) return(mean(1 - exp(rowSums(log1p(-p)))) - target)    # without death
    S   <- rep(1, nrow(cd))                                                 # alive and in the study at the start of week k
    tot <- rep(0, nrow(cd))                                                 # P(dropout until week k)
    for (k in seq_len(t.end)) {
      tot <- tot + S * (1 - q) * p[, k]
      S   <- S * (1 - q) * (1 - p[, k])
    }
    mean(tot) - target
  }, c(-20, 5))$root
}

# intercepts for all scenarios (one row per scenario), targets from GlobalVars.R
calibrate_scenarios <- function(scen) {
  # reference population with the same distribution as in the data generation (frailty_sd = 1)
  # generator given explicitly: after rng_stream() the session uses L'Ecuyer-CMRG, the intercepts must not change
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  cd <- tibble(age = runif(2e5, 65, 80), sex = rbinom(2e5, 1, 0.5), Z = rnorm(2e5))  # Z = frailty

  scen %>%
    rowwise() %>%                          # calibration for each scenario separately
    mutate(
      Y_intercept     = calib_Y(cd, eff_frailty_y, target_falls),
      death_intercept = calib_logit(COEF$d_age * (cd$age - COEF$age_cent) + COEF$d_sex * (cd$sex - COEF$sex_cent) +
                                      eff_frailty_death * cd$Z, target_death[[mortality_level]]),
      drop_intercept  = {
        lp_drop <- COEF$cens_age * (cd$age - COEF$age_cent) + COEF$cens_sex * (cd$sex - COEF$sex_cent) +
          eff_frailty_drop * cd$Z
        # weekly death probability in the control arm (frailty constant without training), 0 without death
        q_death <- plogis(death_intercept + COEF$d_age * (cd$age - COEF$age_cent) +
                            COEF$d_sex * (cd$sex - COEF$sex_cent) + eff_frailty_death * cd$Z)
        if (eff_fall_drop == 0) calib_drop(lp_drop, target_dropout, q_death)                                       # scenario 1, 2, 3, 5
        else calib_drop_fall(cd, lp_drop, eff_fall_drop, Y_intercept, eff_frailty_y, target_dropout, q = q_death)  # scenario 4
      },
      train_intercept = calib_logit(COEF$train_age * (cd$age - COEF$age_cent) + COEF$train_sex * (cd$sex - COEF$sex_cent) +
                                      eff_frailty_train * cd$Z, target_train, t.end = 1)
    ) %>%
    ungroup()
}


#------------------- Section 4: True values (Monte Carlo) ----------------------------------##


# simulate a very large trial per scenario without dropout (complete data, only death stops falls). Then
#   - mean number of falls per arm  -> True_LogRateRatio
#   - fit the same models as in the analysis (LWYY, Ghosh-Lin, adjusted for age, sex,
#     log_base_falls on complete data -> True_Beta_LWYY, True_Beta_GL)
# The true value of a method is what it estimates with complete data and very large n.
# Both arms start at the same state of the stream (same frailty, age, sex) -> less Monte Carlo noise.

simulate_true <- function(n_pts_true = 50000,
                          rng,                   # stream state (rng_stream(true_stream_id(causal_id, block)))
                          ...) {                 # all other DGP parameters -> sim_dag()
  bind_rows(
    sim_dag(n = n_pts_true, rng = rng, arm_value = 0, ...),                                # control
    sim_dag(n = n_pts_true, rng = rng, arm_value = 1, ...) %>%                             # training
      mutate(ID = ID + n_pts_true)  
  )
}

# Prepare counting-process data for true model
# complete data -> one row per interval between falls 
# k falls in week t are placed at t-1 + j/(k+d), j = 1..k (d = 1 if death at the end of week t) minus tie_offset(ID) (no ties between persons)
prep_data_true <- function(dat, t.end = 52) {
  base <- dat %>%
    filter(t == 0) %>%
    transmute(ID, arm, age, sex, log_base_falls = log1p(L1))   # log(1 + falls in the year before)
  
  wk <- dat %>%                                         # all study weeks without baseline week (t=0)
    filter(t >= 1) %>%
    group_by(ID) %>%
    mutate(death_time = suppressWarnings(min(t[Death == 1]))) %>%   # Inf = no death
    ungroup() %>%
    filter(t <= death_time)                                           # falls of the death week count
  
  falls <- wk %>%
    filter(Yfull > 0) %>%            # only for weeks with minimum 1 fall
    mutate(d = as.integer(t == death_time)) %>%   # d=1, if death at end of week
    tidyr::uncount(Yfull, .id = "j", .remove = FALSE) %>%           # k falls -> k rows
    transmute(ID, stop = t - 1 + j / (Yfull + d) - tie_offset(ID), status = 1L)   # no ties between persons
  
  ends <- wk %>%
    group_by(ID) %>%
    summarise(stop   = min(first(death_time), t.end) + 0.001,   # shortly after falls at the same time (as in prep_counting)
              status = ifelse(is.finite(first(death_time)), 3L, 0L), .groups = "drop") 
  
  bind_rows(falls, ends) %>%
    arrange(ID, stop, status != 1L) %>%          # at equal times: fall before end of study
    group_by(ID) %>%
    mutate(start = lag(stop, default = 0)) %>%
    ungroup() %>%
    filter(stop > start) %>%                     # drops a zero-length end row (fall in week 52)
    left_join(base, by = "ID") %>%
    as.data.frame()
}

# one block of the true value simulation for scenario p
# both arms with the same stream state, once with and once without death
true_block <- function(p, block, n_chunk, rng) {
  
  # only the death intercept differs between the calls, only the columns needed later are kept
  sim_both <- function(death_intercept) {
    simulate_true(
      n_pts_true        = n_chunk,
      rng               = rng,                    # same stream state in both arms
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
    tidyr::complete(arm, q, fill = list(y = 0))   # quarter without falls -> 0 
  deaths <- dat %>%                               # deaths per arm (death is absorbing)
    filter(t == max(t)) %>%
    group_by(arm) %>%
    summarise(d = sum(Death), .groups = "drop")
  ag <- prep_data_true(dat)  # counting process data for gl and lwyy
  rm(dat)
  
  # same models as in the analysis, fitted on the complete data of this block
  fit_lwyy <- coxph(Surv(start, stop, status == 1) ~ arm + age + sex + log_base_falls,
                    data = ag, ties = "breslow")
  # Armijo line search as in fit_models() (gl_ctrl): without it Newton-Raphson diverged (NaN) in a block with one
  # very frail person (144 falls, supp causal_id 12 block 158); normal blocks unchanged (difference < 1e-7)
  fit_gl   <- recreg(Event(start, stop, status) ~ arm + age + sex + log_base_falls + cluster(ID),
                     data = ag, cause = 1, death.code = 3, cens.code = 0,
                     control = list(tol = 1e-12, backtrack = "armijo"))

  # model-based mean falls per arm of the LWYY fit (model_mcf, functions.R) = limit of the LWYY falls
  # on complete data, the same kind of true value as True_Beta_LWYY (death ends the risk set)
  # standardised over the persons of the block (both arms have the same baseline -> arm 0 is enough)
  mcf_lwyy <- model_mcf(fit_lwyy, ag %>% filter(arm == 0) %>% distinct(ID, age, sex, log_base_falls))
  rm(ag)
  
  # without death (hypothetical)
  # mortality "none" (death_intercept = -Inf): the run above has no death already,
  # same stream + same parameters -> sim_both(-Inf) would give exactly the same data
  # -> reuse the falls of the run above instead of simulating again 
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
         beta_gl   = unname(coef(fit_gl)["arm"]),
         mcf_lwyy_ctrl = mcf_lwyy[1],
         mcf_lwyy_int  = mcf_lwyy[2])
}
