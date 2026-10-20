###############################################################################
### Function ##################################################################
### - Section 1: Data preparation
### - Section 2: Function for analyse
### - Section 3: weights and MCF
### - Section 4: Bootstrap
###############################################################################

#------------- Section 1: Data preparation -----------------------------------

# Ties: falls of different persons at exavtly the same time
# Different handled by coxph and mets 
# Every fall time of a person is moved by the same small person-specific offset (0, eps) to the left
# - deterministic, not related to the arm (golden ratio hash of ID)
# - same offset for all falls of a person -> order within a person stays the same
# - risk set is constant within a week (death/censoring only at t + 0.001) -> Breslow estimate unchanged
tie_offset <- function(ID, eps = 0.01) eps * ((ID * 0.6180339887498949) %% 1)


# Prepare counting-process data for true model
# complete data -> one row per interval between falls 
# k falls in week t are placed at t-1 + j/(k+d), j = 1..k (d = 1 if death at the end of week t) minus tie_offset(ID) (no ties between persons)
prep_data_true <- function(dat, t.end = 52) {
  base <- dat %>%
    filter(t == 0) %>%
    transmute(ID, arm, age, sex, log_base_falls = log1p(L1))   # log(1 + falls in the year before)
  
  wk <- dat %>%
    filter(t >= 1) %>%
    group_by(ID) %>%
    mutate(death_time = suppressWarnings(min(t[Death == 1]))) %>%   # Inf = no death
    ungroup() %>%
    filter(t <= death_time)                                           # falls of the death week count
  
  falls <- wk %>%
    filter(Yfull > 0) %>%
    mutate(d = as.integer(t == death_time)) %>%
    tidyr::uncount(Yfull, .id = "j", .remove = FALSE) %>%
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



# Counting-process data for LWYY and Ghosh-Lin (one data set for both models)
# Week with k falls is split into k sub-intervals (+1 if death/dropout at the end of the week (for recreg))
# status (coded for recreg; LWYY uses status == 1, all other codes are no event):
# 1 = fall, 3 = death, 0 = censoring (dropout or end of study, only in the last row of a person),
# 9 = interval without event (recreg counts every row with status 0 as censoring)
# Death and censoring lie 0.001 after the end of the week (falls at the same time count first)
# Last row is a fall (fall at the very end of week 52): a short censoring row is added.
# The column t is kept to merge the IPCW weight of week t.
# data = "observed": only observed weeks (n_falls, drop, death)
# data = "full"    : complete data for the oracle M0 (n_falls_full, death_full)
prep_counting <- function(dat, data = c("observed", "full")) {
  data <- match.arg(data)
  d <- if (data == "observed") {
    dat %>% filter(observed) %>% mutate(k = n_falls, term_death = death, term_drop = drop)
  } else {
    dat %>% mutate(k = n_falls_full, term_death = death_full, term_drop = 0L)
  }
  cp <- d %>%
    arrange(ID, t) %>%
    mutate(m = pmax(1L, as.integer(k + term_death + term_drop))) %>%  # number of falls and death or dropout event if nessesary
    tidyr::uncount(m, .id = "j", .remove = FALSE) %>%  # duplicate rows for every event
    mutate(start = t - 1 + (j - 1) / m,                # divide in parts
           stop  = t - 1 + j / m) %>%
    group_by(ID) %>%
    mutate(last_row = row_number() == n()) %>%
    ungroup() %>%
    mutate(status = case_when(j <= k                   ~ 1L,   # fall
                              j == m & term_death == 1 ~ 3L,   # death at the end of the week
                              last_row                 ~ 0L,   # dropout or end of study = censoring
                              TRUE                     ~ 9L))  # no event
  
  # last row is a fall: short censoring row at the end of study (administrative censoring)
  extra <- cp %>%
    filter(last_row & status == 1L) %>%
    mutate(start = stop, status = 0L)
  
  bind_rows(cp, extra) %>%
    mutate(stop = ifelse(status %in% c(0L, 3L), stop + 0.001, stop),        # ties: falls first
           stop = ifelse(status == 1L, stop - tie_offset(ID), stop)) %>%    # no ties between persons
    arrange(ID, start, stop) %>%
    group_by(ID) %>%
       mutate(start = lag(stop, default = 0)) %>%                              # next row starts at the moved fall
       ungroup() %>%
    dplyr::select(-k, -m, -j, -term_death, -term_drop, -last_row) %>%
    as.data.frame()
}


#--------------- Section 2: Function for analyse --------------------------------



# extract result out of cox model (LWYY) or recreg model (Ghosh-Lin)
# true_beta: the true value that belongs to this method (True_Beta_LWYY or True_Beta_GL)
# no rounding here -> round only when reporting (Monte Carlo errors are small)
extract_model_results <- function(model, model_name, true_beta) {
  
  est <- unname(coef(model)["arm"])
  se  <- sqrt(vcov(model)["arm", "arm"])         # robust (cluster) SE
  
  tibble(
    Model      = model_name,
    True_Beta  = true_beta,
    Estimate   = est,                            # log rate ratio (LWYY) / log mean ratio (GL)
    Rate_Ratio = exp(est),
    SE         = se,
    CI_low     = est - qnorm(0.975) * se,
    CI_high    = est + qnorm(0.975) * se,
    P_Value    = 2 * pnorm(-abs(est / se)),      # Wald test
    Reject     = P_Value < 0.05,                 # power (H1) / type I error (H0)
    Bias       = est - true_beta,
    Covered    = CI_low <= true_beta & true_beta <= CI_high
  )
}

# fit models
# Odat_all: complete dataset without dropout
# Odat_w: obseved data set with ipcw weights (with compute weights)
# same argument order as mcf_estimate(): fit_models(Odat_all, Odat_w)
# fit LWYY and Gosh-Lin model with and without proxies

fit_models <- function(Odat_all, Odat_w){
    cp      <- prep_counting(Odat_w, "observed")
    cp_full <- prep_counting(Odat_all, "full")
    
    f_lwyy <- Surv(start, stop, status == 1) ~ arm + age + sex + log_base_falls + cluster(ID)
    f_gl   <- Event(start, stop, status) ~ arm + age + sex + log_base_falls + cluster(ID)
    
    # GL: cens.model = ~ strata(arm) -> internal censoring Kaplan-Meier per arm
    list(
      "LWYY 0: Full data (oracle M0)" = coxph(f_lwyy, data = cp_full, ties = "breslow"),
      "LWYY 1: Naive (Unweighted)"    = coxph(f_lwyy, data = cp, ties = "breslow"),
      "LWYY 3: IPCW sym"              = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_sym),
      "LWYY 4: IPCW sym_train"        = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_symtrain),
      "LWYY 5: IPCW oracle"           = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_oracle),
      "GL 0: Full data (oracle M0)"   = recreg(f_gl, data = cp_full, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), control = list(tol = 1e-12)),
      "GL 1: Naive (Unweighted)"      = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), control = list(tol = 1e-12)),
      "GL 3: IPCW sym"                = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_sym, control = list(tol = 1e-12)),
      "GL 4: IPCW sym_train"          = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_symtrain, control = list(tol = 1e-12)),
      "GL 5: IPCW oracle"             = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_oracle, control = list(tol = 1e-12))  # stricter convergence (default tol 1e-6 stops too early)
    )
}

# ----------------- Section 3: Weights and mean number of falls -------------

# Compute ipw weights
# IPCW weights per arm (for 02_analyse_weekly.R. 03_mehrfach analyse and bootstrap)
# iptwm fitted on weeks without death, lag_ipw shifts weigths by one week (falls of week t are observed before a dropout at its end)
# w_* stabilised weights, wm_* for marginal mean (only time in numerator)
# dat: observed weekly data (observed == TRUE) with tstop, log_base_falls, fall_rate_lag, train_rate_lag, frailty_t


compute_weights <- function(dat){
  dat %>%
    mutate(tstart = t - 1) %>%    
    split(.$arm) %>%
    lapply(function(d_all){
      d_all <- d_all %>% arrange(ID, t) %>% as.data.frame()
      d_fit <- d_all %>% filter(death == 0)          # weeks in which a dropout is possible
      shift <- function(fit)lag_ipw(fit, d_fit, d_all) # shift ipcw by one week
      
      w_sym <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
        numerator = ~ age + sex + log_base_falls + ns(tstop, df = 3), 
        denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + fall_rate_lag,   
        id = ID, tstart = tstart,  # is ignored with family "binomial" 
        timevar = tstop,
        type = "cens",
        data = d_fit))
      
      w_oracle <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                           numerator = ~ age + sex + log_base_falls + ns(tstop, df = 3), 
                           denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + frailty_t,   
                           id = ID, tstart = tstart,  # is ignored with family "binomial" 
                           timevar = tstop,
                           type = "cens",
                           data = d_fit))
      # unstabilisied weights
      wm_sym <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                     numerator = ~ ns(tstop, df = 3), 
                     denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + fall_rate_lag,   
                     id = ID, tstart = tstart,  # is ignored with family "binomial" 
                     timevar = tstop,
                     type = "cens",
                     data = d_fit))
      wm_oracle <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                        numerator = ~ ns(tstop, df = 3), 
                        denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + frailty_t,   
                        id = ID, tstart = tstart,  # is ignored with family "binomial" 
                        timevar = tstop,
                        type = "cens",
                        data = d_fit))
      
      
      if(d_all$arm[1] == 1){
        w_symtrain <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                                  numerator   = ~ age + sex + log_base_falls + ns(tstop, df = 3),
                                  denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + fall_rate_lag + train_rate_lag,
                                  id = ID, tstart = tstart, timevar = tstop, type = "cens", data = d_fit))
        wm_symtrain <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                                   numerator   = ~ ns(tstop, df = 3),
                                   denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + fall_rate_lag + train_rate_lag,
                                   id = ID, tstart = tstart, timevar = tstop, type = "cens", data = d_fit))
      } else {                                       # control arm never trains -> same as sym
        w_symtrain  <- w_sym
        wm_symtrain <- wm_sym
      }
      d_all %>% mutate(w_sym = w_sym, w_symtrain = w_symtrain, w_oracle = w_oracle,
                       wm_sym = wm_sym, wm_symtrain = wm_symtrain, wm_oracle = wm_oracle)
    }) %>%
    bind_rows() %>%
    arrange(ID, t)
    }
  


# Weights from ipwtm() shifted by one week
#   ipwtm() returns cumulative weights including the current week (probability of no dropout
#   at the end of week t). The falls of week t happen before a dropout at the end of week t,
#   so the weight for week t is the ipwtm weight of week t-1 (week 1: weight 1).
#   ipwtm() is fitted on weeks without death (d_fit); a death week gets the weight of the
#   week before, which is exactly what the shift does.
#   ipw_fit: result of ipwtm(..., data = d_fit)
#   d_fit  : data used in ipwtm (weeks without death, sorted by ID, t)
#   d_all  : all observed weeks of the same arm (sorted by ID, t) -> weights aligned with d_all
lag_ipw <- function(ipw_fit, d_fit, d_all) {
  d_fit %>%
    dplyr::select(ID, t) %>%
    mutate(w_cum = ipw_fit$ipw.weights) %>%
    right_join(d_all %>% dplyr::select(ID, t), by = c("ID", "t")) %>%
    arrange(ID, t) %>%
    group_by(ID) %>%
    mutate(w = lag(w_cum, default = 1)) %>%   # shift by one entry
    ungroup() %>%
    pull(w)
}


# Mean number of falls up to t.end per arm (nonparametric Ghosh-Lin / Cook-Lawless, weighted)
#   mu_a = sum_t S_a(t-1) * dR_a(t)
#   dR_a(t): weighted mean number of falls in week t among persons at risk (Nelson-Aalen)
#   S_a(t) : weighted Kaplan-Meier for death (death at the end of week t)
#   Use UNSTABILISED weights (or weights whose numerator depends only on arm and time):
#   a numerator with covariates would change the covariate mix of the pseudo-population.
# d: weekly data, one row per person and week at risk, columns arm, t, falls, death, weight
weighted_mcf <- function(d, w = "w", falls = "n_falls", death = "death", t.end = 52) {
  d %>%
    mutate(.w = .data[[w]], .y = .data[[falls]], .dth = .data[[death]]) %>%
    group_by(arm, t) %>%
    summarise(dR = sum(.w * .y) / sum(.w),          # fallsrate weighted
              hD = sum(.w * .dth) / sum(.w),        # death hazard weighted
              .groups = "drop") %>%
    group_by(arm) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(S_lag = c(1, head(cumprod(1 - hD), -1))) %>%   # alive at start of week t
    summarise(MCF = sum((S_lag * dR)[t <= t.end]),         # mean falls in the real world (Ghosh-Lin)
              R   = sum(dR[t <= t.end]),                   # cumulative rate among survivors
              .groups = "drop")
}


# Mean number of falls up to quarter of year per arm (nonparametric Ghosh-Lin / Cook-Lawless, weighted)

weighted_mcf_q <- function(d, w = "w", falls = "n_falls", death = "death", breaks =c(0,13,26,39,52)) {
  d %>%
    mutate(.w = .data[[w]], .y = .data[[falls]], .dth = .data[[death]]) %>%
    group_by(arm, t) %>%
    summarise(dR = sum(.w * .y) / sum(.w),          # fallsrate weighted
              hD = sum(.w * .dth) / sum(.w),        # death hazard weighted
              .groups = "drop") %>%
    group_by(arm) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(S_lag = c(1, head(cumprod(1 - hD), -1)),        # alive at start of week t
           q     = cut(t, breaks, labels = paste0("Q", seq_len(length(breaks) - 1)))) %>%   # cut at time points
    filter(!is.na(q)) %>%
    group_by(arm, q) %>%
    summarise(MCF = sum(S_lag * dR), .groups = "drop")
}





# Model-based mean number of falls up to t.end per arm (standardisation)
# prediction for every person with arm set to 0 / 1, averaged over all persons
# base_cov: one row per person with the covariates of the model (age, sex, log_base_falls)
# Ghosh-Lin (recreg): mu_0(t) * exp(X beta)     -> real world with death (compare with MCF_*)
# LWYY (coxph):       Lambda_0(t) * exp(X beta) -> hypothetical world without death (compare with MCF_*_NoDeath),
# only valid if death is independent of the falls given the covariates
# log(MCF_1 / MCF_0) is exactly the coefficient of arm; new is only the level per arm
# basehaz with centered or not???
model_mcf <- function(fit, base_cov, t.end = 52) {
  if (inherits(fit, "coxph")) {                              # LWYY
    H0   <- basehaz(fit, centered = TRUE)                    # Lambda_0(t) at the centering values zbar = fit$means, Breslow estimator
    H0_t <- max(H0$hazard[H0$time <= t.end])                 # (centered = FALSE can over-/underflow, see ?basehaz)
    b    <- coef(fit)
    zbar <- fit$means[names(b)]
    sapply(0:1, function(a) {
      X <- as.matrix(base_cov %>% mutate(arm = a) %>% dplyr::select(all_of(names(b))))
      H0_t * mean(exp(sweep(X, 2, zbar) %*% b))              # Lambda_0(tau; zbar) * mean(exp((X_i - zbar) b))
    })
  } else {                                                   # Ghosh-Lin
    sapply(0:1, function(a)
      mean(predict(fit, newdata = base_cov %>% mutate(arm = a), times = t.end, se = FALSE)$cumhaz))
  }
}

# fits data prepared by fit_models()
mcf_estimate <- function(Odat_all, Odat_w, fits){
  # non-parametric estimate
  mcf_list <- list(
    "MCF 0: Full data (oracle M0)" = weighted_mcf(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full"), # w1 = 1: unweighted, Estimate without dropout
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
}


#----------- Section 4: Bootstrap -------------------------------------

# Draw bootstrap sample
boot_sample <- function(dat){
  draw <- dat %>%
    distinct(ID, arm) %>%  # one line per patient
    group_by(arm) %>%          # seperate for each arm
    slice_sample(prop = 1, replace = TRUE) %>%  # draw with replacemant same number of patient per arm
    ungroup() %>%
    mutate(ID_new = row_number()) %>%           # new ID, if one unit is sampled twice
    dplyr::select(ID, ID_new)
  
  draw <- draw %>%
    inner_join(dat, by ="ID", relationship = "many-to-many") %>%  # join all weeks
    mutate(ID=ID_new) %>%
    dplyr::select(-ID_new) %>%
    arrange(ID,t)
}

# Estimate weights and models
boot_estimate <- function(Odat_all){
  Odat_w  <- compute_weights(Odat_all %>% filter(observed)) 
  fits    <- fit_models(Odat_all, Odat_w)                      
  mcf_est <- mcf_estimate(Odat_all, Odat_w, fits)             
  c(sapply(fits,    function(f)  unname(coef(f)["arm"])),      
    sapply(mcf_est, function(mu) log(mu[2] / mu[1])))
  
}


# Run bootstrap

run_bootstrap <- function(Odat_all, n_boot){
  est <- lapply(seq_len(n_boot), function(b)
    tryCatch(boot_estimate(boot_sample(Odat_all)), error = function(e) NULL))
  est <- est[!vapply(est, is.null, logical(1))]                                   # save in matrix
  if (length(est) == 0) stop("all Bootstrap-samples failed")
  est <- do.call(rbind, est) 
  
  tibble(
    Model        = colnames(est),
    Boot_SE      = apply(est, 2, sd, na.rm = TRUE),                               # Bootstrap SE
    Boot_CI_low  = apply(est, 2, quantile, probs = 0.025, na.rm = TRUE),          # Confidence interval
    Boot_CI_high = apply(est, 2, quantile, probs = 0.975, na.rm = TRUE),
    n_boot_ok    = colSums(!is.na(est))                                           # boot ok
  )
}


#----------- Section 5: Cluster helpers -------------------------------

# save an rds file atomically: first write to a temporary file, then rename
# if a job is killed during saveRDS (time limit), only the .tmp file is broken -> the final file
# either exists completely or not at all (a restart then simply computes it again)
save_rds_atomic <- function(object, file) {
  tmp <- paste0(file, ".tmp", Sys.getpid())
  saveRDS(object, tmp)
  if (!file.rename(tmp, file)) stop("could not rename ", tmp, " to ", file)
  invisible(file)
}

# check that existing results were produced with the same code and settings
# (on the cluster existing files are skipped after a restart -> without this check old results
#  would be mixed with new ones after a change of the code or of the settings)
# sig:  list with everything the results depend on
# file: settings file, one per scenario (-> array tasks never write the same file)
# file missing -> sig is saved; same sig -> ok; different sig -> stop (nothing is deleted automatically)
check_settings <- function(sig, file, what = "Results") {
  if (file.exists(file)) {
    old <- readRDS(file)
    if (!isTRUE(all.equal(old, sig))) {
      stop(what, " in ", dirname(file), " were produced with other code/settings (see ", basename(file), ").\n",
           "  Delete the old files and ", basename(file), " before starting a new run.")
    }
  } else {
    save_rds_atomic(sig, file)
  }
  invisible(TRUE)
}


