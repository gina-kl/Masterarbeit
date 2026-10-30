###############################################################################
### Functions for the analysis (02, 03, 06) ###################################
### - Section 1: Data preparation
### - Section 2: Models (LWYY, Ghosh-Lin)
### - Section 3: Weights
### - Section 4: Mean number of falls (MCF)
### - Section 5: Bootstrap
### - Section 6: Cluster helpers
### functions for the data generation and the true values: functions_dgp.R
###############################################################################

#------------- Section 1: Data preparation -----------------------------------

# Ties: falls of different persons at exactly the same time
# Different handled by coxph and mets
# Every fall time of a person is moved by the same small person-specific offset (0, eps) to the left
# - deterministic, not related to the arm (golden ratio hash of ID)
# - same offset for all falls of a person -> order within a person stays the same
# - risk set is constant within a week (death/censoring only at t + 0.001) -> Breslow estimate unchanged
tie_offset <- function(ID, eps = 0.01) eps * ((ID * 0.6180339887498949) %% 1)


# variables for the analysis, computed from the simulated data (same in 02, 03 and the bootstrap)
prep_odat <- function(dat) {
  dat %>%
    mutate(log_base_falls = log1p(base_risk),                       # log(1 + falls in the year before)
           tstop          = t,                                      # time variable for the weight models
           log_falls_lag  = log1p(cum_falls_lag)) %>%               # log(1 + observed falls up to week t-1)
    arrange(ID, t)
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


#--------------- Section 2: Models (LWYY, Ghosh-Lin) --------------------------------

# Wald confidence interval, test and comparison with the true value
# (model rows: robust SE of the model, MCF rows: robust SE from mets, see mcf_mets())
wald_results <- function(model_name, est, se, true_beta) {
  tibble(
    Model      = model_name,
    True_Beta  = true_beta,
    Estimate   = est,                            # log rate ratio (LWYY) / log mean ratio (GL, MCF)
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

# extract result out of cox model (LWYY) or recreg model (Ghosh-Lin)
# true_beta: the true value that belongs to this method (True_Beta_LWYY or True_Beta_GL)
extract_model_results <- function(model, model_name, true_beta) {
  est <- unname(coef(model)["arm"])
  se  <- sqrt(vcov(model)["arm", "arm"])         # robust (cluster) SE
  wald_results(model_name, est, se, true_beta)
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

    # recreg optimiser (lava::NR): tol 1e-12 (default 1e-6 stops too early) and Armijo line search.
    # Without line search, Newton-Raphson diverged (NaN) in rare data sets with one very frail
    # person (50-60 falls), mostly GL 0. Armijo halves a step only if it does not decrease the
    # objective enough -> same solution, normal fits unchanged
    gl_ctrl <- list(tol = 1e-12, backtrack = "armijo")

    # GL: cens.model = ~ strata(arm) -> internal censoring Kaplan-Meier per arm
    # 6 / 7: truncated weights (sensitivity analysis, see compute_weights())
    list(
      "LWYY 0: Full data (oracle M0)" = coxph(f_lwyy, data = cp_full, ties = "breslow"),
      "LWYY 1: Naive (Unweighted)"    = coxph(f_lwyy, data = cp, ties = "breslow"),
      "LWYY 3: IPCW sym"              = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_sym),
      "LWYY 4: IPCW sym_train"        = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_symtrain),
      "LWYY 5: IPCW oracle"           = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_oracle),
      "LWYY 6: IPCW sym trunc"        = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_sym_tr),
      "LWYY 7: IPCW sym_train trunc"  = coxph(f_lwyy, data = cp, ties = "breslow", weights = cp$w_symtrain_tr),
      "GL 0: Full data (oracle M0)"   = recreg(f_gl, data = cp_full, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), control = gl_ctrl),
      "GL 1: Naive (Unweighted)"      = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm), control = gl_ctrl),
      "GL 3: IPCW sym"                = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_sym, control = gl_ctrl),
      "GL 4: IPCW sym_train"          = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_symtrain, control = gl_ctrl),
      "GL 5: IPCW oracle"             = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_oracle, control = gl_ctrl),
      "GL 6: IPCW sym trunc"          = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_sym_tr, control = gl_ctrl),
      "GL 7: IPCW sym_train trunc"    = recreg(f_gl, data = cp, cause = 1, death.code = 3, cens.code = 0, cens.model = ~ strata(arm),
                                               weights = cp$w_symtrain_tr, control = gl_ctrl)
    )
}

# ----------------- Section 3: Weights -------------------------------------

# Compute ipw weights
# IPCW weights per arm (for 02_analyse_weekly.R. 03_mehrfach analyse and bootstrap)
# iptwm fitted on weeks without death, lag_ipw shifts weigths by one week (falls of week t are observed before a dropout at its end)
# w_* stabilised weights, wm_* for marginal mean (only time in numerator)
# *_tr: truncated at the trunc_q percentile (sensitivity analysis, only the two proxy models)
# dat: observed weekly data (observed == TRUE) with tstop, log_base_falls, log_falls_lag, train_rate_lag, frailty_t


compute_weights <- function(dat, trunc_q = 0.99){
  dat %>%
    mutate(tstart = t - 1) %>%
    split(.$arm) %>%
    lapply(function(d_all){
      d_all <- d_all %>% arrange(ID, t) %>% as.data.frame()
      d_fit <- d_all %>% filter(death == 0)          # weeks in which a dropout is possible
      shift <- function(fit)lag_ipw(fit, d_fit, d_all) # shift ipcw by one week

      w_sym <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
        numerator = ~ age + sex + log_base_falls + ns(tstop, df = 3),
        denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + log_falls_lag,
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
                     denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + log_falls_lag,
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
                                  denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + log_falls_lag + train_rate_lag,
                                  id = ID, tstart = tstart, timevar = tstop, type = "cens", data = d_fit))
        # unstabilisied weights
        wm_symtrain <- shift(ipwtm(exposure = drop, family = "binomial", link = "logit",
                                   numerator   = ~ ns(tstop, df = 3),
                                   denominator = ~ age + sex + log_base_falls + ns(tstop, df = 3) + log_falls_lag + train_rate_lag,
                                   id = ID, tstart = tstart, timevar = tstop, type = "cens", data = d_fit))
      } else {                                       # control arm never trains -> same as sym
        w_symtrain  <- w_sym
        wm_symtrain <- wm_sym
      }
      d_all %>% mutate(w_sym = w_sym, w_symtrain = w_symtrain, w_oracle = w_oracle,
                       wm_sym = wm_sym, wm_symtrain = wm_symtrain, wm_oracle = wm_oracle,
                       # truncated per arm (sensitivity analysis), only the two proxy models
                       w_sym_tr  = trunc_w(w_sym, trunc_q),  w_symtrain_tr  = trunc_w(w_symtrain, trunc_q),
                       wm_sym_tr = trunc_w(wm_sym, trunc_q), wm_symtrain_tr = trunc_w(wm_symtrain, trunc_q))
    }) %>%
    bind_rows() %>%
    arrange(ID, t)
    }


# truncate weights at the q-th percentile of all person-weeks of one arm (Cole & Hernán 2008), only from above:
# the problem are single very large weights, small weights (< 1) hardly change the variance
trunc_w <- function(w, q) pmin(w, quantile(w, q, names = FALSE))


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


# ----------------- Section 4: Mean number of falls (MCF) ------------------

# Mean number of falls per arm and period (nonparametric Ghosh-Lin / Cook-Lawless, weighted)
#   mu_a = sum_t S_a(t-1) * dR_a(t)
#   dR_a(t): weighted mean number of falls in week t among persons at risk (Nelson-Aalen)
#   S_a(t) : weighted Kaplan-Meier for death (death at the end of week t)
#   Use unstabilised weights (or weights whose numerator depends only on arm and time):
#   a numerator with covariates would change the covariate mix of the pseudo-population.
# d: weekly data, one row per person and week at risk, columns arm, t, falls, death, weight
# breaks: limits of the periods in weeks, c(0, 52) = whole year (one row per arm),
#         c(0, 13, 26, 39, 52) = quarters (period 1-4)
weighted_mcf <- function(d, w = "w", falls = "n_falls", death = "death", breaks = c(0, 52)) {
  d %>%
    mutate(.w = .data[[w]], .y = .data[[falls]], .dth = .data[[death]]) %>%
    group_by(arm, t) %>%
    summarise(dR = sum(.w * .y) / sum(.w),          # fallsrate weighted
              hD = sum(.w * .dth) / sum(.w),        # death hazard weighted
              .groups = "drop") %>%
    group_by(arm) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(S_lag  = c(1, head(cumprod(1 - hD), -1)),       # alive at start of week t
           period = cut(t, breaks, labels = FALSE)) %>%    # number of the period (NA outside the breaks)
    filter(!is.na(period)) %>%
    group_by(arm, period) %>%
    summarise(MCF = sum(S_lag * dR), .groups = "drop")      # mean falls in the real world (Ghosh-Lin)
}


# Model-based mean number of falls up to t.end per arm (standardisation)
# prediction for every person with arm set to 0 / 1, averaged over all persons
# base_cov: one row per person with the covariates of the model (age, sex, log_base_falls)
# Ghosh-Lin (recreg): mu_0(t) * exp(X beta)     -> real world with death (compare with MCF_*)
# LWYY (coxph):       Lambda_0(t) * exp(X beta) -> compare with the limit of LWYY on complete data (MCF_LWYY_*,
#                     true_block() in functions_dgp.R); equals the world without death (MCF_*_NoDeath) only if
#                     death is independent of the falls given the covariates (here not: death depends on frailty)
# log(MCF_1 / MCF_0) is exactly the coefficient of arm; new is only the level per arm
# basehaz(centered = TRUE): Lambda_0 at the means of the covariates, centered = FALSE can over-/underflow
model_mcf <- function(fit, base_cov, t.end = 52) {
  if (inherits(fit, "coxph")) {                              # LWYY
    H0   <- basehaz(fit, centered = TRUE)                    # Lambda_0(t) at the centering values zbar = fit$means, Breslow estimator
    H0_t <- max(H0$hazard[H0$time <= t.end])                 # (centered = FALSE can over-/underflow, see ?basehaz)
    b    <- coef(fit)
    zbar <- fit$means[names(b)]
    sapply(0:1, function(a) {
      X <- as.matrix(base_cov %>% mutate(arm = a) %>% dplyr::select(all_of(names(b)))) # base_cov = ID, age, sex, log_base_falls
      H0_t * mean(exp(sweep(X, 2, zbar) %*% b))              # Lambda_0(tau; zbar) * mean(exp((X_i - zbar) b))
    })
  } else {                                                   # Ghosh-Lin
    sapply(0:1, function(a)
      mean(predict(fit, newdata = base_cov %>% mutate(arm = a), times = t.end, se = FALSE)$cumhaz))
  }
}

# mean number of falls per arm: nonparametric (weighted_mcf) and, if fits are given,
# model-based from all LWYY and GL fits of fit_models() (model_mcf)
# result: list, one element c(control, intervention) per version
mcf_estimate <- function(Odat_all, Odat_w, fits = NULL){
  # non-parametric estimate
  mcf_list <- list(
    "MCF 0: Full data (oracle M0)" = weighted_mcf(Odat_all %>% mutate(w1 = 1), "w1", "n_falls_full", "death_full"), # w1 = 1: unweighted, Estimate without dropout
    "MCF 1: Naive (Unweighted)"    = weighted_mcf(Odat_w %>% mutate(w1 = 1), "w1"),
    "MCF 3: IPCW sym"              = weighted_mcf(Odat_w, "wm_sym"),
    "MCF 4: IPCW sym_train"        = weighted_mcf(Odat_w, "wm_symtrain"),
    "MCF 5: IPCW oracle"           = weighted_mcf(Odat_w, "wm_oracle"),
    "MCF 6: IPCW sym trunc"        = weighted_mcf(Odat_w, "wm_sym_tr"),
    "MCF 7: IPCW sym_train trunc"  = weighted_mcf(Odat_w, "wm_symtrain_tr")
  )
  base_cov <- Odat_all %>% filter(t == 1) %>% dplyr::select(ID, age, sex, log_base_falls)   # one row per person
  c(
    lapply(mcf_list, function(mm) c(mm$MCF[mm$arm == 0], mm$MCF[mm$arm == 1])),
    if (!is.null(fits)) setNames(lapply(fits, function(f) tryCatch(model_mcf(f, base_cov), error = function(e) c(NA, NA))),
                                 paste0("MCF-", names(fits)))                                  # "MCF-LWYY ..." / "MCF-GL ..."
  )
}


# Nonparametric mean number of falls up to t.end per arm with mets (Ghosh & Lin 2000), with robust SE
# same estimator as weighted_mcf(): the risk set is constant within a week, unweighted the results are identical,
# weighted they differ by ~1e-5 (tie_offset moves single rows into the next week)
# weights are treated as known (as for the robust SE of the weighted LWYY / GL models)
# cp: counting process data from prep_counting(), w: weight column (NULL = unweighted)
# status 1 = fall (cause), 3 = death (terminal event), strata(arm) -> one estimate per arm
mets_mcf <- function(cp, w = NULL, t.end = 52) {
  wt  <- if (is.null(w)) rep(1, nrow(cp)) else cp[[w]]
  fit <- recurrent_marginal(Event(start, stop, status) ~ strata(arm) + cluster(ID), data = cp,
                            cause = 1, death.code = 3, weights = wt)
  s <- summary(fit, times = t.end)$pbaseci   # one data frame per arm (0, 1)
  c(mu_0 = s[[1]]$mean, mu_1 = s[[2]]$mean, se_0 = s[[1]]$se, se_1 = s[[2]]$se)
}

# nonparametric MCF with mets for every version (same names as mcf_estimate()):
# - mean falls per arm, check against weighted_mcf (unweighted identical, weighted ~1e-5 because of tie_offset)
# - robust SE of the log ratio log(mu_1 / mu_0), delta method, arms independent: se^2 = (se_1 / mu_1)^2 + (se_0 / mu_0)^2
mcf_mets <- function(Odat_all, Odat_w) {
  cp      <- prep_counting(Odat_w, "observed")
  cp_full <- prep_counting(Odat_all, "full")
  m <- list(
    "MCF 0: Full data (oracle M0)" = mets_mcf(cp_full),
    "MCF 1: Naive (Unweighted)"    = mets_mcf(cp),
    "MCF 3: IPCW sym"              = mets_mcf(cp, "wm_sym"),
    "MCF 4: IPCW sym_train"        = mets_mcf(cp, "wm_symtrain"),
    "MCF 5: IPCW oracle"           = mets_mcf(cp, "wm_oracle"),
    "MCF 6: IPCW sym trunc"        = mets_mcf(cp, "wm_sym_tr"),
    "MCF 7: IPCW sym_train trunc"  = mets_mcf(cp, "wm_symtrain_tr")
  )
  tibble(Model                 = names(m),
         MCF_Control_mets      = sapply(m, `[[`, "mu_0"),
         MCF_Intervention_mets = sapply(m, `[[`, "mu_1"),
         SE_log_ratio          = sapply(m, function(x) sqrt((x[["se_1"]] / x[["mu_1"]])^2 + (x[["se_0"]] / x[["mu_0"]])^2)))
}


#----------- Section 5: Bootstrap (06_bootstrap_cluster.R) -------------------------------------

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

# estimates of one bootstrap sample: weights are estimated again, then only the nonparametric MCF
# (LWYY / GL have a robust SE) -> log ratio per version
boot_estimate_mcf <- function(Odat_all) {
  Odat_w <- compute_weights(Odat_all %>% filter(observed))
  mcf    <- mcf_estimate(Odat_all, Odat_w)                   # without fits -> only nonparametric rows
  sapply(mcf, function(mu) log(mu[2] / mu[1]))
}

# summary over the bootstrap samples (est: matrix, one row per successful sample)
summarise_boot <- function(est) {
  tibble(
    Model        = colnames(est),
    Boot_SE      = apply(est, 2, sd, na.rm = TRUE),                               # Bootstrap SE
    Boot_CI_low  = apply(est, 2, quantile, probs = 0.025, na.rm = TRUE),          # percentile confidence interval
    Boot_CI_high = apply(est, 2, quantile, probs = 0.975, na.rm = TRUE),
    n_boot_ok    = colSums(!is.na(est))                                           # successful samples
  )
}


#----------- Section 6: Cluster helpers -------------------------------

# save an rds file atomically: first write to a temporary file, then rename
# if a job is killed during saveRDS (time limit), only the .tmp file is broken -> the final file
# either exists completely or not at all (a restart then simply computes it again)
save_rds_atomic <- function(object, file) {
  tmp <- paste0(file, ".tmp", Sys.getpid())
  saveRDS(object, tmp)
  if (!file.rename(tmp, file)) stop("could not rename ", tmp, " to ", file)
  invisible(file)
}

