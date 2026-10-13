###############################################################################
### Function ##################################################################
### - Section 1: Data preparation
### - Section 2: Function for analyse
###############################################################################

############ Section 1: Data preparation

# Prepare counting-process data for true model
# complete data -> one row per interval between falls 
# k falls in week t are placed at t-1 + j/(k+d), j = 1..k (d = 1 if death at the end of week t),
# i.e. inside week t as in expand_counting(); status: 1 = fall, 3 = death, 0 = end of study
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
    transmute(ID, stop = t - 1 + j / (Yfull + d), status = 1L)
  
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
    mutate(m = pmax(1L, as.integer(k + term_death + term_drop))) %>%
    tidyr::uncount(m, .id = "j", .remove = FALSE) %>%
    mutate(start = t - 1 + (j - 1) / m,
           stop  = t - 1 + j / m) %>%
    group_by(ID) %>%
    mutate(last_row = row_number() == n()) %>%
    ungroup() %>%
    mutate(status = case_when(j <= k                   ~ 1L,   # fall
                              j == m & term_death == 1 ~ 3L,   # death at the end of the week
                              last_row                 ~ 0L,   # dropout or end of study = censoring
                              TRUE                     ~ 9L))  # no event
  
  # last row is a fall: short censoring row at the end of study
  extra <- cp %>%
    filter(last_row & status == 1L) %>%
    mutate(start = stop, status = 0L)
  
  bind_rows(cp, extra) %>%
    mutate(stop = ifelse(status %in% c(0L, 3L), stop + 0.001, stop)) %>%   # ties: falls first
    arrange(ID, start, stop) %>%
    dplyr::select(-k, -m, -j, -term_death, -term_drop, -last_row) %>%
    as.data.frame()
}


############ Section 2: Function for analyse




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


############ Section 3: Weights and mean number of falls

# Compute ipw weights
# IPVW weights per arm (for 02_analyse_weekly.R. 03_mehrfach analyse and bootstrap)
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
#   ipwtm() returns cumulative weights INCLUDING the current week (probability of no dropout
#   at the end of week t). The falls of week t happen BEFORE a dropout at the end of week t,
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
    mutate(w = lag(w_cum, default = 1)) %>%
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
    #  .dth instead of .d: R completes argument names that are a unique abbreviation,
    #   ".d" was read as mutate's first argument ".data" -> error "Can't subset .data outside of a data mask"
    mutate(.w = .data[[w]], .y = .data[[falls]], .dth = .data[[death]]) %>%
    group_by(arm, t) %>%
    summarise(dR = sum(.w * .y) / sum(.w),          # falls per person at risk in week t
              hD = sum(.w * .dth) / sum(.w),        # death hazard in week t   
              .groups = "drop") %>%
    group_by(arm) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(S_lag = c(1, head(cumprod(1 - hD), -1))) %>%   # alive at start of week t
    summarise(MCF = sum((S_lag * dR)[t <= t.end]),         # mean falls in the real world (Ghosh-Lin)
              R   = sum(dR[t <= t.end]),                   # cumulative rate among survivors
              .groups = "drop")
}






# Model-based mean number of falls up to t.end per arm (standardisation)
# prediction for every person with arm set to 0 / 1, averaged over all persons
# base_cov: one row per person with the covariates of the model (age, sex, log_base_falls)
# Ghosh-Lin (recreg): mu_0(t) * exp(X beta)     -> real world with death (compare with MCF_*)
# LWYY (coxph):       Lambda_0(t) * exp(X beta) -> hypothetical world without death (compare with MCF_*_NoDeath),
# only valid if death is independent of the falls given the covariates
# log(MCF_1 / MCF_0) is exactly the coefficient of arm; new is only the level per arm
# LWYY (coxph): Lambda_0(t) * exp(X beta) -> hypothetical world without death (compare with MCF_*_NoDeath),
model_mcf <- function(fit, base_cov, t.end = 52) {
  if (inherits(fit, "coxph")) {                              # LWYY
    H0   <- basehaz(fit, centered = FALSE)                   # cumulative baseline rate Lambda_0(t)
    H0_t <- max(H0$hazard[H0$time <= t.end])
    b    <- coef(fit)
    sapply(0:1, function(a) {
      X <- as.matrix(base_cov %>% mutate(arm = a) %>% dplyr::select(all_of(names(b))))
      H0_t * mean(exp(X %*% b))
    })
  } else {                                                   # Ghosh-Lin
    sapply(0:1, function(a)
      mean(predict(fit, newdata = base_cov %>% mutate(arm = a), times = t.end, se = FALSE)$cumhaz))
  }
}


