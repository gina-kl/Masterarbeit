###############################################################################
### Function ##################################################################
### - Section 1: Function for split events
### - Section 2: Function for analyse
###############################################################################

# load packages
library(tidyverse)

############ Section 1: Function for split events
# Split weeks with k falls in k+1 intervalls (at random times within this week)
#   status_col: "status" or "status_full"
#   n_col:      number of falls ("Yrep" or "Yfull")
expand_falls <- function(dat, status_col, n_col) {
  
  dat$n_ev <- ifelse(dat[[status_col]] == 1, dat[[n_col]], 0)  # falls per column
  
  no_ev <- dat %>% filter(n_ev == 0)         # no events
  
  ev_long <- dat %>%
    filter(n_ev > 0) %>%                       # only for columns with event
    uncount(n_ev + 1) %>%                      # k fall rows + 1 remaining row
    group_by(ID, t) %>%                        # patient per week
    mutate(
      ftime  =  c(first(tstart) +
                    sort(sample(1:999, n() - 1)) / 1000 * (first(tstop) - first(tstart)),
                  first(tstop)),  # Draw using a fine grid so that the elements aren't too close together
      tstart = c(first(tstart), head(ftime, -1)),                                   # start where last interval stops
      tstop  = ftime,
      !!status_col := c(rep(1, n() - 1), 0)    # falls, then rest of week without event
    ) %>%
    ungroup() %>%
    dplyr::select(-ftime)
  
  bind_rows(no_ev, ev_long) %>%
    dplyr::select(-n_ev) %>%                      # drop help column
    arrange(ID, tstart, tstop) %>%
    as.data.frame()
}



############ Section 2: Function for analyse
# extract result out of cox model (LWYY) or recreg model (Ghosh-Lin)
extract_model_results <- function(model, model_name, true_beta, true_beta_nodeath = NA) {

  # extract values
  log_hr <- unname(coef(model)["arm"])
  robust_se <- sqrt(vcov(model)["arm", "arm"])
  p_val <- 2 * pnorm(-abs(log_hr / robust_se))     # Wald test

  # summary in tibble
  tibble(
    Model = model_name,
    True_Beta = true_beta,
    Log_HR = round(log_hr, 4),
    Hazard_Ratio = round(exp(log_hr), 2),  
    Robust_SE = round(robust_se, 4),
    P_Value = round(p_val, 4),
    Significant = ifelse(p_val < 0.05, "Yes", "No"),
    Bias = round(log_hr - true_beta, 4),
    Bias_NoDeath = round(log_hr - true_beta_nodeath, 4)
  )
}


# Prepare counting-process data for Ghosh-Lin
# recreg counts every row with status 0 as censoring. Every
# row without a fall would be a censoring. Therefore new status:
#   1 = fall, 3 = death, 0 = censoring (only in the last row), 9 = week without event
# status_col: "status" (observed data) or "status_full" (complete data)
prep_gl_data <- function(dat, status_col) {

  dat$st <- dat[[status_col]]

  # remove week of dropout (no fall possible) 
  dat <- dat %>% filter(st != 2)

  # mark last row per patient
  dat <- dat %>%
    arrange(ID, tstart) %>%
    group_by(ID) %>%
    mutate(last_row = row_number() == n()) %>%
    ungroup()

  # new status
  dat$status_gl <- 9                # no event
  dat$status_gl[dat$st == 1] <- 1   # fall
  dat$status_gl[dat$st == 3] <- 3   # death
  dat$status_gl[dat$st == 0 & dat$last_row] <- 0 #censoring


  # censoring and death shortly after the falls of the same week (avoids ties)
  shift <- dat$status_gl %in% c(0, 3)
  dat$tstop[shift] <- dat$tstop[shift] + 0.001

  as.data.frame(dat)
}

# Truncate IPW weights (same as trunc in ipwtm):
# weights below the level-percentile and above the (1 - level)-percentile
# are set to these percentiles, e.g. level = 0.01 -> 1st and 99th percentile
trunc_weights <- function(w, level = 0.01) {
  q <- quantile(w, probs = c(level, 1 - level), na.rm = TRUE)
  pmin(pmax(w, q[1]), q[2])
}
