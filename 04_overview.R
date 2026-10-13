###############################################################################
### Overview of scenarios and simulated data
### - Section 1: Table of scenario parameters and true values
### - Section 2: Check distribution in one data set
###############################################################################

library(tidyverse)
library(dplyr)
library(gt)

sim_log        <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)


#####################################################################
#----------- Section 1: Scenario parameters and true values --------#
#####################################################################

scenario_table <- sim_log %>%
  filter(sim_run == 1) %>%
  dplyr::select(scenario_id, scenario_name, n_pts, true_beta_arm, drop_beta_arm,
         mortality_level, eff_frailty_y, eff_frailty_cens,
         eff_frailty_drop,                                       # NEU 07.10.2026
         eff_frailty_train,
         Y_intercept, drop_intercept, death_intercept, train_intercept) %>%       # NEU 07.10.2026: calibrated intercepts
  left_join(
    true_estimands %>%
      dplyr::select(Scenario_Name, Mortality, True_Beta_Arm,
             True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath),
    by = c("scenario_name"   = "Scenario_Name",
           "mortality_level" = "Mortality",
           "true_beta_arm"   = "True_Beta_Arm")
  ) %>%
  mutate(
    scenario = sprintf("S%02d", scenario_id),                    # short name for plots
    method   = ifelse(mortality_level == "high", "Ghosh-Lin", "LWYY")   # NEU 07.10.2026: "none" -> LWYY
  )

print(scenario_table, width = Inf)
#write_csv(scenario_table, "./results/01_scenario_overview.csv")

# # save as html
# scenario_gt <- scenario_table %>%
#   dplyr::select(scenario, scenario_name, method, n_pts, true_beta_arm, drop_beta_arm,
#                 mortality_level, report_prob, eff_frailty_y, eff_frailty_cens,
#                 eff_frailty_drop,                                    
#                 eff_frailty_train, True_RateRatio, True_LogRateRatio,
#                 True_LogRateRatio_NoDeath) %>%
#   gt() %>%
#   fmt_number(columns = c(true_beta_arm, drop_beta_arm, eff_frailty_y, eff_frailty_cens,
#                          eff_frailty_drop,
#                          eff_frailty_train, True_RateRatio, True_LogRateRatio,
#                          True_LogRateRatio_NoDeath),
#              decimals = 3) %>%
#   fmt_percent(columns = report_prob, decimals = 0) %>%
#   cols_label(
#     scenario                  = "Szenario",
#     scenario_name             = "Name",
#     method                    = "Methode",
#     n_pts                     = "n",
#     true_beta_arm             = "β (DGP)",
#     drop_beta_arm             = "β Drop",
#     mortality_level           = "Mortalität",
#     report_prob               = "Report-Wkt.",
#     eff_frailty_y             = "Frailty Y",
#     eff_frailty_cens          = "Frailty Tod",       
#     eff_frailty_drop          = "Frailty Dropout",   
#     eff_frailty_train         = "Frailty Train",
#     True_RateRatio            = "RR",
#     True_LogRateRatio         = "log RR",
#     True_LogRateRatio_NoDeath = "log RR (ohne Tod)"
#   ) %>%
#   tab_spanner(label = "Frailty-Effekte",
#               columns = starts_with("eff_frailty")) %>%
#   tab_spanner(label = "Wahre Schätzgrößen",
#               columns = c(True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath)) %>%
#   tab_header(title = "Simulationsszenarien")
# 

#gtsave(scenario_gt, "res_tables/scenario_table.html")



#####################################################################
#----- Section 2Check of the distributions in one data set ---------#
#####################################################################

# needs: Odat_all (weekly data of one simulated study, complete data until death)

library(tidyverse)
Odat_all <- current_file <- "./results/Odat_scen30_run0001.rds"   # 2_Informative, beta = -0.4, beta_C = -1, high mortality
Odat_all <- readRDS(current_file)

# one row per person
per <- Odat_all %>%
  group_by(ID) %>%
  summarise(arm          = first(arm),
            age          = first(age),
            sex          = first(sex),
            base_risk    = first(base_risk),               # falls in the year before (L1)
            frailty      = first(frailty),
            falls_full   = sum(n_falls_full),              # all falls (complete data)
            falls_obs    = sum(n_falls, na.rm = TRUE),     # observed falls (until dropout)
            dropped      = first(ever_dropped),
            died         = first(ever_died),
            dropout_time = first(dropout_time),
            death_time   = first(death_time),
            train_share  = mean(Train),                    # share of weeks with training
            .groups = "drop") %>%
  mutate(Arm     = factor(arm, 0:1, c("Control", "Training")),
         frail_q = factor(ntile(frailty, 4), 1:4, c("Q1 (least frail)", "Q2", "Q3", "Q4 (most frail)")))

#--- Table 1: overview per arm ----------------------------------------------------------
tab_overview <- per %>%
  group_by(Arm) %>%
  summarise(n               = n(),
            age_mean        = mean(age),
            male_pct        = 100 * mean(sex),
            L1_0_pct        = 100 * mean(base_risk == 0),
            L1_1_pct        = 100 * mean(base_risk == 1),
            L1_2plus_pct    = 100 * mean(base_risk >= 2),
            fall_1plus_pct  = 100 * mean(falls_full >= 1),      # target control arm: ~30 %
            fall_2plus_pct  = 100 * mean(falls_full >= 2),
            falls_mean_full = mean(falls_full),
            falls_mean_obs  = mean(falls_obs),
            falls_max       = max(falls_full),
            dropout_pct     = 100 * mean(dropped),              # target control arm: ~20 %
            death_pct       = 100 * mean(died),
            train_share     = mean(train_share)) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
print(as.data.frame(tab_overview))

#--- Table 2: how strongly do few persons dominate? -------------------------------------
# share of all falls caused by the x % persons with the most falls
tab_concentration <- per %>%
  group_by(Arm) %>%
  arrange(desc(falls_full), .by_group = TRUE) %>%
  summarise(top1_pct  = 100 * sum(falls_full[seq_len(ceiling(0.01 * n()))]) / sum(falls_full),
            top5_pct  = 100 * sum(falls_full[seq_len(ceiling(0.05 * n()))]) / sum(falls_full),
            top10_pct = 100 * sum(falls_full[seq_len(ceiling(0.10 * n()))]) / sum(falls_full)) %>%
  mutate(across(where(is.double), ~ round(.x, 1)))
cat("\nShare of all falls (%) caused by the top 1 / 5 / 10 % of persons:\n")
print(as.data.frame(tab_concentration))

#--- Table 3: events per frailty quartile ------------------------------------------------
tab_frailty <- per %>%
  group_by(Arm, frail_q) %>%
  summarise(fall_1plus_pct = 100 * mean(falls_full >= 1),
            falls_mean     = mean(falls_full),
            dropout_pct    = 100 * mean(dropped),
            death_pct      = 100 * mean(died),
            train_share    = mean(train_share), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
print(as.data.frame(tab_frailty))

#--- Figure 1: number of falls per person ------------------------------------------------
fig_falls <- per %>%
  mutate(falls_cat = factor(pmin(falls_full, 5), 0:5, c(0:4, "5+"))) %>%
  count(Arm, falls_cat) %>%
  group_by(Arm) %>% mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(falls_cat, pct, fill = Arm)) +
  geom_col(position = "dodge") +
  labs(x = "Falls per person in 52 weeks (complete data)", y = "Persons (%)",
       title = "Distribution of the number of falls")

#--- Figure 2: falls in the year before (L1) as proxy -------------------------------------
fig_L1 <- per %>%
  mutate(L1_cat = factor(pmin(base_risk, 4), 0:4, c(0:3, "4+"))) %>%
  group_by(Arm, L1_cat) %>%
  summarise(falls_mean = mean(falls_full), n = n(), .groups = "drop") %>%
  ggplot(aes(L1_cat, falls_mean, colour = Arm, group = Arm)) +
  geom_line() + geom_point(aes(size = n)) +
  labs(x = "Falls in the year before the study (L1)", y = "Mean falls in the study",
       size = "Persons", title = "Proxy L1: falls in the study by falls in the year before")

#--- Figure 3: falls, dropout, death and training by frailty quartile ---------------------
fig_frailty <- tab_frailty %>%
  pivot_longer(c(fall_1plus_pct, dropout_pct, death_pct), names_to = "measure", values_to = "pct") %>%
  mutate(measure = recode(measure, fall_1plus_pct = ">= 1 fall", dropout_pct = "Dropout", death_pct = "Death")) %>%
  ggplot(aes(frail_q, pct, colour = Arm, group = Arm)) +
  geom_line() + geom_point() +
  facet_wrap(~ measure, scales = "free_y") +
  labs(x = "Frailty quartile", y = "Persons (%)", title = "Events by frailty quartile") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#--- Figure 4: cumulative dropout and death over time -------------------------------------
fig_time <- expand_grid(per %>% select(Arm, dropout_time, death_time), week = 1:52) %>%
  group_by(Arm, week) %>%
  summarise(Dropout = 100 * mean(!is.na(dropout_time) & dropout_time <= week),
            Death   = 100 * mean(!is.na(death_time)   & death_time   <= week), .groups = "drop") %>%
  pivot_longer(c(Dropout, Death), names_to = "event", values_to = "pct") %>%
  ggplot(aes(week, pct, colour = Arm)) +
  geom_step() +
  facet_wrap(~ event, scales = "free_y") +
  labs(x = "Week", y = "Cumulative share (%)", title = "Dropout and death over time")

#--- Figure 5: training participation (training arm only) --------------------------------
fig_train <- per %>%
  filter(arm == 1) %>%
  ggplot(aes(frailty, train_share)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, formula = y ~ x) +
  labs(x = "Frailty (baseline)", y = "Share of weeks with training",
       title = "Training participation by frailty (training arm)")

print(fig_falls); print(fig_L1); print(fig_frailty); print(fig_time); print(fig_train)

