###############################################################################
### Overview of scenarios and simulated data
### - Section 1: Table of scenario parameters and true values
### - Section 2: Table with key figures of the data (one run per scenario)
### - Section 3: Plots
###############################################################################

library(tidyverse)
library(dplyr)

sim_log        <- read_csv("./results/00_simulation_log.csv", show_col_types = FALSE)
true_estimands <- read_csv("./results/00_true_estimands.csv", show_col_types = FALSE)


#####################################################################
#----------- Section 1: Scenario parameters and true values --------#
#####################################################################

scenario_table <- sim_log %>%
  filter(sim_run == 1) %>%
  dplyr::select(scenario_id, scenario_name, n_pts, true_beta_arm, drop_beta_arm,
         mortality_level, eff_frailty_y, eff_frailty_cens,
         eff_frailty_train, eff_frailty_contact, report_prob) %>%
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
    method   = ifelse(mortality_level == "low", "LWYY", "Ghosh-Lin")
  )

print(scenario_table, width = Inf)
write_csv(scenario_table, "./results/01_scenario_overview.csv")


#####################################################################
#----------- Section 2: Key figures of the simulated data ----------#
#####################################################################

# load first run of each scenario
first_runs <- sim_log %>% filter(sim_run == 1)

all_data <- NULL
for (i in 1:nrow(first_runs)) {
  dat_i <- readRDS(first_runs$file_path[i])
  dat_i$scenario <- sprintf("S%02d", first_runs$scenario_id[i])
  all_data <- bind_rows(all_data, dat_i)
}

all_data <- all_data %>%
  filter(t > 0) %>%
  mutate(group = ifelse(arm == 1, "Training", "Control"))

# one row per patient
patients <- all_data %>%
  group_by(scenario, group, ID) %>%
  summarise(
    frailty      = first(frailty),
    dropped      = first(ever_dropped),
    died         = first(ever_died),
    weeks_obs    = sum(observed),                     # weeks observed in study
    falls_full   = sum(Yfull),                        # all falls (complete data)
    falls_obs    = sum(Yrep[observed], na.rm=TRUE),               # falls observed in study
    train_rate   = sum(Train[observed], na.rm=TRUE)   / weeks_obs,
    contact_rate = sum(Contact[observed], na.rm=TRUE) / weeks_obs,
    .groups = "drop"
  ) %>%
  # proxy: training in the training group, answered calls in the control group
  mutate(proxy_rate = ifelse(group == "Training", train_rate, contact_rate))

data_table <- patients %>%
  group_by(scenario, group) %>%
  summarise(
    n                   = n(),
    dropout_pct         = round(100 * mean(dropped), 1),
    death_pct           = round(100 * mean(died), 1),
    weeks_obs_mean      = round(mean(weeks_obs), 1),
    falls_per_pt_full   = round(mean(falls_full), 2),
    falls_per_pt_obs    = round(mean(falls_obs), 2),
    share_falls_obs_pct = round(100 * sum(falls_obs) / sum(falls_full), 1),
    proxy_rate_mean     = round(mean(proxy_rate), 2),
    cor_proxy_frailty   = round(cor(proxy_rate, frailty), 2),   # correlation of proxy
    .groups = "drop"
  )

print(data_table, n = Inf, width = Inf)
write_csv(data_table, "./results/01_data_overview.csv")


#####################################################################
#----------- Section 3: Plots --------------------------------------#
#####################################################################

colors_arm <- c("Control" = "#D55E00", "Training" = "#0072B2")

# Plot 1: dropout and death per scenario and arm
data_table %>%
  dplyr::select(scenario, group, dropout_pct, death_pct) %>%
  pivot_longer(c(dropout_pct, death_pct), names_to = "type", values_to = "percent") %>%
  mutate(type = ifelse(type == "dropout_pct", "Dropout", "Death")) %>%
  ggplot(aes(x = scenario, y = percent, fill = group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = colors_arm) +
  labs(title = "Dropout and Death per sceanrio", x = "Scenario", y = "percentage (%)", fill = "arm") +
  theme_minimal(base_size = 12)

# Plot 2: patients still in study over time
all_data %>%
  filter(observed) %>%
  group_by(scenario, group, t) %>%
  summarise(n_active = n(), .groups = "drop") %>%
  group_by(scenario, group) %>%
  mutate(prop_active = n_active / first(n_active)) %>%
  ggplot(aes(x = t, y = prop_active, color = group)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ scenario, ncol = 4) +
  scale_color_manual(values = colors_arm) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Time in study", x = "week", y = "perc active", color = "arm") +
  theme_minimal(base_size = 11)

# Plot 3: mean cumulative falls, complete vs. observed data
# per week: falls / patients at risk, then summed up over the weeks
mcf_full <- all_data %>%
  group_by(scenario, group, t) %>%
  summarise(rate = sum(Yfull, na.rm = TRUE) / n(), .groups = "drop") %>%
  mutate(data = "vollständig")

mcf_obs <- all_data %>%
  filter(observed) %>%
  group_by(scenario, group, t) %>%
  summarise(rate = sum(Yrep, na.rm = TRUE) / n(), .groups = "drop") %>%
  mutate(data = "beobachtet")

bind_rows(mcf_full, mcf_obs) %>%
  group_by(scenario, group, data) %>%
  arrange(t) %>%
  mutate(mcf = cumsum(rate)) %>%
  ggplot(aes(x = t, y = mcf, color = group, linetype = data)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ scenario, ncol = 4) +
  scale_color_manual(values = colors_arm) +
  labs(title = "Mittlere kumulierte Stürze", x = "Woche", y = "Stürze pro Patient",
       color = "Gruppe", linetype = "Daten") +
  theme_minimal(base_size = 11)

# Plot 4: distribution of falls per patient (complete data)
patients %>%
  ggplot(aes(x = falls_full, fill = group)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~ scenario, ncol = 4) +
  scale_fill_manual(values = colors_arm) +
  labs(title = "Falls per patient in 52 weeks (full data))",
       x = "Number of Falls", y = "Patients", fill = "arm") +
  theme_minimal(base_size = 11)

# Plot 5: proxy vs. frailty (how good is the proxy?)
patients %>%
  ggplot(aes(x = frailty, y = proxy_rate, color = group)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(se = FALSE, method = "loess") +
  facet_wrap(~ scenario, ncol = 4) +
  scale_color_manual(values = colors_arm) +
  labs(title = "Proxy (Training- bzw. Kontaktrate) vs. Frailty",
       x = "Frailty (unbeobachtet)", y = "Proxy-Rate pro Woche", color = "Gruppe") +
  theme_minimal(base_size = 11)
