###############################################################################
### Overview of the data generation over all runs (no single data sets)
### - Section 1: Table of scenario parameters and true values
### - Section 2: Realised data per scenario and arm (mean over the runs)
### - Section 3: Mean number of falls against the true values
### - Section 4: Events per frailty quartile
### - Section 5: Figures
### Input: 00_simulation_log.csv with the summary of each data set (summarise_odat, functions_dgp.R)
###        and 00_true_estimands.csv per block
###############################################################################

library(tidyverse)
library(gt)
source("GlobalVars.R")   # target values

# blocks to show: c("main"), c("supp"), c("main", "supp"), c("sens") or c("lokal") (as in 05_ergebnisgrafik.R)
blocks   <- c("mort")
res_dirs <- c(main = "./results_cluster/ergebnisse_main",
              supp = "./results_cluster/ergebnisse_supp",
              mort = "./results_cluster/ergebnisse_mort",
              sens = "./results_cluster/ergebnisse_sens",
              lokal = "./results")[blocks]
out_dir  <- file.path("./05_results", paste(blocks, collapse = "_"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# one file per block (scenario_id is unique over all blocks)
sim_log        <- map(file.path(res_dirs, "00_simulation_log.csv"), read_csv, show_col_types = FALSE) %>% list_rbind() # combine to one data frame
true_estimands <- map(file.path(res_dirs, "00_true_estimands.csv"), read_csv, show_col_types = FALSE) %>%
  list_rbind() %>%
  distinct(Scenario_Name, Mortality, True_Beta_Arm, .keep_all = TRUE)   # only first row of rows with this combination, keep all other rows
if (!"ctrl_fall1" %in% names(sim_log)) stop("no summaries in 00_simulation_log.csv: the data were generated before summarise_odat")

n_runs <- max(sim_log$sim_run)
scenarios <- sim_log %>%
  distinct(scenario_id, scenario_name, effect, drop_beta_arm, mortality_level) %>%
  arrange(scenario_id) %>%
  mutate(scenario = sprintf("S%02d", scenario_id),
         label    = sprintf("%s: %s, %s, β_C = %g, Mortalität %s", scenario, scenario_name, effect, drop_beta_arm, mortality_level),
         label    = factor(label, levels = label))

# long format: one row per run, arm and measure, e.g. ctrl_fall1 -> arm = Kontrolle, measure = fall1
runs_long <- sim_log %>%
  dplyr::select(scenario_id, sim_run, matches("^(ctrl|int)_")) %>%
  pivot_longer(-c(scenario_id, sim_run), names_to = c("arm", "measure"), names_pattern = "^(ctrl|int)_(.*)$") %>%
  mutate(arm = factor(arm, levels = c("ctrl", "int"), labels = c("Kontrolle", "Training")))

# mean and standard deviation over the runs (one row per scenario_id, arm and measure), check of dgp
data_stats <- runs_long %>%
  group_by(scenario_id, arm, measure) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  left_join(scenarios %>% dplyr::select(scenario_id, scenario, label), by = "scenario_id")


#----------- Section 1: Scenario parameters and true values --------#

scenario_table <- sim_log %>%
  filter(sim_run == 1) %>%
  dplyr::select(scenario_id, scenario_name, effect, n_pts, true_beta_arm, max_frailty_red, drop_beta_arm,
                mortality_level, eff_frailty_y, eff_frailty_death, eff_frailty_drop, eff_fall_drop, eff_frailty_train,
                Y_intercept, drop_intercept, death_intercept, train_intercept) %>%
  left_join(true_estimands %>%
              dplyr::select(Scenario_Name, Mortality, True_Beta_Arm, True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath, True_Beta_GL, True_Beta_LWYY),
            by = c("scenario_name" = "Scenario_Name", "mortality_level" = "Mortality", "true_beta_arm" = "True_Beta_Arm")) %>%
  arrange(scenario_id) %>%
  mutate(scenario = sprintf("S%02d", scenario_id), .before = 1) %>%
  dplyr::select(-scenario_id)

scenario_gt <- scenario_table %>%
  gt() %>%
  fmt_number(columns = c(Y_intercept, drop_intercept, death_intercept, train_intercept,
                         True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath), decimals = 3) %>%
  cols_label(scenario = "Szenario", scenario_name = "Name", effect = "Effekt", n_pts = "n", true_beta_arm = "β",
             max_frailty_red = "max. Frailty-Senkung", drop_beta_arm = "β_C", mortality_level = "Mortalität",
             eff_frailty_y = "γ_Y", eff_frailty_death = "γ_D", eff_frailty_drop = "γ_C", eff_fall_drop = "γ_F",
             eff_frailty_train = "γ_Tr", Y_intercept = "α_Y", drop_intercept = "α_C", death_intercept = "α_D",
             train_intercept = "α_Tr", True_RateRatio = "RR", True_LogRateRatio = "log RR",
             True_LogRateRatio_NoDeath = "log RR (ohne Tod)", True_Beta_GL = "true beta GL", True_Beta_LWYY = "true beta LWYY") %>%
  tab_spanner(label = "Effekte der Frailty", columns = c(eff_frailty_y, eff_frailty_death, eff_frailty_drop, eff_fall_drop, eff_frailty_train)) %>%
  tab_spanner(label = "Kalibrierte Intercepts", columns = c(Y_intercept, drop_intercept, death_intercept, train_intercept)) %>%
  tab_spanner(label = "Wahre Werte", columns = c(True_RateRatio, True_LogRateRatio, True_LogRateRatio_NoDeath)) %>%
  tab_header(title = "Simulationsszenarien") %>%
  tab_footnote(
    footnote = "γ_Y: Effect frailty on falls, γ_D: Eff. frailty on Death, γ_C: Eff. frailty on dropout,  γ_F: Eff. falls on dropout, γ_Tr: Effect frailty on Training"
  )

gtsave(scenario_gt, file.path(out_dir, "04_scenario_table.html"))


#----- Section 2: Realised data per scenario and arm ---------------#

pct_measures <- c("male", "L1_0", "L1_1", "L1_2plus", "fall1", "fall2", "fall5", "top10", "dropout", "death", "train_w1", "train_share")

overview <- data_stats %>%
  filter(measure %in% c("n", "age_mean", pct_measures, "falls_mean", "falls_obs_mean", "falls_max")) %>%
  mutate(mean = ifelse(measure %in% pct_measures, 100 * mean, mean)) %>%
  dplyr::select(label, arm, measure, mean) %>%
  pivot_wider(names_from = measure, values_from = mean) %>%
  arrange(label, arm) %>%
  dplyr::select(label, arm, n, age_mean, male, L1_0, L1_1, L1_2plus, fall1, fall2, fall5, falls_mean, falls_obs_mean,
                falls_max, top10, dropout, death, train_w1, train_share)

overview_gt <- overview %>%
  gt(groupname_col = "label", rowname_col = "arm") %>%
  fmt_number(columns = c(n, falls_max), decimals = 0) %>%
  fmt_number(columns = c(age_mean, male, L1_0, L1_1, L1_2plus, fall1, fall2, fall5, top10, dropout, death, train_w1, train_share), decimals = 1) %>%
  fmt_number(columns = c(falls_mean, falls_obs_mean), decimals = 3) %>%
  cols_label(n = "n", age_mean = "Alter", male = "Männer %", L1_0 = "0", L1_1 = "1", L1_2plus = "≥ 2",
             fall1 = "≥ 1", fall2 = "≥ 2", fall5 = "≥ 5", falls_mean = "vollständig", falls_obs_mean = "beobachtet",
             falls_max = "Maximum", top10 = "Top 10 %", dropout = "Dropout %", death = "Tod %",
             train_w1 = "Woche 1", train_share = "Jahr") %>%
  tab_spanner(label = "Stürze im Vorjahr (%)", columns = c(L1_0, L1_1, L1_2plus)) %>%
  tab_spanner(label = "Stürze in der Studie (%)", columns = c(fall1, fall2, fall5)) %>%
  tab_spanner(label = "Stürze je Person", columns = c(falls_mean, falls_obs_mean, falls_max)) %>%
  tab_spanner(label = "Training (%)", columns = c(train_w1, train_share)) %>%
  tab_header(title = "Realisierte Daten je Szenario und Arm",
             subtitle = paste0("Mittelwert über ", n_runs, " Läufe; Maximum = mittleres Maximum je Lauf; Top 10 % = Anteil aller Stürze der 10 % Personen mit den meisten Stürzen")) %>%
  tab_source_note(sprintf("Zielwerte im Kontrollarm: %g %% mit ≥ 1 Sturz (gerechnet ohne Tod), %g %% Dropout, Tod %s %%. Trainingsarm: %g %% Training in Woche 1.",
                          100 * target_falls, 100 * target_dropout, paste(100 * target_death, collapse = " / "), 100 * target_train))

gtsave(overview_gt, file.path(out_dir, "04_data_overview.html"))
write_csv(overview, file.path(out_dir, "04_data_overview.csv"))



#----- Section 3: Mean number of falls against the true values -----#

# check that data and true values belong together: difference in units of its Monte Carlo SE (|z| > 3: look at it)

truth_long <- true_estimands %>%
  dplyr::select(Scenario_Name, Mortality, True_Beta_Arm, MCF_Control, MCF_Intervention) %>%
  pivot_longer(c(MCF_Control, MCF_Intervention), names_to = "arm", values_to = "true_mean") %>%
  mutate(arm = recode(arm, MCF_Control = "Kontrolle", MCF_Intervention = "Training"))

truth_check <- sim_log %>%
  distinct(scenario_id, scenario_name, mortality_level, true_beta_arm) %>%
  inner_join(truth_long, by = c("scenario_name" = "Scenario_Name", "mortality_level" = "Mortality", "true_beta_arm" = "True_Beta_Arm"),
             relationship = "many-to-many") %>%
  inner_join(data_stats %>% filter(measure == "falls_mean") %>% mutate(arm = as.character(arm)), by = c("scenario_id", "arm")) %>%
  mutate(mcse = sd / sqrt(n_runs), z = (mean - true_mean) / mcse) %>%
  arrange(scenario_id, arm) %>%
  dplyr::select(label, arm, mean, true_mean, mcse, z)

truth_gt <- truth_check %>%
  gt(groupname_col = "label", rowname_col = "arm") %>%
  fmt_number(columns = c(mean, true_mean, mcse), decimals = 3) %>%
  fmt_number(columns = z, decimals = 1) %>%
  cols_label(mean = "Mittel der Daten", true_mean = "wahrer Wert", mcse = "Monte-Carlo-SE", z = "z") %>%
  tab_header(title = "Mittlere Sturzzahl: simulierte Daten gegen wahre Werte",
             subtitle = paste0("Vollständige Daten, ", n_runs, " Läufe; z = Abweichung in Monte-Carlo-SE"))

gtsave(truth_gt, file.path(out_dir, "04_truth_check.html"))



#----- Section 4: Events per frailty quartile ----------------------#


frailty_stats <- data_stats %>%
  filter(str_detect(measure, "_fq[1-4]$")) %>%
  mutate(quartile = str_extract(measure, "[1-4]$"),
         quartile = factor(quartile, 1:4, c("Q1 (am wenigsten frail)", "Q2", "Q3", "Q4 (am meisten frail)")),
         measure  = sub("_fq[1-4]$", "", measure),
         mean     = ifelse(measure == "falls_mean", mean, 100 * mean))

frailty_table <- frailty_stats %>%
  dplyr::select(label, arm, quartile, measure, mean) %>%
  pivot_wider(names_from = measure, values_from = mean) %>%
  arrange(label, arm, quartile) %>%
  dplyr::select(label, arm, quartile, fall1, falls_mean, dropout, death, train_share)

frailty_gt <- frailty_table %>%
  gt(groupname_col = "label") %>%
  fmt_number(columns = c(fall1, dropout, death, train_share), decimals = 1) %>%
  fmt_number(columns = falls_mean, decimals = 2) %>%
  cols_label(arm = "Arm", quartile = "Frailty-Quartil", fall1 = "≥ 1 Sturz %", falls_mean = "Stürze je Person",
             dropout = "Dropout %", death = "Tod %", train_share = "Training %") %>%
  tab_header(title = "Ereignisse je Quartil der Baseline-Frailty", subtitle = paste0("Mittelwert über ", n_runs, " Läufe"))

gtsave(frailty_gt, file.path(out_dir, "04_frailty_quartiles.html"))



#----- Section 5: Figures ------------------------------------------#


theme_ov <- list(theme_minimal(base_size = 10),
                 theme(panel.grid.minor = element_blank(), legend.position = "top",
                       strip.text = element_text(face = "bold")),
                 scale_color_manual(values = c("Kontrolle" = "#D55E00", "Training" = "#0072B2")),
                 scale_fill_manual(values = c("Kontrolle" = "#D55E00", "Training" = "#0072B2")))
width_ov <- min(40, 3 + 2.6 * nrow(scenarios))   # one column per scenario; meant for main or sens, supp has too many columns

#--- Figure 1: falls per person ----------------------------------------------------------
fig_falls <- data_stats %>%
  filter(measure %in% c("fall1", "fall2", "fall5")) %>%
  mutate(falls = factor(measure, c("fall1", "fall2", "fall5"), c("≥ 1", "≥ 2", "≥ 5"))) %>%
  ggplot(aes(falls, 100 * mean, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", 100 * mean)), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.6, color = "grey30") +
  facet_wrap(~ scenario, nrow = 1) +
  labs(x = "Stürze je Person in 52 Wochen (vollständige Daten)", y = "Personen (%)", fill = NULL,
       title = "Verteilung der Sturzzahl je Arm", subtitle = paste0("Mittelwert über ", n_runs, " Läufe")) +
  theme_ov

#--- Figure 2: events per frailty quartile ----------------------------------------------

fig_frailty <- frailty_stats %>%
  filter(measure %in% c("fall1", "dropout", "train_share", if (max(sim_log$ctrl_death) > 0) "death")) %>%
  mutate(measure = recode(measure, fall1 = "≥ 1 Sturz (%)", dropout = "Dropout (%)", death = "Tod (%)", train_share = "Training (%)"),
         quartile = factor(as.integer(quartile))) %>%
  ggplot(aes(quartile, mean, color = arm, group = arm)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_grid(rows = vars(measure), cols = vars(scenario), scales = "free_y") +
  labs(x = "Quartil der Baseline-Frailty (1 = am wenigsten frail)", y = NULL, color = NULL,
       title = "Ereignisse je Frailty-Quartil", subtitle = paste0("Mittelwert über ", n_runs, " Läufe")) +
  theme_ov

#--- Figure 3: cumulative dropout and death over the year -------------------------------
fig_time <- data_stats %>%
  filter(str_detect(measure, "^(drop|death)_w[0-9]+$")) %>%
  mutate(week  = as.integer(str_extract(measure, "[0-9]+$")),
         event = ifelse(str_detect(measure, "^drop"), "Dropout", "Tod")) %>%
  filter(event == "Dropout" | max(sim_log$ctrl_death) > 0) %>%
  ggplot(aes(week, 100 * mean, color = arm)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_grid(rows = vars(event), cols = vars(scenario)) +
  scale_x_continuous(breaks = c(13, 26, 39, 52)) +
  labs(x = "Woche", y = "Kumulierter Anteil (%)", color = NULL,
       title = "Dropout und Tod im Studienverlauf", subtitle = paste0("Mittelwert über ", n_runs, " Läufe")) +
  theme_ov

#--- Figure 4: variation of the realised values between the runs -------------------------
fig_runs <- runs_long %>%
  filter(measure %in% c("fall1", "dropout")) %>%
  left_join(scenarios %>% dplyr::select(scenario_id, scenario), by = "scenario_id") %>%
  mutate(measure = recode(measure, fall1 = "≥ 1 Sturz", dropout = "Dropout")) %>%
  ggplot(aes(100 * value, fill = arm)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.55) +
  facet_grid(rows = vars(measure), cols = vars(scenario), scales = "free_y") +
  labs(x = "Anteil je Lauf (%)", y = "Läufe", fill = NULL,
       title = "Schwankung zwischen den Läufen", subtitle = paste0(n_runs, " Läufe je Szenario")) +
  theme_ov

figures <- list(falls = fig_falls, frailty = fig_frailty, time = fig_time, runs = fig_runs)
for (nm in names(figures)) {
  print(figures[[nm]])
  ggsave(file.path(out_dir, sprintf("04_%s.png", nm)), figures[[nm]], width = width_ov, height = 5.5, dpi = 200, bg = "white")
}
