###############################################################################
### Results plots and table
### - Plot 1: Bias of the effect estimate per model (informative scenario), one plot per method
### - Plot 2: Bias of the mean number of falls per arm (informative scenario)
### - Plot 3: Coverage of the 95% CI: robust SE vs bootstrap (informative scenario)
### - Table : performance measures per scenario and model
### Input: ./results/00_aggregated_results.csv (from 03_mehrfach_analyse.R)
###############################################################################

library(tidyverse)
library(gt)

agg <- read_csv("./results/00_aggregated_results.csv", show_col_types = FALSE) %>%
  mutate(across(where(is.double), ~ ifelse(is.nan(.x), NA, .x)))    # NaN (e.g. coverage of MCF rows) -> NA
n_runs <- max(agg$n_runs)
if (!dir.exists("./res_tables")) dir.create("./res_tables")

# model labels without method prefix (LWYY / GL / MCF / MCF-LWYY / MCF-GL)
model_levels <- c("0: Full data (oracle M0)", "1: Naive (Unweighted)",
                  "3: IPCW sym", "4: IPCW sym_train", "5: IPCW oracle")
model_labels <- c("Vollständige Daten", "Naiv (ungewichtet)",
                  "IPCW sym", "IPCW sym + Training", "IPCW Oracle")

# highlight the two proxy models, rest in grey
model_colors <- c("Vollständige Daten"  = "#8a8984",
                  "Naiv (ungewichtet)"  = "#8a8984",
                  "IPCW sym"            = "#009E73",
                  "IPCW sym + Training" = "#CC79A7",
                  "IPCW Oracle"         = "#8a8984")

# Method column from 03 (Method = text before the first blank of Model)
method_labels <- c("LWYY"     = "LWYY (Tod = Zensierung), log-Rate-Ratio",
                   "GL"       = "Ghosh-Lin, log-Verhältnis der mittleren Sturzzahl",
                   "MCF"      = "Mittlere Sturzzahl nichtparametrisch (weighted_mcf), log-Verhältnis")

# common labels for the facets of the informative scenario
prep_plot <- function(dat) {
  dat %>%
    mutate(
      model      = factor(sub("^\\S+ ", "", Model), levels = model_levels, labels = model_labels),
      model      = fct_rev(model),                                  # first model on top
      # beta_C = extra arm effect on dropout; with frailty the dropout differs between the arms also for beta_C = 0
      dropout    = ifelse(drop_beta_arm == 0, "β_C = 0", "β_C = −1"),
      beta_label = factor(ifelse(true_beta_arm_dgp == 0, "H0 (β = 0)", "H1 (β = −0,4)"),
                          levels = c("H0 (β = 0)", "H1 (β = −0,4)")),
      mort_label = factor(mortality_level, levels = c("none", "high"),
                          labels = c("keine Mortalität", "hohe Mortalität"))
    )
}

theme_res <- list(
  theme_minimal(base_size = 11),
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines")),
  labs(caption = "β_C: zusätzlicher Armeffekt auf den Dropout (−1: Kontrolle bricht häufiger ab)")
)


#--- Plot 1: bias of the effect estimate, one plot per method ---------------------------
# each method against its own true value (True_Beta_LWYY, True_Beta_GL, True_LogRateRatio)
plot_bias <- function(method) {
  agg %>%
    filter(scenario_name == "2_Informative", Method == method) %>%
    prep_plot() %>%
    mutate(lower = mean_bias - 1.96 * mcse_bias,                   # Monte Carlo interval (95%)
           upper = mean_bias + 1.96 * mcse_bias) %>%
    ggplot(aes(x = mean_bias, y = model, color = model)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.6, orientation = "y") +
    geom_point(size = 2.5) +
    facet_grid(dropout ~ beta_label + mort_label) +
    scale_color_manual(values = model_colors, guide = "none") +
    labs(title    = paste("Bias im informativen Szenario:", method_labels[[method]]),
         subtitle = paste0("Punkt = mittlerer Bias, Linie = ±1,96 Monte-Carlo-SE (", n_runs, " Läufe pro Szenario)"),
         x = "Bias (Schätzer − wahrer Wert, log-Skala)", y = NULL) +
    theme_res
}

for (m in c("LWYY", "GL", "MCF")) {
  p <- plot_bias(m)
  print(p)
  ggsave(sprintf("./results/05_bias_%s.png", m), p, width = 16, height = 5.5, dpi = 200, bg = "white")
}


#--- Plot 2: bias of the mean number of falls per arm -----------------------------------
# shows which arm drives the bias of the ratio (e.g. sym_train corrects only the training arm)
# MCF: nonparametric, real world with death; MCF-LWYY: compared with the truth without death
plot_falls <- function(method) {
  agg %>%
    filter(scenario_name == "2_Informative", Method == method) %>%
    prep_plot() %>%
    dplyr::select(model, dropout, beta_label, mort_label,
                  Training = bias_mcf_intervention, Kontrolle = bias_mcf_control) %>%
    pivot_longer(c(Training, Kontrolle), names_to = "arm", values_to = "bias") %>%
    ggplot(aes(x = bias, y = model, color = arm)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(dropout ~ beta_label + mort_label) +
    scale_color_manual(values = c("Kontrolle" = "#D55E00", "Training" = "#0072B2")) +
    labs(title    = paste("Bias der mittleren Sturzzahl bis Woche 52 je Arm:", method_labels[[method]]),
         subtitle = paste0(n_runs, " Läufe pro Szenario"),
         x = "Bias (Stürze pro Person)", y = NULL, color = "Arm") +
    theme_res +
    theme(legend.position = "top")
}

for (m in c("MCF", "LWYY", "GL")) {
  p <- plot_falls(m)
  print(p)
  ggsave(sprintf("./results/05_bias_falls_per_arm_%s.png", m), p, width = 16, height = 5.5, dpi = 200, bg = "white")
}


#--- Plot 3: coverage of the 95% CI --------------------------------------------------------
# robust (sandwich) SE only for LWYY / GL; bootstrap (Wald with bootstrap SE, percentile) where available
# error bars: ±1.96 Monte Carlo SE of the coverage
plot_coverage <- function(method) {
  dat <- agg %>%
    filter(scenario_name == "2_Informative", Method == method) %>%
    prep_plot() %>%
    transmute(model, dropout, beta_label, mort_label, n_runs, n_boot_runs,
              `Robuster SE`          = coverage,
              `Bootstrap-SE (Wald)`  = boot_coverage_wald,
              `Bootstrap Perzentil`  = boot_coverage_perc) %>%
    pivot_longer(c(`Robuster SE`, `Bootstrap-SE (Wald)`, `Bootstrap Perzentil`),
                 names_to = "ci", values_to = "cov") %>%
    filter(!is.na(cov)) %>%
    mutate(n    = ifelse(ci == "Robuster SE", n_runs, n_boot_runs),
           mcse = sqrt(cov * (1 - cov) / n))
  if (nrow(dat) == 0) return(NULL)                                 # e.g. MCF without bootstrap

  ggplot(dat, aes(x = cov, y = model, color = ci)) +
    geom_vline(xintercept = 0.95, color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = cov - 1.96 * mcse, xmax = cov + 1.96 * mcse), width = 0, linewidth = 0.5,
                  orientation = "y", position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    facet_grid(dropout ~ beta_label + mort_label) +
    scale_color_manual(values = c("Robuster SE" = "#8a8984", "Bootstrap-SE (Wald)" = "#0072B2",
                                  "Bootstrap Perzentil" = "#E69F00")) +
    labs(title    = paste("Überdeckung des 95%-KI im informativen Szenario:", method_labels[[method]]),
         subtitle = "Linie = ±1,96 Monte-Carlo-SE; Bootstrap nur in den Läufen mit Bootstrap (n_boot_runs)",
         x = "Überdeckung", y = NULL, color = NULL) +
    theme_res +
    theme(legend.position = "top")
}

for (m in c("LWYY", "GL", "MCF")) {
  p <- plot_coverage(m)
  if (is.null(p)) next
  print(p)
  ggsave(sprintf("./results/05_coverage_%s.png", m), p, width = 16, height = 5.5, dpi = 200, bg = "white")
}


#--- Table: performance measures per scenario and model -----------------------------------
results_gt <- agg %>%
  mutate(scenario_label = sprintf("S%02d – %s, β = %s, drop_beta_arm = %s, Mortalität: %s", scenario_id, scenario_name,
                                  true_beta_arm_dgp, drop_beta_arm, mortality_level)) %>%
  dplyr::select(scenario_label, Model, n_runs, true_beta, mean_estimate,
                mean_bias, mcse_bias, empirical_SE, mean_model_SE, SE_ratio, RMSE, coverage, power,
                n_boot_runs, mean_boot_SE, boot_SE_ratio, boot_coverage_wald, boot_coverage_perc,
                mean_mcf_control, true_mcf_control, mean_mcf_intervention, true_mcf_intervention) %>%
  gt(groupname_col = "scenario_label") %>%
  fmt_number(columns = c(true_beta, mean_estimate, mean_bias, mcse_bias, empirical_SE,
                         mean_model_SE, SE_ratio, RMSE, mean_boot_SE, boot_SE_ratio), decimals = 3) %>%
  fmt_percent(columns = c(coverage, power, boot_coverage_wald, boot_coverage_perc), decimals = 1) %>%
  fmt_number(columns = c(mean_mcf_control, true_mcf_control,
                         mean_mcf_intervention, true_mcf_intervention), decimals = 3) %>%
  sub_missing(missing_text = "–") %>%                              # NA: no SE / CI for MCF rows, no bootstrap
  cols_merge(columns = c(mean_bias, mcse_bias), pattern = "{1} ({2})") %>%
  cols_merge(columns = c(mean_mcf_control, true_mcf_control), pattern = "{1} ({2})") %>%
  cols_merge(columns = c(mean_mcf_intervention, true_mcf_intervention), pattern = "{1} ({2})") %>%
  cols_label(
    Model                 = "Modell",
    n_runs                = "Läufe",
    true_beta             = "β wahr",
    mean_estimate         = "β̂ (Mittel)",
    mean_bias             = "Bias (MCSE)",
    empirical_SE          = "Emp. SE",
    mean_model_SE         = "Robuster SE",
    SE_ratio              = "SE-Ratio",
    RMSE                  = "RMSE",
    coverage              = "Überdeckung",
    power                 = "Ablehnung",
    n_boot_runs           = "Läufe",
    mean_boot_SE          = "SE",
    boot_SE_ratio         = "SE-Ratio",
    boot_coverage_wald    = "Überd. Wald",
    boot_coverage_perc    = "Überd. Perzentil",
    mean_mcf_control      = "Kontrolle",
    mean_mcf_intervention = "Training"
  ) %>%
  tab_spanner(label = "Bootstrap",
              columns = c(n_boot_runs, mean_boot_SE, boot_SE_ratio, boot_coverage_wald, boot_coverage_perc)) %>%
  tab_spanner(label = "Mittlere Sturzzahl bis Woche 52: geschätzt (wahr)",
              columns = c(mean_mcf_control, mean_mcf_intervention)) %>%
  tab_header(title = "Simulationsergebnisse je Szenario und Modell") %>%
  tab_footnote(
    footnote = "Ablehnung: Power unter H1, Fehler 1. Art unter H0 (zweiseitiger Wald-Test, robuster SE).",
    locations = cells_column_labels(columns = power)
  ) %>%
  tab_footnote(
    footnote = "Nur für MCF-Zeilen; MCF-LWYY wird mit der Wahrheit ohne Tod verglichen.",
    locations = cells_column_spanners(spanners = "Mittlere Sturzzahl bis Woche 52: geschätzt (wahr)")
  ) %>%
  tab_options(row_group.font.weight = "bold")

gtsave(results_gt, "./res_tables/00_aggregated_results.html")


#--- Compact table: true value and mean estimate of the main methods per scenario ----------
# one row per scenario, mean effect estimate for LWYY / GL / MCF (52 weeks, no quarters)
# columns: true value, complete data (reference for the Monte Carlo noise), naive and the three IPCW versions
# cells are coloured by |bias| (hidden bias columns)
scenario_labels <- c("1_NoFrailty"    = "Keine Frailty",
                     "2_Informative"  = "Informativ",
                     "3_NoInfDropout" = "Dropout nicht informativ",
                     "4_ObsDropout"   = "Dropout über L1",
                     "5_NoIndirect"   = "Ohne indirekten Effekt")
model_short <- c("0" = "Vollst.", "1" = "Naiv", "3" = "sym", "4" = "sym+Train", "5" = "Oracle")

key_long <- agg %>%
  filter(Method %in% c("LWYY", "GL", "MCF")) %>%
  mutate(model = model_short[str_extract(Model, "\\d(?=:)")])

# true value per method (True_Beta_LWYY, True_Beta_GL, True_LogRateRatio) as extra column "wahr"
key_true <- key_long %>%
  group_by(scenario_id, scenario_name, true_beta_arm_dgp, drop_beta_arm, mortality_level, Method) %>%
  summarise(mean_estimate = first(true_beta), .groups = "drop") %>%
  mutate(model = "wahr")

# columns "LWYY_Naiv" (mean estimate) and "LWYY_Naiv_bias" (only for the colour)
key_dat <- bind_rows(key_long, key_true) %>%
  mutate(col = paste(Method, model, sep = "_")) %>%
  dplyr::select(scenario_id, scenario_name, true_beta_arm_dgp, drop_beta_arm, mortality_level, col,
                est = mean_estimate, bias = mean_bias) %>%
  pivot_wider(names_from = col, values_from = c(est, bias), names_glue = "{col}_{.value}") %>%
  rename_with(~ sub("_est$", "", .x)) %>%
  mutate(hypothesis = ifelse(true_beta_arm_dgp == 0, "H0: β = 0", "H1: β = −0,4"),
         scenario   = sprintf("S%02d", scenario_id),
         szenario   = scenario_labels[scenario_name],
         dropout    = ifelse(drop_beta_arm == 0, "0", "−1"),
         mort       = recode(mortality_level, none = "keine", high = "hoch")) %>%
  arrange(scenario_id)

key_cols  <- as.vector(outer(c("LWYY", "GL", "MCF"), c("wahr", model_short), paste, sep = "_") %>% t())   # per method: true value, then estimates
true_cols <- grep("_wahr$", key_cols, value = TRUE)
est_cols  <- setdiff(key_cols, true_cols)
mcse_med  <- median(agg$mcse_bias[agg$Method %in% c("LWYY", "GL", "MCF")], na.rm = TRUE)
bias_pal  <- function(x)                                                       # colour by |bias|
  scales::col_numeric(c("#FFFFFF", "#F4A582", "#CA0020"), domain = c(0, 0.12))(pmin(abs(x), 0.12))

key_gt <- key_dat %>%
  dplyr::select(hypothesis, scenario, szenario, dropout, mort, all_of(key_cols), all_of(paste0(est_cols, "_bias"))) %>%
  gt(groupname_col = "hypothesis") %>%
  fmt_number(columns = all_of(key_cols), decimals = 3)
for (cl in est_cols) {                                                         # estimate cell gets the colour of its bias
  key_gt <- key_gt %>% data_color(columns = all_of(paste0(cl, "_bias")), target_columns = all_of(cl), fn = bias_pal)
}
key_gt <- key_gt %>%
  cols_hide(columns = ends_with("_bias")) %>%
  tab_style(style = list(cell_fill(color = "#EEEEEE"), cell_text(weight = "bold")),   # true values grey
            locations = cells_body(columns = all_of(true_cols))) %>%
  tab_spanner(label = "LWYY (log-Rate-Ratio)",    columns = starts_with("LWYY_")) %>%
  tab_spanner(label = "Ghosh-Lin (log-Ratio μ)",  columns = starts_with("GL_")) %>%
  tab_spanner(label = "MCF nichtparametrisch",    columns = starts_with("MCF_")) %>%
  cols_label(.list = c(list(scenario = "", szenario = "Szenario", dropout = "β_C", mort = "Mortalität"),
                       setNames(as.list(sub("^[^_]+_", "", key_cols)), key_cols))) %>%
  tab_header(title    = "Wahrer Wert und mittlerer Effektschätzer je Szenario",
             subtitle = paste0(n_runs, " Läufe pro Szenario; grau = wahrer Wert, übrige Spalten = mittlerer Schätzer (log-Skala)")) %>%
  tab_source_note(paste0("Typischer Monte-Carlo-SE des Bias (Median): ", sprintf("%.3f", mcse_med),
                         ". Vollst. = vollständige Daten ohne Dropout (Referenz); sym+Train = Training nur im Trainingsarm im Dropout-Modell.")) %>%
  tab_source_note("β_C: zusätzlicher Armeffekt auf den Dropout (−1: Kontrolle bricht häufiger ab). Farbe: |Bias| = |Schätzer − wahrer Wert| (weiß = 0, rot ≥ 0,12).") %>%
  tab_options(row_group.font.weight = "bold", table.font.size = px(12), data_row.padding = px(2))

gtsave(key_gt, "./res_tables/00_key_results.html")
