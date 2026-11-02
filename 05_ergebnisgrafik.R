###############################################################################
### Results plots and table
### - Plot 1: Bias of the effect estimate per model (informative scenario), one plot per method
### - Plot 2: Bias of the mean number of falls per arm (informative scenario)
### - Plot 3: Coverage of the 95% CI: robust SE vs bootstrap (informative scenario)
### - Plot 4: Time-varying treatment effect: log ratio of the falls per quarter (H1 scenarios)
### - Plot 5: Bias of the time-varying effect per quarter and over the whole year (H1 scenarios)
### - Table : performance measures per scenario and model (interactive HTML with column filters)
### Input: 00_aggregated_results.csv per block (from 03_multiple_analyse_cluster.R), see blocks / res_dirs
###############################################################################

library(tidyverse)
library(gt)
library(reactable)
library(htmlwidgets)

# blocks to evaluate: c("main"), c("supp") or c("main", "supp") (together)
blocks <- c("supp")

# one folder per block (cluster archive unpacked)
# scenario_id is given on the complete grid (scenarios.R) -> the blocks can simply be stacked
res_dirs <- c(main = "./results_cluster",
              supp = "./results_cluster/ergebnisse_supp",
              mort = "./results_cluster/ergebnisse_mort")[blocks]
res_files <- file.path(res_dirs, "00_aggregated_results.csv")
# all plots and tables in ./05_results, one subfolder per choice of blocks (main, supp, main_supp)
out_dir <- file.path("./05_results", paste(blocks, collapse = "_"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!all(file.exists(res_files))) stop("missing: ", paste(res_files[!file.exists(res_files)], collapse = ", "))

agg <- map(res_files, read_csv, show_col_types = FALSE) %>%
  list_rbind() %>%
  mutate(across(where(is.double), ~ ifelse(is.nan(.x), NA, .x)))   
if (anyDuplicated(agg[c("scenario_id", "Model")])) stop("scenario in more than one block folder")
n_runs <- max(agg$n_runs)
# realised dropout per arm only in results since the restructuring (older runs: NA)
for (v in c("mean_drop_ctrl", "mean_drop_int")) if (!v %in% names(agg)) agg[[v]] <- NA_real_

# model labels without method prefix (LWYY / GL / MCF / MCF-LWYY / MCF-GL)
# 8 / 9 / 10: sym + synthetic proxy P(t) = FrailtyT(t) + noise (sd 2 / 1 / 0.5, GlobalVars.R: proxy_sd),
# only 2_Informative without mortality, LWYY and MCF only; order weak -> strong proxy -> oracle
model_levels <- c("0: Full data (oracle M0)", "1: Naive (Unweighted)",
                  "3: IPCW sym", "4: IPCW sym_train",
                  "8: IPCW sym+P sd=2", "9: IPCW sym+P sd=1", "10: IPCW sym+P sd=0.5",
                  "5: IPCW oracle",
                  "6: IPCW sym trunc", "7: IPCW sym_train trunc")
model_labels <- c("Vollständige Daten", "Naiv (ungewichtet)",
                  "IPCW sym", "IPCW sym + Training",
                  "IPCW sym + P (SD 2)", "IPCW sym + P (SD 1)", "IPCW sym + P (SD 0,5)",
                  "IPCW Oracle",
                  "IPCW sym (trunkiert)", "IPCW sym + Training (trunkiert)")

#treatment effect (scenarios.R: effect_def)
effect_labels <- c(none = "H0 (β = 0)", low = "H1 niedrig (β = −0,15)", high = "H1 hoch (β = −0,4)")

# scenario names (scenarios.R: scenarios_def), used in plot 4 and in the compact table
scenario_labels <- c("1_NoFrailty"    = "Keine Frailty",
                     "2_Informative"  = "Informativ",
                     "3_NoInfDropout" = "Dropout nicht informativ",
                     "4_ObsDropout"   = "Dropout über Sturzhistorie",
                     "5_NoIndirect"   = "Ohne indirekten Effekt")

# highlight sym and sym + Training (truncated: lighter) and the P models (blue, darker = stronger proxy), rest in grey
model_colors <- c("Vollständige Daten"              = "#8a8984",
                  "Naiv (ungewichtet)"              = "#8a8984",
                  "IPCW sym"                        = "#009E73",
                  "IPCW sym + Training"             = "#CC79A7",
                  "IPCW sym + P (SD 2)"             = "#9ECAE1",
                  "IPCW sym + P (SD 1)"             = "#4292C6",
                  "IPCW sym + P (SD 0,5)"           = "#08519C",
                  "IPCW Oracle"                     = "#8a8984",
                  "IPCW sym (trunkiert)"            = "#7FCEB9",
                  "IPCW sym + Training (trunkiert)" = "#E6BCD3")

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
      dropout    = factor(ifelse(drop_beta_arm == 0, "β_C = 0", "β_C = −1"), levels = c("β_C = 0", "β_C = −1")),
      beta_label = factor(effect_labels[effect], levels = effect_labels),   # none / low / high
      mort_label = factor(mortality_level, levels = c("none", "high", "very_high"),
                          labels = c("keine Mortalität", "hohe Mortalität", "sehr hohe Mortalität"))
    )
}

theme_res <- list(
  theme_minimal(base_size = 11),
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),                        # row labels (β_C, mortality) readable
        panel.spacing = unit(1, "lines")),
  labs(caption = "β_C: Armeffekt auf den Dropout (log-OR je Woche); Kontrolle immer ca. 20 %, −1: Trainingsarm bricht seltener ab")
)

# plots 1-3: one facet column per effect, mortality (if more than one level) and β_C in the rows
# main: 2 columns x 2 rows, main + supp: 3 columns x 4 rows
facet_inf <- function(dat) {
  rows <- if (n_distinct(dat$mort_label) > 1) vars(mort_label, dropout) else vars(dropout)
  facet_grid(rows = rows, cols = vars(beta_label))
}
dims_inf   <- agg %>% filter(scenario_name == "2_Informative") %>%
  summarise(n_eff = n_distinct(effect), n_mort = n_distinct(mortality_level))
width_inf  <- 3.5 + 4 * dims_inf$n_eff
height_inf <- 1.5 + 2.2 * 2 * dims_inf$n_mort

# save via a temporary file: if the PNG is open in another program (viewer, OneDrive sync),
# ggsave would fail and leave its device open, so the following plots would land in that file
save_plot <- function(file, plot, width, height) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, plot, width = width, height = height, dpi = 200, bg = "white")
  if (!file.copy(tmp, file, overwrite = TRUE)) warning("could not write ", file, " (open in another program?)")
  unlink(tmp)
}


#--- Plot 1: bias of the effect estimate, one plot per method ---------------------------
# each method against its own true value (True_Beta_LWYY, True_Beta_GL, True_LogRateRatio)
plot_bias <- function(method) {
  dat <- agg %>%
    filter(scenario_name == "2_Informative", Method == method) %>%
    prep_plot() %>%
    mutate(lower = mean_bias - 1.96 * mcse_bias,                   # Monte Carlo interval (95%)
           upper = mean_bias + 1.96 * mcse_bias)

  ggplot(dat, aes(x = mean_bias, y = model, color = model)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.6, orientation = "y") +
    geom_point(size = 2.5) +
    facet_inf(dat) +
    scale_color_manual(values = model_colors, guide = "none") +
    labs(title    = paste("Bias im informativen Szenario:", method_labels[[method]]),
         subtitle = paste0("Punkt = mittlerer Bias, Linie = ±1,96 Monte-Carlo-SE (", n_runs, " Läufe pro Szenario)"),
         x = "Bias (Schätzer − wahrer Wert, log-Skala)", y = NULL) +
    theme_res
}

for (m in c("LWYY", "GL", "MCF")) {
  p <- plot_bias(m)
  print(p)
  save_plot(file.path(out_dir, sprintf("05_bias_%s.png", m)), p, width = width_inf, height = height_inf)
}


#--- Plot 2: bias of the mean number of falls per arm -----------------------------------
# shows which arm drives the bias of the ratio (e.g. sym_train corrects only the training arm)
# relative bias in % of the true mean: bias of the log ratio ≈ rel. bias training − rel. bias control,
# so the grey line between the two arms is about the bias of the ratio (x 100), a short line = small ratio bias
# MCF: nonparametric, real world with death; MCF-LWYY: compared with the limit of LWYY on complete data (MCF_LWYY_*)
plot_falls <- function(method) {
  dat <- agg %>%
    filter(scenario_name == "2_Informative", Method == method) %>%
    prep_plot() %>%
    mutate(Training  = 100 * bias_mcf_intervention / true_mcf_intervention,
           Kontrolle = 100 * bias_mcf_control / true_mcf_control)

  ggplot(dat, aes(y = model)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_segment(aes(x = Kontrolle, xend = Training, yend = model), color = "grey80", linewidth = 1.2) +
    geom_point(aes(x = Kontrolle, color = "Kontrolle"), size = 2.8) +
    geom_point(aes(x = Training,  color = "Training"),  size = 2.8) +
    facet_inf(dat) +
    scale_color_manual(values = c("Kontrolle" = "#D55E00", "Training" = "#0072B2")) +
    labs(title    = paste0("Relativer Bias der mittleren Sturzzahl bis Woche 52 je Arm\n", method_labels[[method]]),
         subtitle = paste0("Graue Linie = Abstand der Arme ≈ Bias des log-Verhältnisses × 100; ",
                           n_runs, " Läufe pro Szenario"),
         x = "Relativer Bias (% des wahren Mittelwerts)", y = NULL, color = "Arm") +
    theme_res +
    theme(legend.position = "top",
          strip.text.y = element_text(angle = 0),                        # β_C readable
          panel.grid.major.y = element_line(color = "grey95"))
}

for (m in c("MCF", "LWYY", "GL")) {
  p <- plot_falls(m)
  print(p)
  save_plot(file.path(out_dir, sprintf("05_bias_falls_per_arm_%s.png", m)), p, width = width_inf, height = height_inf)
}


#--- Plot 3: coverage of the 95% CI --------------------------------------------------------
# robust SE: sandwich SE for LWYY / GL, robust SE from mets for the MCF (mcf_mets in functions.R);
# bootstrap (Wald with bootstrap SE, percentile) where available (06_bootstrap_cluster.R, only MCF)
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
    facet_inf(dat) +
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
  save_plot(file.path(out_dir, sprintf("05_coverage_%s.png", m)), p, width = width_inf, height = height_inf)
}


#--- Plot 4: time-varying treatment effect ----------------------------------------------------
# training lowers the frailty over time -> the effect grows during the year
# log ratio of the mean number of falls IN each quarter (MCF-Q1..Q4, nonparametric) against the true value per quarter
# dashed line: true log ratio over the whole year (52 weeks) = one constant effect, as LWYY / GL / MCF estimate it
# only H1 (under H0 there are no quarter rows); complete data left out (lies on the truth, only Monte Carlo noise)
quarter_labels <- c(Q1 = "Woche 1–13", Q2 = "14–26", Q3 = "27–39", Q4 = "40–52")
time_colors    <- c("Wahrer Wert"           = "black",
                    "Naiv (ungewichtet)"    = "#8a8984",
                    "IPCW sym"              = "#009E73",
                    "IPCW sym + Training"   = "#CC79A7",
                    "IPCW sym + P (SD 2)"   = "#9ECAE1",    # P models only in 2_Informative without mortality
                    "IPCW sym + P (SD 1)"   = "#4292C6",
                    "IPCW sym + P (SD 0,5)" = "#08519C",
                    "IPCW Oracle"           = "#E69F00")

time_dat <- agg %>%
  filter(str_detect(Method, "^MCF-Q"), effect != "none") %>%
  prep_plot() %>%
  mutate(model    = as.character(model),
         quarter  = factor(quarter_labels[sub("^MCF-", "", Method)], levels = quarter_labels),
         dropout  = factor(dropout, levels = c("β_C = 0", "β_C = −1")),
         szenario = factor(scenario_labels[scenario_name], levels = scenario_labels)) %>%
  filter(model != "Vollständige Daten")
# true value per quarter as an extra "model" -> line and legend entry in black
time_true <- time_dat %>%
  distinct(scenario_id, szenario, effect, mortality_level, beta_label, mort_label, dropout, quarter,
           mean_estimate = true_beta) %>%
  mutate(model = "Wahrer Wert", mcse_bias = NA_real_)
# true value over the whole year (Method "MCF", 52 weeks)
year_true <- agg %>%
  filter(Method == "MCF", scenario_id %in% time_dat$scenario_id) %>%
  distinct(scenario_id, true_year = true_beta) %>%
  inner_join(time_true %>% distinct(scenario_id, szenario, effect, mortality_level, dropout), by = "scenario_id")

# one plot per effect x mortality (main + supp: up to 4), columns = scenarios
plot_time_effect <- function(eff, mort) {
  dat <- bind_rows(time_dat, time_true) %>%
    filter(effect == eff, mortality_level == mort) %>%
    mutate(model = factor(model, levels = names(time_colors)))
  yr <- year_true %>% filter(effect == eff, mortality_level == mort)
  pd <- position_dodge(width = 0.4)

  ggplot(dat, aes(x = quarter, y = mean_estimate, color = model, group = model)) +
    geom_hline(data = yr, aes(yintercept = true_year), linetype = "dashed", color = "grey40") +
    geom_line(position = pd, linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_estimate - 1.96 * mcse_bias, ymax = mean_estimate + 1.96 * mcse_bias),
                  width = 0, linewidth = 0.5, position = pd, na.rm = TRUE) +
    geom_point(aes(size = model == "Wahrer Wert"), position = pd) +
    facet_grid(rows = vars(dropout), cols = vars(szenario)) +
    scale_color_manual(values = time_colors) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.8), guide = "none") +
    labs(title    = paste0("Zeitlich veränderlicher Behandlungseffekt: Sturzzahl je Quartal\n",
                           dat$beta_label[1], ", ", dat$mort_label[1]),
         subtitle = paste0("Schwarz = wahrer Effekt im Quartal, gestrichelt = wahrer Effekt über das ganze Jahr; ",
                           "Linie = ±1,96 Monte-Carlo-SE; ", n_runs, " Läufe pro Szenario"),
         x = "Quartal", y = "log-Verhältnis der mittleren Sturzzahl im Quartal (Training / Kontrolle)", color = NULL) +
    theme_res +
    theme(legend.position = "top",
          strip.text.y = element_text(angle = 0),
          panel.grid.major.y = element_line(color = "grey92"))
}

time_combos <- time_dat %>% distinct(effect, mortality_level)    # empty if the results contain no H1 scenario
for (i in seq_len(nrow(time_combos))) {
  eff  <- time_combos$effect[i]
  mort <- time_combos$mortality_level[i]
  n_cols_time <- time_dat %>% filter(effect == eff, mortality_level == mort) %>% distinct(szenario) %>% nrow()
  p <- plot_time_effect(eff, mort)
  print(p)
  save_plot(file.path(out_dir, sprintf("05_time_effect_MCF_%s_%s.png", eff, mort)), p,
            width = 3 + 3.5 * n_cols_time, height = 6.5)
}


#--- Plot 5: bias of the time-varying effect per quarter --------------------------------------
# bias (estimate − true value) of the log ratio per quarter (MCF-Q1..Q4) and over the whole year (MCF, 52 weeks)
# shows whether the bias grows during the year (more dropout -> more selection)
# complete data in light grey = reference for the Monte Carlo noise; same facets and files as plot 4
bias_time_colors <- c("Vollständige Daten" = "#c9c8c4", time_colors[names(time_colors) != "Wahrer Wert"])
year_label       <- "Jahr (1–52)"

bias_time_dat <- agg %>%
  filter(Method %in% c("MCF", paste0("MCF-Q", 1:4)), effect != "none") %>%
  prep_plot() %>%
  mutate(model    = factor(as.character(model), levels = names(bias_time_colors)),
         period   = ifelse(Method == "MCF", year_label, quarter_labels[sub("^MCF-", "", Method)]),
         period   = factor(period, levels = c(quarter_labels, year_label)),
         szenario = factor(scenario_labels[scenario_name], levels = scenario_labels))

plot_time_bias <- function(eff, mort) {
  dat <- bias_time_dat %>% filter(effect == eff, mortality_level == mort)
  pd  <- position_dodge(width = 0.5)

  ggplot(dat, aes(x = period, y = mean_bias, color = model, group = model)) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 4.5, linetype = "dotted", color = "grey60") +        # quarters | whole year
    geom_line(data = filter(dat, period != year_label), position = pd, linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_bias - 1.96 * mcse_bias, ymax = mean_bias + 1.96 * mcse_bias),
                  width = 0, linewidth = 0.5, position = pd) +
    geom_point(size = 1.8, position = pd) +
    facet_grid(rows = vars(dropout), cols = vars(szenario)) +
    scale_color_manual(values = bias_time_colors) +
    labs(title    = paste0("Bias des zeitlich veränderlichen Behandlungseffekts (MCF) je Quartal und über das ganze Jahr\n",
                           dat$beta_label[1], ", ", dat$mort_label[1]),
         subtitle = paste0("Bias = Schätzer − wahrer Wert im Zeitraum (log-Verhältnis der mittleren Sturzzahl); ",
                           "Linie = ±1,96 Monte-Carlo-SE; ", n_runs, " Läufe pro Szenario"),
         x = "Zeitraum (Wochen)", y = "Bias (log-Skala)", color = NULL) +
    theme_res +
    theme(legend.position = "top",
          panel.grid.major.y = element_line(color = "grey92"))
}

for (i in seq_len(nrow(time_combos))) {                          # same effect x mortality combinations as plot 4
  eff  <- time_combos$effect[i]
  mort <- time_combos$mortality_level[i]
  n_cols_time <- bias_time_dat %>% filter(effect == eff, mortality_level == mort) %>% distinct(szenario) %>% nrow()
  p <- plot_time_bias(eff, mort)
  print(p)
  save_plot(file.path(out_dir, sprintf("05_time_bias_MCF_%s_%s.png", eff, mort)), p,
            width = 3 + 3.5 * n_cols_time, height = 6.5)
}


#--- Interactive table: performance measures per scenario and model -----------------------
# reactable instead of one long gt table: filter field above each column (e.g. "S10" in Nr., "LWYY" in Methode)
tab_dat <- agg %>%
  arrange(scenario_id) %>%
  transmute(Nr          = sprintf("S%02d", scenario_id),
            Szenario    = scenario_labels[scenario_name],
            Effekt      = effect,
            beta_C      = drop_beta_arm,
            Drop_K      = mean_drop_ctrl,                 # realised dropout per arm
            Drop_T      = mean_drop_int,
            Mortalitaet = mortality_level,
            Methode     = Method,
            Modell      = sub("^\\S+ ", "", Model),
            n_runs, true_beta, mean_estimate, mean_bias, mcse_bias, empirical_SE, mean_model_SE, SE_ratio, RMSE,
            coverage, power, mean_mcf_control, true_mcf_control, mean_mcf_intervention, true_mcf_intervention) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))

results_tbl <- reactable(
  tab_dat,
  filterable = TRUE, searchable = TRUE, defaultPageSize = 35, compact = TRUE, striped = TRUE,
  defaultColDef = colDef(na = "–"),
  columns = list(
    beta_C                = colDef(name = "β_C"),
    Drop_K                = colDef(name = "Dropout Kontrolle", format = colFormat(percent = TRUE, digits = 1)),
    Drop_T                = colDef(name = "Dropout Training", format = colFormat(percent = TRUE, digits = 1)),
    Mortalitaet           = colDef(name = "Mortalität"),
    Modell                = colDef(minWidth = 160),
    n_runs                = colDef(name = "Läufe"),
    true_beta             = colDef(name = "β wahr"),
    mean_estimate         = colDef(name = "β̂ (Mittel)"),
    mean_bias             = colDef(name = "Bias"),
    mcse_bias             = colDef(name = "MCSE"),
    empirical_SE          = colDef(name = "Emp. SE"),
    mean_model_SE         = colDef(name = "Robuster SE"),
    SE_ratio              = colDef(name = "SE-Ratio"),
    coverage              = colDef(name = "Überdeckung", format = colFormat(percent = TRUE, digits = 1)),
    power                 = colDef(name = "Ablehnung", format = colFormat(percent = TRUE, digits = 1)),
    mean_mcf_control      = colDef(name = "MCF Kontrolle"),
    true_mcf_control      = colDef(name = "MCF Kontrolle wahr"),
    mean_mcf_intervention = colDef(name = "MCF Training"),
    true_mcf_intervention = colDef(name = "MCF Training wahr")
  )
)

# without pandoc the JavaScript files are saved in the folder "lib" next to the HTML file
saveWidget(results_tbl, file.path(normalizePath(out_dir), "00_aggregated_results.html"),
           selfcontained = rmarkdown::pandoc_available(), libdir = "lib", title = "Simulationsergebnisse")


#--- Compact table: true value and mean estimate of the main methods per scenario ----------
# one row per scenario, mean effect estimate for LWYY / GL / MCF (52 weeks, no quarters)
# columns: true value, complete data (reference for the Monte Carlo noise), naive and the three IPCW versions
# cells are coloured by |bias| (hidden bias columns)
bias_pal <- function(x)                                                       # colour by |bias|
  scales::col_numeric(c("#FFFFFF", "#F4A582", "#CA0020"), domain = c(0, 0.12))(pmin(abs(x), 0.12))
model_short <- c("0" = "Vollst.", "1" = "Naiv", "3" = "sym", "4" = "sym+Train",
                 "8" = "sym+P 2", "9" = "sym+P 1", "10" = "sym+P 0,5", "5" = "Oracle",
                 "6" = "sym trunk.", "7" = "sym+Train trunk.")

key_long <- agg %>%
  filter(Method %in% c("LWYY", "GL", "MCF")) %>%
  mutate(model = model_short[str_extract(Model, "\\d+(?=:)")])      # \\d+: two digits for model 10

# true value per method (True_Beta_LWYY, True_Beta_GL, True_LogRateRatio) as extra column "wahr"
key_true <- key_long %>%
  group_by(scenario_id, scenario_name, effect, true_beta_arm_dgp, drop_beta_arm, mortality_level, Method) %>%   
  summarise(mean_estimate = first(true_beta), .groups = "drop") %>%
  mutate(model = "wahr")

# columns "LWYY_Naiv" (mean estimate) and "LWYY_Naiv_bias" (only for the colour)
key_dat <- bind_rows(key_long, key_true) %>%
  mutate(col = paste(Method, model, sep = "_")) %>%
  dplyr::select(scenario_id, scenario_name, effect, true_beta_arm_dgp, drop_beta_arm, mortality_level, col,   
                est = mean_estimate, bias = mean_bias) %>%
  pivot_wider(names_from = col, values_from = c(est, bias), names_glue = "{col}_{.value}") %>%
  rename_with(~ sub("_est$", "", .x)) %>%
  mutate(hypothesis = effect_labels[effect],                          
         scenario   = sprintf("S%02d", scenario_id),
         szenario   = scenario_labels[scenario_name],
         dropout    = ifelse(drop_beta_arm == 0, "0", "−1"),
         mort       = recode(mortality_level, none = "keine", high = "hoch")) %>%
  arrange(scenario_id)

key_cols  <- as.vector(outer(c("LWYY", "GL", "MCF"), c("wahr", model_short), paste, sep = "_") %>% t())   # per method: true value, then estimates
key_cols  <- intersect(key_cols, names(key_dat))                  # results without truncated models (older runs)
true_cols <- grep("_wahr$", key_cols, value = TRUE)
est_cols  <- setdiff(key_cols, true_cols)
mcse_med  <- median(agg$mcse_bias[agg$Method %in% c("LWYY", "GL", "MCF")], na.rm = TRUE)

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
                         ". Vollst. = vollständige Daten ohne Dropout (Referenz); sym+Train = Training nur im Trainingsarm im Dropout-Modell; trunk. = Gewichte je Arm am 99. Perzentil gekürzt.",
                         " sym+P x = sym + synthetischer Proxy P(t) = Frailty(t) + Fehler mit SD x (nur informatives Szenario ohne Mortalität, nur LWYY und MCF).")) %>%
  tab_source_note("β_C: Armeffekt auf den Dropout (log-OR je Woche); Kontrolle immer ca. 20 %, −1: Trainingsarm bricht seltener ab. Farbe: |Bias| = |Schätzer − wahrer Wert| (weiß = 0, rot ≥ 0,12).") %>%
  tab_options(row_group.font.weight = "bold", table.font.size = px(12), data_row.padding = px(2))

gtsave(key_gt, file.path(out_dir, "00_key_results.html"))
