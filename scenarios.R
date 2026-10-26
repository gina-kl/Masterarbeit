############################################################################
### Scenarios: complete design and blocks ##################################
############################################################################
# needed in 01 and 03 and cluster_jobs/n_scenarios.R (number of array tasks for submit_all.sh)
# needs: dplyr, tidyr, tibble and GlobalVars.R (scenario_blocks, causal_vars)
# scenario_id and causal_id are given on the COMPLETE grid, the blocks are filtered afterwards
# -> a scenario keeps its id, file names and seeds, whichever blocks are run


# Scenarios
#  - 1_NoFrailty: frailty without any effect (proportional rates, non-informative dropout)
#  - 2_Informative: frailty -> falls, death, dropout, training; training -> frailty
#  - 3_NoInfDropout: as 2, but dropout does not depend on frailty
#  - 4_ObsDropout: as 2, but dropout depends only on the observed falls up to t-1, log(1 + falls) (log_falls_lag),
#                  NEU: MAR given the observed history -> IPCW sym is correctly specified (positive control)
#  - 5_NoIndirect: as 2, but training does not reduce frailty (no indirect effect via frailty)
scenarios_def <- tribble(
  ~scenario_name,    ~eff_frailty_y, ~eff_frailty_death, ~eff_frailty_drop, ~eff_L1_drop, ~eff_fall_drop, ~eff_frailty_train, ~train_frailty_eff,
  "1_NoFrailty",     0.0,            0.0,               0.0,               0.0,          0.0,            0.0,                0.0,
  "2_Informative",   1.2,            1.2,               0.8,               0.0,          0.0,           -1.5,                0.02,
  "3_NoInfDropout",  1.2,            1.2,               0.0,               0.0,          0.0,           -1.5,                0.02,
  "4_ObsDropout",    1.2,            1.2,               0.0,               0.0,          1.5,           -1.5,                0.02,   # NEU: eff_fall_drop instead of eff_L1_drop = 1
  "5_NoIndirect",    1.2,            1.2,               0.8,               0.0,          0.0,           -1.5,                0.0
)

# treatment effect: direct effect (true_beta_arm) and maximal frailty reduction by training (max_frailty_red)
# none: H0 (no direct effect, training does not change frailty either)
# low:  realistic, both parts scaled down by the same factor, marginal log RR about -0.26 in 2_Informative
# high: former effect, marginal log RR about -0.58 in 2_Informative
effect_def <- tribble(
  ~effect, ~true_beta_arm, ~max_frailty_red,
  "none",  0.0,            0.0,
  "low",   -0.15,          0.18,       #log RR about -0.26 (RR 0.77)
  "high",  -0.4,           0.5
)

scenario_grid <- expand_grid(
  n_pts = 1800,
  effect = effect_def$effect,            # none (H0), low, high 
  drop_beta_arm = c(0.0, -1.0),          # 0 = same Dropout-Rate, -1.0 = control group drops out more often
  mortality_level = c("none", "high")    # "none" = no death, "high" = non-negligible mortality
) %>%
  left_join(effect_def, by = "effect") %>%
  cross_join(scenarios_def) %>%
  filter(scenario_name %in% c("1_NoFrailty", "2_Informative") |  # scenario 1 and 2 complete, others only under H1
           effect != "none") %>%
  mutate(train_frailty_eff = ifelse(effect == "none", 0, train_frailty_eff),   # strict H0: no effect via training either
         scenario_id = row_number())

# causal scenario: same causal parameters -> same complete data (same seed) and same truth
# numbering in order of first appearance on the complete grid
causal_ids <- scenario_grid %>%
  distinct(across(all_of(causal_vars))) %>%
  mutate(causal_id = row_number())
scenario_grid <- scenario_grid %>% left_join(causal_ids, by = causal_vars)

# blocks
#   main: core question, no mortality: 2_Informative x {H0, low} x drop_beta_arm,
#         3_NoInfDropout x low x drop_beta_arm (6 scenarios)
#   supp: all others (1_NoFrailty, 4_ObsDropout (validation), 5_NoIndirect, high effect, mortality "high")
scenario_grid <- scenario_grid %>%
  mutate(block = ifelse(mortality_level == "none" &
                          ((scenario_name == "2_Informative" & effect %in% c("none", "low")) |
                             (scenario_name == "3_NoInfDropout" & effect == "low")),
                        "main", "supp"))

base_scenarios <- scenario_grid %>% filter(block %in% scenario_blocks)
if (nrow(base_scenarios) == 0) stop("no scenarios in block(s) ", paste(scenario_blocks, collapse = ", "))
