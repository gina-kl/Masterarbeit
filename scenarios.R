############################################################################
### Scenarios: complete design and blocks ##################################
############################################################################
# needed in 01 and 03 and cluster_jobs/n_scenarios.R (number of array tasks for submit_all.sh)
# needs: dplyr, tidyr, tibble and GlobalVars.R (scenario_blocks, causal_vars)
# scenario_id and causal_id are given on the complete grid, the blocks are filtered afterwards
# -> a scenario keeps its id, file names and seeds, whichever blocks are run


# Scenarios
#  - 1_NoFrailty: frailty without any effect (proportional rates, non-informative dropout)
#  - 2_Informative: frailty -> falls, death, dropout, training; training -> frailty
#  - 3_NoInfDropout: as 2, but dropout does not depend on frailty
#  - 4_ObsDropout: as 2, but dropout depends only on the observed falls up to t-1,
#                  MAR given the observed history -> IPCW sym is correctly specified (positive control)
#  - 5_NoIndirect: as 2, but training does not reduce frailty (no indirect effect via frailty)

scenarios_def <- tribble(
  ~scenario_name,    ~eff_frailty_y, ~eff_frailty_death, ~eff_frailty_drop, ~eff_fall_drop, ~eff_frailty_train, ~train_frailty_eff,
  "1_NoFrailty",     0.0,            0.0,               0.0,               0.0,            0.0,                0.0,
  "2_Informative",   1.2,            0.6,               0.8,               0.0,           -1.0,                0.02,
  "3_NoInfDropout",  1.2,            0.6,               0.0,               0.0,           -1.0,                0.02,
  "4_ObsDropout",    1.2,            0.6,               0.0,               1.5,           -1.0,                0.02,
  "5_NoIndirect",    1.2,            0.6,               0.8,               0.0,           -1.0,                0.0
)

# treatment effect: direct effect (true_beta_arm) and maximal frailty reduction by training (max_frailty_red)
# none: H0 (no direct effect, training does not change frailty either)
# low:  realistic, marginal log RR about -0.29 in 2_Informative, at least one fall: 35 % -> about 28 %
# high: former effect, stronger than in exercise trials
effect_def <- tribble(
  ~effect, ~true_beta_arm, ~max_frailty_red, 
  "none",  0.0,            0.0,
  "low",   -0.15,          0.18,       
  "high",  -0.4,           0.5
)

scenario_grid <- expand_grid(
  n_pts = 1800,
  effect = effect_def$effect,            # none (H0), low, high 
  drop_beta_arm = c(0.0, -1.0),          # arm effect on dropout (log-OR per week): control always ~20 %, -1.0 = training arm drops out less often
  mortality_level = c("none", "high")    # "none" = no death, "high" = non-negligible mortality
) %>%
  left_join(effect_def, by = "effect") %>%
  cross_join(scenarios_def) %>%
  filter(scenario_name %in% c("1_NoFrailty", "2_Informative") |  # scenario 1 and 2 complete, others only under H1
           effect != "none") %>%
  mutate(train_frailty_eff = ifelse(effect == "none", 0, train_frailty_eff),   # strict H0: no effect via training either
         scenario_id = row_number())

# add very high mortality for main scenarios
mort_rows <- scenario_grid %>%
  filter(mortality_level == "none" &    # each combination only once
           ((scenario_name == "2_Informative" & effect %in% c("none", "low")) |
              (scenario_name == "4_ObsDropout" & effect == "low"))) %>%
  mutate(mortality_level = "very_high", scenario_id = max(scenario_grid$scenario_id) + row_number()) # add at the end with new scenario ids

scenario_grid <- bind_rows(scenario_grid, mort_rows)

# sensitivity analyses, only beta_C = -1 (unbalanced dropout)
#  - 2_Informative_Y08: weaker effect of frailty on falls (H0 and low effect)
#  - 2_Informative_Tr0: training does not depend on frailty, train_rate_lag is no proxy (low effect)
# own names: the true values are matched by scenario name, mortality and beta
# appended at the end (scenario ids and streams are not changed)
sens_base <- scenario_grid %>%
  filter(scenario_name == "2_Informative", mortality_level == "none", drop_beta_arm == -1)
sens_rows <- bind_rows(
  sens_base %>% filter(effect %in% c("none", "low")) %>%
    mutate(scenario_name = "2_Informative_Y08", eff_frailty_y = 0.8),
  sens_base %>% filter(effect == "low") %>%
    mutate(scenario_name = "2_Informative_Tr0", eff_frailty_train = 0)) %>%
  mutate(scenario_id = max(scenario_grid$scenario_id) + row_number())

scenario_grid <- bind_rows(scenario_grid, sens_rows)

# different strength of informative dropout (block sens; only for beta_C = -1, no effect and low effect, no mortality)
# - 2_Informative_Drop04: eff_frailty_drop = 0.4
# - 2_Informative_Drop12: eff_frailty_drop = 1.2
# same causal scenarios as S6 / S20 -> paired comparison
drop_base <- scenario_grid %>%
  filter(scenario_name == "2_Informative", mortality_level == "none", drop_beta_arm == -1)
drop_rows <- bind_rows(
    drop_base %>% filter(effect %in% c("none", "low"))%>%
      mutate(scenario_name = "2_Informative_Drop04", eff_frailty_drop = 0.4),
  drop_base %>% filter(effect %in% c("none", "low"))%>%
    mutate(scenario_name = "2_Informative_Drop12", eff_frailty_drop = 1.2)) %>%
  mutate(scenario_id = max(scenario_grid$scenario_id) + row_number())  

scenario_grid <- bind_rows(scenario_grid, drop_rows)

# stronger effect of frailty on death (block mort; very high mortality, beta_C = -1, H0 and low effect)
#  - 2_Informative_Death12: eff_frailty_death = 1.2 (as strong as on the falls, instead of 0.6)
# eff_frailty_death is a causal variable -> new causal_id: new complete data and new true values
death_rows <- scenario_grid %>%
  filter(scenario_name == "2_Informative", mortality_level == "very_high", drop_beta_arm == -1) %>%   
  mutate(scenario_name = "2_Informative_Death12", eff_frailty_death = 1.2,
         scenario_id = max(scenario_grid$scenario_id) + row_number())                                  
scenario_grid <- bind_rows(scenario_grid, death_rows)

# causal scenario: same causal parameters -> same complete data (same seed) and same truth
# numbering in order of first appearance on the complete grid
causal_ids <- scenario_grid %>%
  distinct(across(all_of(causal_vars))) %>%
  mutate(causal_id = row_number())
scenario_grid <- scenario_grid %>% left_join(causal_ids, by = causal_vars)


# blocks
#   main: core question, no mortality: 2_Informative x {H0, low} x drop_beta_arm,
#         4_ObsDropout x low x drop_beta_arm (positive control) (6 scenarios)
#   supp: reference / all others (1_NoFrailty, 3_NoInfDropout (negative control), 5_NoIndirect,
#         high effect, mortality "high")
#   mort: very high mortality (15 %): the main scenarios and 2_Informative_Death12 (eff_frailty_death 1.2) (8 scenarios)
#   sens: sensitivity analyses of 2_Informative with beta_C = -1: Y08, Tr0 and
#         Drop04 / Drop12 (eff_frailty_drop 0.4 / 1.2) (7 scenarios)
scenario_grid <- scenario_grid %>%
  mutate(block = ifelse(mortality_level == "none" &
                          ((scenario_name == "2_Informative" & effect %in% c("none", "low")) |
                             (scenario_name == "4_ObsDropout" & effect == "low")),
                        "main", "supp"))
scenario_grid <- scenario_grid %>%
  mutate(block = ifelse(mortality_level == "very_high", "mort", block),          # also 2_Informative_Death12
         block = ifelse(scenario_name %in% c("2_Informative_Y08", "2_Informative_Tr0",
                                             "2_Informative_Drop04", "2_Informative_Drop12"), "sens", block))

base_scenarios <- scenario_grid %>% filter(block %in% scenario_blocks) # scenario_blocks = main /supp/mort/sens

if (nrow(base_scenarios) == 0) stop("no scenarios in block(s) ", paste(scenario_blocks, collapse = ", "))
