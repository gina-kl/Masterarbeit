# =============================================================================
# time_check.R – Laufzeit- und Speichertest, rechnet die Werte für die Job-Skripte hoch
#
# Ausführen: im Projektordner (neben 01_generate_weekly_cluster.R, functions.R,
# GlobalVars.R, cluster_setup.R) in RStudio "Source" klicken.
# NEU: run_on_cluster wird hier auf FALSE gesetzt (cluster_setup.R erkennt den Cluster sonst
# automatisch, z. B. in einer salloc-Sitzung) – kein manuelles Umstellen mehr nötig.
# Dauer: je nach Rechner ca. 2-10 Minuten. Es wird nichts gespeichert.
# =============================================================================

run_on_cluster <- FALSE   # NEU: überschreibt die automatische Erkennung in cluster_setup.R

# ---- Geplante Größen auf dem Cluster (bei Bedarf anpassen) -------------------
n_sim_cluster      <- 1000   # Läufe pro Szenario
n_true_cluster     <- 1e6    # Patienten pro Arm für die wahren Werte
cores_data         <- 16     # --cpus-per-task in job_data.sh
cores_true         <- 48     # --cpus-per-task in job_true.sh
cores_analyse      <- 32     # --cpus-per-task in job_analyse.sh
puffer             <- 1.5    # Sicherheitsaufschlag (Cluster-Kerne können langsamer sein, Schwankungen)
n_reps             <- 3      # wie oft jede Messung wiederholt wird (Mittelwert)
worker_overhead_mb <- 400    # NEU: jeder PSOCK-Worker lädt die Pakete selbst (wird unten nicht gemessen)

# ---- Setup aus 01_generate laden, OHNE die Simulation zu starten -------------
# (führt nur den Teil bis zur Szenario-Definition aus: Pakete, Funktionen, COEF, base_scenarios)
gen_lines <- readLines("01_generate_weekly_cluster.R", encoding = "UTF-8")
i_setup   <- grep("^run_true <- ", gen_lines)[1]
eval(parse(text = gen_lines[1:(i_setup - 1)]), envir = globalenv())
if (run_on_cluster) stop("Für den Zeittest in cluster_setup.R run_on_cluster <- FALSE setzen")

# true_block() und unique_causal_scenarios definieren (stehen im Teil "true")
i1 <- grep("unique_causal_scenarios <- base_scenarios", gen_lines)[1]
i2 <- grep('cat\\("Set up a parallel backend for true values', gen_lines)[1]
eval(parse(text = gen_lines[i1:(i2 - 1)]), envir = globalenv())

library(ipw)       # für die Analyse (compute_weights)
library(splines)

# Hilfsfunktionen
timeit <- function(expr_fun, reps = n_reps) {
  gc(reset = TRUE)
  t <- system.time(for (r in seq_len(reps)) res <- expr_fun(r))["elapsed"] / reps
  list(sec = unname(t), mem_mb = sum(gc()[, 6]), res = res)   # gc()[, 6] = maximal genutzter Speicher (MB)
}
fmt_time <- function(sec) {
  m <- ceiling(sec / 60)
  sprintf("%02d:%02d:00", m %/% 60, m %% 60)
}

# Aufwendigstes Szenario: informativ, mit Mortalität, H1
p_heavy <- base_scenarios %>%
  filter(scenario_name == "2_Informative", mortality_level == "high", true_beta_arm == -0.4) %>%
  slice(1)

sim_one <- function(seed) simulate_scenario(
  seed = seed, n_pts = p_heavy$n_pts, true_beta_arm = p_heavy$true_beta_arm,
  drop_beta_arm = p_heavy$drop_beta_arm, Y_intercept = p_heavy$Y_intercept,
  drop_intercept = p_heavy$drop_intercept, death_intercept = p_heavy$death_intercept,
  eff_frailty_y = p_heavy$eff_frailty_y, eff_frailty_death = p_heavy$eff_frailty_death,
  eff_frailty_drop = p_heavy$eff_frailty_drop, eff_L1_drop = p_heavy$eff_L1_drop,
  eff_frailty_train = p_heavy$eff_frailty_train, train_intercept = p_heavy$train_intercept,
  train_frailty_eff = p_heavy$train_frailty_eff)

# ---- 1) Ein Datensatz --------------------------------------------------------
cat("\n[1/4] Datensatz erzeugen ...\n")
m_data <- timeit(function(r) sim_one(seed = 900000 + r))
Odat_test <- m_data$res
cat(sprintf("      %.1f s pro Datensatz, max. Speicher %.0f MB\n", m_data$sec, m_data$mem_mb))

# ---- 2) Ein Block der wahren Werte -------------------------------------------
cat("[2/4] Ein Block der wahren Werte (", n_chunk_true, " Patienten pro Arm) ...\n", sep = "")
p_true <- unique_causal_scenarios %>%
  filter(mortality_level == "high", eff_frailty_y > 0, true_beta_arm == -0.4) %>% slice(1)
m_true <- timeit(function(r) true_block(p_true, block = r, n_chunk = n_chunk_true, seed = 1e9 + r), reps = 1)
cat(sprintf("      %.1f s pro Block, max. Speicher %.0f MB\n", m_true$sec, m_true$mem_mb))

# ---- 3) Eine Analyse (wie analyse_one_run, ohne Bootstrap) --------------------
cat("[3/4] Eine Analyse (Gewichte, 10 Modelle, MCF) ...\n")
Odat_all <- Odat_test %>%
  mutate(log_base_falls = log1p(base_risk), tstop = t,
         fall_rate_lag  = cum_falls_lag / pmax(t - 1, 1)) %>%
  arrange(ID, t)
m_ana <- timeit(function(r) {
  Odat_w <- compute_weights(Odat_all %>% filter(observed))
  fits   <- fit_models(Odat_all, Odat_w)
  mcf_estimate(Odat_all, Odat_w, fits)
})
cat(sprintf("      %.1f s pro Lauf, max. Speicher %.0f MB\n", m_ana$sec, m_ana$mem_mb))

# ---- 4) Eine Bootstrap-Stichprobe --------------------------------------------
cat("[4/4] Eine Bootstrap-Stichprobe ...\n")
set.seed(1)
m_boot <- timeit(function(r) tryCatch(boot_estimate(boot_sample(Odat_all)),   # wie in run_bootstrap:
                                      error = function(e) NULL))              # Fehler stoppen den Test nicht
cat(sprintf("      %.1f s pro Bootstrap-Stichprobe, max. Speicher %.0f MB\n", m_boot$sec, m_boot$mem_mb))

# ---- Hochrechnung --------------------------------------------------------------
n_causal     <- nrow(unique_causal_scenarios)
n_blocks_cl  <- n_true_cluster / n_chunk_true
n_boot_runs  <- min(boot_max_run, n_sim_cluster)

sec_data     <- n_sim_cluster * m_data$sec / cores_data * puffer
sec_true     <- ceiling(n_causal * n_blocks_cl / cores_true) * m_true$sec * puffer
sec_ana      <- n_sim_cluster * m_ana$sec / cores_analyse * puffer
sec_ana_boot <- (n_sim_cluster * m_ana$sec + n_boot_runs * n_boot * m_boot$sec) / cores_analyse * puffer

# NEU: + Overhead pro Worker (gemessen wird nur der Speicher im Hauptprozess)
mem_per_core <- function(mb) paste0(ceiling((mb + worker_overhead_mb) * 1.5 / 500) * 500, "mb")   # +50 %, auf 500 MB gerundet

cat("\n=================== EMPFEHLUNG FÜR DIE JOB-SKRIPTE ===================\n")
cat(sprintf("(Hochrechnung: n_sim = %d, n_true = %g, n_boot = %d, Puffer x%.1f)\n\n",
            n_sim_cluster, n_true_cluster, n_boot, puffer))
cat(sprintf("job_data.sh    --cpus-per-task=%d  --time=%s  --mem-per-cpu=%s\n",
            cores_data, fmt_time(sec_data), mem_per_core(m_data$mem_mb)))
cat(sprintf("job_true.sh    --cpus-per-task=%d  --time=%s  --mem-per-cpu=%s   (%d Blöcke)\n",
            cores_true, fmt_time(sec_true), mem_per_core(m_true$mem_mb), n_causal * n_blocks_cl))
cat(sprintf("job_analyse.sh --cpus-per-task=%d  --time=%s  --mem-per-cpu=%s\n",
            cores_analyse, fmt_time(max(sec_ana, sec_ana_boot)), mem_per_core(max(m_ana$mem_mb, m_boot$mem_mb))))
cat(sprintf("               (Szenarien ohne Bootstrap: %s, mit Bootstrap: %s)\n",
            fmt_time(sec_ana), fmt_time(sec_ana_boot)))

if (max(sec_data, sec_true, sec_ana_boot) > 72 * 3600) {
  cat("\nACHTUNG: mindestens ein Job läge über dem Maximum von 72 h ->",
      "mehr Kerne anfordern (bis 96 in Queue cpu) oder die Arbeit weiter aufteilen.\n")
}
cat("=======================================================================\n")
cat("Hinweis: lokal gemessen. Cluster-Kerne können schneller oder langsamer sein –\n",
    "am genauesten ist derselbe Test in einer salloc-Sitzung (dev_cpu) auf dem Cluster.\n")
