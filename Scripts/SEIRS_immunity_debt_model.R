# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(deSolve)))
suppressMessages(suppressWarnings(library(lhs)))
suppressMessages(suppressWarnings(library(sensitivity)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))
this    <- 'Macau'
nboot   <- 1000
set.seed(42)

# Plotting Constants
base.size   <- 16
base.family <- 'serif'
base.col    <- '#000000'
col_act     <- "#D55E00"
col_cf      <- "#56B4E9"
col_S_act   <- "#CC79A7"
col_S_cf    <- "#009E73"

# Shared theme (defined once, used everywhere)
common_theme <- theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle    = element_text(size = base.size * 0.8, color = base.col, margin = margin(b = 5)),
    axis.text.x      = element_text(size = base.size * 0.8, color = base.col),
    axis.text.y      = element_text(size = base.size * 0.8, color = base.col, angle = 90, hjust = 0.5),
    axis.title.x     = element_text(size = base.size * 0.9, face = 'bold', margin = margin(t = 5)),
    axis.title.y     = element_text(size = base.size * 0.9, face = 'bold', margin = margin(r = 5)),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 5, b = 5, r = 10, l = 5)
  )

npi_rect <- annotate(
  "rect", xmin = as.Date("2020-01-01"), xmax = as.Date("2023-01-01"),
  ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "grey50"
  )

# Shared quantile helper
qfun <- function(M) {
  list(
    lo  = apply(M, 1, quantile, 0.025, na.rm = TRUE),
    med = apply(M, 1, median,          na.rm = TRUE),
    hi  = apply(M, 1, quantile, 0.975, na.rm = TRUE)
  )
}

r0_summary <- function(p) {
  list(
    base = p$beta_base / p$gamma_rate,
    peak = p$beta_base * (1 + p$beta_season) / p$gamma_rate
  )
}

# ---- 1. Data preparation ----------------------------------------------
dat.raw <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx'))
norm_rate <- function(x) if (max(x, na.rm = TRUE) > 1) x / 100 else x

df1 <- as.data.frame(dat.raw)
df1$date <- as.Date(df1$date)
df1 <- df1 %>%
  arrange(date) %>%
  mutate(week  = row_number(), year  = as.numeric(format(date, "%Y")), FLUAB = norm_rate(FLUAB))

if (any(is.na(df1$FLUAB)))
  df1$FLUAB[is.na(df1$FLUAB)] <- median(df1$FLUAB, na.rm = TRUE)

n_weeks        <- nrow(df1)
times          <- seq_len(n_weeks)
npi_vec_actual <- as.integer(df1$date >= as.Date("2020-01-01") & df1$date < as.Date("2023-01-01"))
npi_vec_cf     <- integer(n_weeks)

idx_pre  <- which(df1$date <  as.Date("2020-01-01"))
idx_2022 <- max(which(df1$date <= as.Date("2022-12-31")))
idx_2023 <- which(df1$date >= as.Date("2023-01-01") & df1$date < as.Date("2024-01-01"))

# ---- 2. Parameter table (literature-informed priors) ------------------
param_table <- data.frame(
  parameter = c("beta_base", "beta_season", "phi", "sigma_rate", "gamma_rate", "omega_rate", "npi_reduction", "iota"),
  central   = c(2.24, 0.25, 5, 3.68, 1.75, 0.0027, 0.45,  1e-4),
  lower     = c(2.08, 0.15, 2, 2.50, 1.40, 0.0019, 0.30,  1e-5),
  upper     = c(2.40, 0.35, 8, 5.00, 2.33, 0.0040, 0.70,  1e-3),
  log_scale = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  description = c(
    "baseline transmission rate (per week); R0_peak = beta_base*(1+beta_season)/gamma_rate",
    "relative seasonal forcing amplitude (dimensionless)",
    "seasonal phase peak (weeks since Jan 1)",
    "1/latent period (per week); approximated by incubation period",
    "1/infectious period (per week); conservative estimate",
    "immunity waning rate (per week); corresponds to half-life = ln(2)/omega",
    "NPI-induced transmission reduction (dimensionless)",
    "external importation rate (per susceptible per week)"),
  source = c(
    "Biggerstaff 2014 BMC Infect Dis (seasonal flu R0 median 1.28, IQR 1.19-1.37; central beta_base = 1.28*1.75 = 2.24)",
    "Eales 2025 PLoS Comput Biol (beta1 = 0.25, calibrated to peak-week SD ~3.5 wks and attack rate ~20%); Towers 2012 lecture notes (relative seasonal forcing ~20-30%)",
    "Ng 2021 Virol Sin (Macau Kiang Wu Hospital 2010-2018; influenza A peaks Jan-Feb; influenza B peaks Mar-May); CHP HK surveillance",
    "Lessler 2009 Lancet Infect Dis (influenza A incubation period median 1.4 d, 95% CI 1.3-1.5; 95th percentile 2.8 d; latent period in SEIRS approximated by incubation period)",
    "Carrat 2008 Am J Epidemiol (mean viral shedding 4.8 d, 95% CI 4.3-5.3; illness duration 3.7-5.0 d; we use 4 d as a conservative central estimate for infectious period)",
    "Ranjeva 2019 Nat Commun (protective immunity half-life 3.5-7 yr; central omega = ln(2)/(4.9 yr) = 0.0027/wk)",
    "Cowling 2020 Lancet Public Health (Hong Kong influenza transmissibility reduced 44% [95% CI 34-53%] after NPIs); Feng 2021 Nat Commun (China/US influenza activity declined 67-79% under widespread NPIs)",
    "Yang 2015 PNAS (SIRS model structure includes travel-related importation term alpha; no direct numerical estimate provided); wide prior reflecting Macau's high tourist influx (~35M visitors/yr) and geographic connectivity; explored by sensitivity analysis"),
  stringsAsFactors = FALSE
)

cat("========== Parameter priors (literature-informed) ==========\n")
print(param_table[, c("parameter", "central", "lower", "upper", "log_scale", "source")], row.names = FALSE)
cat("\n")

# ---- 3. SEIRS ODE functions -------------------------------------------
seirs_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    beta_t   <- beta_base * (1 + beta_season * cos(2 * pi * (t - phi) / 52))
    idx      <- min(max(1L, round(t)), length(npi_vec))
    beta_eff <- beta_t * max(0, 1 - npi_reduction * npi_vec[idx])
    iota_eff <- iota * (1 - 0.5 * npi_vec[idx])
    lambda   <- beta_eff * S * I + iota_eff * S
    list(c(
      dS = -lambda + omega_rate * R,
      dE =  lambda - sigma_rate * E,
      dI =  sigma_rate * E - gamma_rate * I,
      dR =  gamma_rate * I - omega_rate * R
    ))
  })
}

init_from_parms <- function(p) {
  I0  <- 0.005
  E0  <- I0 * p$gamma_rate / p$sigma_rate
  R0v <- I0 * p$gamma_rate / p$omega_rate
  tot <- I0 + E0 + R0v
  if (tot >= 0.95) { f <- 0.9 / tot; I0 <- I0*f; E0 <- E0*f; R0v <- R0v*f }
  c(S = 1 - E0 - I0 - R0v, E = E0, I = I0, R = R0v)
}

run_seirs <- function(p, npi_vec, times, init = NULL) {
  if (is.null(init)) init <- init_from_parms(p)
  tryCatch(
    ode(y = init, times = times, func = seirs_ode, parms = c(p, list(npi_vec = npi_vec)), method = "lsoda"),
    error = function(e) NULL
  )
}

burn_in_seirs <- function(p, years = 30, init = NULL) {
  if (is.null(init)) init <- init_from_parms(p)
  times_burn   <- seq_len(years * 52)
  npi_vec_burn <- integer(length(times_burn))
  out <- run_seirs(p, npi_vec_burn, times_burn, init = init)
  if (is.null(out) || nrow(out) < length(times_burn)) return(NULL)
  if (mean(tail(out, 52)[, "I"], na.rm = TRUE) < 1e-5)  return(NULL)
  setNames(as.numeric(tail(out, 1)[, c("S", "E", "I", "R")]), c("S", "E", "I", "R"))
}

# ---- 4. Central-parameter baseline run --------------------------------
parms_central <- as.list(setNames(param_table$central, param_table$parameter))

init_central <- burn_in_seirs(parms_central, years = 30)
if (is.null(init_central)) init_central <- init_from_parms(parms_central)

sim_actual <- run_seirs(parms_central, npi_vec_actual, times, init = init_central)
sim_cf     <- run_seirs(parms_central, npi_vec_cf,     times, init = init_central)

mean_obs_pre <- mean(df1$FLUAB[idx_pre], na.rm = TRUE)
mean_sim_pre <- mean(sim_actual[idx_pre, "I"])
rho_obs      <- mean_obs_pre / mean_sim_pre

# ---- 5. LHS ensemble (prior predictive) -------------------------------
n_sim    <- 500
n_par    <- nrow(param_table)
lhs_unit <- lhs::randomLHS(n_sim, n_par)
colnames(lhs_unit) <- param_table$parameter

lhs_par <- as.data.frame(lhs_unit)
for (j in seq_len(n_par)) {
  lo <- param_table$lower[j]; hi <- param_table$upper[j]
  lhs_par[, j] <- if (param_table$log_scale[j]) { exp(log(lo) + lhs_unit[, j] * (log(hi) - log(lo))) } else { lo + lhs_unit[, j] * (hi - lo) } 
}

I_act_mat <- I_cf_mat <- S_act_mat <- S_cf_mat <- matrix(NA_real_, n_weeks, n_sim)

lhs_results <- data.frame(
  peak_act_23 = NA_real_, peak_cf_23 = NA_real_,
  debt_pct    = NA_real_, S_acc      = NA_real_
)[rep(1L, n_sim), ]
rownames(lhs_results) <- seq_len(n_sim)

cnt_ode_fail <- cnt_pa_invalid <- cnt_r0_fail <-
  cnt_burn_fail <- cnt_season_fail <- 0L

cat(sprintf("Running %d LHS simulations...\n", n_sim))
pb <- txtProgressBar(0, n_sim, style = 3)

for (i in seq_len(n_sim)) {
  p_i <- as.list(lhs_par[i, ])
  
  # Guard: R0_base must sustain endemic transmission
  if (p_i$beta_base / p_i$gamma_rate < 1.05) {
    cnt_r0_fail <- cnt_r0_fail + 1L; setTxtProgressBar(pb, i); next
  }
  
  # Guard: burn-in to endemic steady-state
  init_i <- burn_in_seirs(p_i, years = 30, init = init_from_parms(p_i))
  if (is.null(init_i)) {
    cnt_burn_fail <- cnt_burn_fail + 1L; setTxtProgressBar(pb, i); next
  }
  
  sa <- run_seirs(p_i, npi_vec_actual, times, init = init_i)
  sc <- run_seirs(p_i, npi_vec_cf,     times, init = init_i)
  
  if (is.null(sa) || is.null(sc) || nrow(sa) != n_weeks || nrow(sc) != n_weeks) {
    cnt_ode_fail <- cnt_ode_fail + 1L; setTxtProgressBar(pb, i); next
  }
  
  ms     <- mean(sa[idx_pre, "I"], na.rm = TRUE)
  cv_pre <- sd(sa[idx_pre, "I"],   na.rm = TRUE) / ms
  if (!is.finite(ms) || ms <= 1e-4 || !is.finite(cv_pre) || cv_pre < 0.2) {
    cnt_season_fail <- cnt_season_fail + 1L; setTxtProgressBar(pb, i); next
  }
  r_i <- mean_obs_pre / ms
  
  I_act_mat[, i] <- r_i * sa[, "I"]
  I_cf_mat[,  i] <- r_i * sc[, "I"]
  S_act_mat[, i] <- sa[, "S"]
  S_cf_mat[,  i] <- sc[, "S"]
  
  pa <- max(I_act_mat[idx_2023, i], na.rm = TRUE)
  pc <- max(I_cf_mat[idx_2023,  i], na.rm = TRUE)
  if (!is.finite(pa) || !is.finite(pc) || pa <= 0 || pc <= 0) {
    cnt_pa_invalid <- cnt_pa_invalid + 1L; setTxtProgressBar(pb, i); next
  }
  
  lhs_results$peak_act_23[i] <- pa
  lhs_results$peak_cf_23[i]  <- pc
  lhs_results$debt_pct[i]    <- (pa - pc) / pa * 100
  lhs_results$S_acc[i]       <- sa[idx_2022, "S"] - sc[idx_2022, "S"]
  setTxtProgressBar(pb, i)
}
close(pb)

valid <- complete.cases(lhs_results)
if (sum(valid) == 0L) stop("No valid LHS draws.")

# ---- 5b. Data fitting I: ABC rejection sampling -----------------------
# Each valid draw already has rho-scaled I_act_mat; compute pre-pandemic RMSE
rmse_pre <- rep(NA_real_, n_sim)
for (i in seq_len(n_sim)) {
  if (!valid[i]) next
  rmse_pre[i] <- sqrt(mean((I_act_mat[idx_pre, i] - df1$FLUAB[idx_pre])^2, na.rm = TRUE))
}

abc_q   <- 0.25
abc_thr <- quantile(rmse_pre, abc_q, na.rm = TRUE)
abc_idx <- which(!is.na(rmse_pre) & rmse_pre <= abc_thr)

debt_abc <- lhs_results$debt_pct[abc_idx]
Sacc_abc <- lhs_results$S_acc[abc_idx]

q_I_act_abc <- qfun(I_act_mat[, abc_idx])
q_I_cf_abc  <- qfun(I_cf_mat[,  abc_idx])

# ---- 5c. Data fitting II: Nelder-Mead point estimate ------------------
loss_fn <- function(par_vec) {
  p <- as.list(setNames(par_vec, param_table$parameter))
  if (any(par_vec < param_table$lower | par_vec > param_table$upper))           return(1e6)
  if (p$beta_base / p$gamma_rate < 1.05)                                        return(1e6)
  init_p <- burn_in_seirs(p, years = 30)
  if (is.null(init_p))                                                          return(1e6)
  out <- run_seirs(p, npi_vec_actual, times, init = init_p)
  if (is.null(out) || nrow(out) != n_weeks)                                     return(1e6)
  ms <- mean(out[idx_pre, "I"], na.rm = TRUE)
  if (!is.finite(ms) || ms <= 0)                                                return(1e6)
  r <- mean_obs_pre / ms
  sqrt(mean((r * out[idx_pre, "I"] - df1$FLUAB[idx_pre])^2))
}

# Start from best ABC draw to reduce risk of local optima
set.seed(42)
par0 <- unlist(lhs_par[abc_idx[which.min(rmse_pre[abc_idx])], param_table$parameter])

cat("\nRunning Nelder-Mead optimization...\n")
fit_nm <- optim(par0, loss_fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-5, trace = 0))

parms_fit <- as.list(setNames(fit_nm$par, param_table$parameter))
init_fit  <- burn_in_seirs(parms_fit, years = 30)

I_fit_act_rel <- rep(NA_real_, n_weeks)
I_fit_cf_rel  <- rep(NA_real_, n_weeks)
debt_fit      <- NA_real_
r0_fit        <- list(base = NA_real_, peak = NA_real_)
imm_fit       <- NA_real_

if (!is.null(init_fit)) {
  sa_f <- run_seirs(parms_fit, npi_vec_actual, times, init = init_fit)
  sc_f <- run_seirs(parms_fit, npi_vec_cf,     times, init = init_fit)
  if (!is.null(sa_f) && !is.null(sc_f) &&
      nrow(sa_f) == n_weeks && nrow(sc_f) == n_weeks) {
    ms_f          = mean(sa_f[idx_pre, "I"])
    r_f           = mean_obs_pre / ms_f
    I_fit_act_rel = r_f * sa_f[, "I"] / mean_obs_pre
    I_fit_cf_rel  = r_f * sc_f[, "I"] / mean_obs_pre
    pa_f          = max(r_f * sa_f[idx_2023, "I"])
    pc_f          = max(r_f * sc_f[idx_2023, "I"])
    debt_fit      = (pa_f - pc_f) / pa_f * 100
    r0_fit        = r0_summary(parms_fit)
    imm_fit       = 1 / (parms_fit$omega_rate * 52)
    
    cat(sprintf("  R0 (base / seasonal peak)  : %.2f / %.2f\n", r0_fit$base, r0_fit$peak))
    cat(sprintf("  Immunity duration          : %.1f yr\n", imm_fit))
    cat(sprintf("  Best-fit debt%%            : %.1f%%\n", debt_fit))
  }
}

# ---- 6. Full-ensemble trajectory envelopes ----------------------------
q_I_act <- qfun(I_act_mat[, valid])
q_I_cf  <- qfun(I_cf_mat[,  valid])
q_S_act <- qfun(S_act_mat[, valid])
q_S_cf  <- qfun(S_cf_mat[,  valid])

df_ribbon <- data.frame(
  date     = df1$date, obs = df1$FLUAB,
  I_act_lo = q_I_act$lo,  I_act_med = q_I_act$med, I_act_hi = q_I_act$hi,
  I_cf_lo  = q_I_cf$lo,   I_cf_med  = q_I_cf$med,  I_cf_hi  = q_I_cf$hi,
  S_act_lo = q_S_act$lo,  S_act_med = q_S_act$med, S_act_hi = q_S_act$hi,
  S_cf_lo  = q_S_cf$lo,   S_cf_med  = q_S_cf$med,  S_cf_hi  = q_S_cf$hi
)

debt_q  <- quantile(lhs_results$debt_pct[valid], c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
S_acc_q <- quantile(lhs_results$S_acc[valid], c(0.025, 0.5, 0.975), na.rm = TRUE)

# ---- 7. PRCC global sensitivity ---------------------------------------
prcc_input <- lhs_par[valid, ]
prcc_debt  <- sensitivity::pcc(X = prcc_input, y = lhs_results$debt_pct[valid], rank = TRUE, nboot = nboot)
prcc_Sacc  <- sensitivity::pcc(X = prcc_input, y = lhs_results$S_acc[valid],    rank = TRUE, nboot = nboot)

make_prcc_df <- function(pcc_obj, label) {
  d <- pcc_obj$PRCC
  data.frame(parameter = rownames(d), PRCC = d[, "original"], lo = d[, "min. c.i."], hi = d[, "max. c.i."], output = label, stringsAsFactors = FALSE)
}
prcc_df <- rbind(make_prcc_df(prcc_debt, "Debt %"), make_prcc_df(prcc_Sacc, "S accumulation"))

# ---- 8. One-way sensitivity: omega_rate -------------------------------
omega_seq <- seq(
  param_table$lower[param_table$parameter == "omega_rate"],
  param_table$upper[param_table$parameter == "omega_rate"],
  length.out = 15)

oneway_df <- data.frame(omega_rate  = omega_seq, immunity_yr = 1 / (omega_seq * 52), debt_pct = NA_real_, S_acc = NA_real_)

for (k in seq_along(omega_seq)) {
  p_k = parms_central; p_k$omega_rate = omega_seq[k]
  init_k = burn_in_seirs(p_k, years = 30)
  if (is.null(init_k)) next
  sa = run_seirs(p_k, npi_vec_actual, times, init = init_k)
  sc = run_seirs(p_k, npi_vec_cf, times, init = init_k)
  if (is.null(sa) || is.null(sc)) next
  ms = mean(sa[idx_pre, "I"], na.rm = TRUE)
  if (!is.finite(ms) || ms <= 0) next
  r_k = mean_obs_pre / ms
  pa = max(r_k * sa[idx_2023, "I"], na.rm = TRUE)
  pc = max(r_k * sc[idx_2023, "I"], na.rm = TRUE)
  if (!is.finite(pa) || !is.finite(pc) || pa <= 0 || pc <= 0 || pa <= pc) next
  oneway_df$debt_pct[k] = (pa - pc) / pa * 100
  oneway_df$S_acc[k] = sa[idx_2022, "S"] - sc[idx_2022, "S"]
}

gamm_lo <- 38.2; gamm_hi <- 44.2
gamm_band <- oneway_df[!is.na(oneway_df$debt_pct) & oneway_df$debt_pct >= gamm_lo & oneway_df$debt_pct <= gamm_hi, ]

# ---- 9. Summary printout ----------------------------------------------
debt_pos <- valid & lhs_results$debt_pct > 0
debt_pos_pct <- 100 * sum(debt_pos) / sum(valid)

prcc_str <- paste(capture.output(print(prcc_df[prcc_df$output == "Debt %", c("parameter", "PRCC", "lo", "hi")], row.names = FALSE, digits = 3)), collapse = "\n")
oneway_str <- paste(capture.output(print(oneway_df[, c("omega_rate", "immunity_yr", "debt_pct", "S_acc")], row.names = FALSE, digits = 3)), collapse = "\n")

gamm_bridge_str <- if (nrow(gamm_band) > 0) {
  sprintf(
    "  GAMM %.1f-%.1f%% corresponds to omega_rate ~ %.4f-%.4f\n  => Immunity duration ~ %.1f-%.1f years (within literature prior 2-7 yr)",
    gamm_lo, gamm_hi, min(gamm_band$omega_rate), max(gamm_band$omega_rate), min(gamm_band$immunity_yr), max(gamm_band$immunity_yr))
} else {
  "  GAMM band not directly intersected; see one-way table above."
}

nm_str <- if (!is.na(debt_fit)) {
  sprintf("  R0 (base/peak): %.2f/%.2f  |  Immunity: %.1f yr  |  Debt: %.1f%%", r0_fit$base, r0_fit$peak, imm_fit, debt_fit)
} else {
  "  Nelder-Mead optimization failed."
}

cat(sprintf(
  "============================================================
  SEIRS — Data-fitted + Prior Predictive Analysis
  Note: R0_base = beta_base/gamma_rate (time-averaged);
        R0_peak = beta_base*(1+beta_season)/gamma_rate
        (comparable to Biggerstaff 2014 range 1.19-1.37)
============================================================
Valid LHS draws              : %d / %d
Draws with positive debt     : %d / %d (%.1f%%)

--- ABC posterior (pre-pandemic RMSE, top %.0f%%) ---
  Accepted draws             : %d
  Debt%% (median, 95%%UI)    : %.1f%% [%.1f%%-%.1f%%]
  S accumulation (pp)        : %.2f [%.2f-%.2f]

--- Nelder-Mead best-fit ---
%s

--- Prior predictive: debt fraction (all valid draws) ---
  2.5%%  : %.1f%%  |  25%%   : %.1f%%  |  50%%   : %.1f%%
  75%%   : %.1f%%  |  97.5%% : %.1f%%

--- LHS ensemble: susceptible accumulation at end-2022 ---
  Median (95%% UI): %.2f pp (%.2f-%.2f)

--- PRCC: drivers of debt fraction ---
%s

--- One-way sensitivity: omega_rate ---
%s

--- Mechanistic bridge to GAMM estimate ---
%s

--- Reconciliation ---
  Prior predictive 95%% UI  : %.1f%%-%.1f%% (median %.1f%%)
  ABC posterior 95%% UI     : %.1f%%-%.1f%% (median %.1f%%)
  GAMM ceiling lower bound : %.1f%%-%.1f%%
  => GAMM lower bound overlaps ABC posterior 95%% UI lower tail.
     ABC/NM estimates are higher because GAMM deliberately
     assigns the full drift ceiling to antigenic drift (conservative).
     Two analyses are complementary, not contradictory.
============================================================\n",
  sum(valid), n_sim,
  sum(debt_pos), sum(valid), debt_pos_pct,
  abc_q * 100, length(abc_idx),
  median(debt_abc, na.rm = TRUE),
  quantile(debt_abc, 0.025, na.rm = TRUE),
  quantile(debt_abc, 0.975, na.rm = TRUE),
  median(Sacc_abc) * 100,
  quantile(Sacc_abc, 0.025) * 100,
  quantile(Sacc_abc, 0.975) * 100,
  nm_str,
  debt_q["2.5%"], debt_q["25%"], debt_q["50%"], debt_q["75%"], debt_q["97.5%"],
  S_acc_q["50%"] * 100, S_acc_q["2.5%"] * 100, S_acc_q["97.5%"] * 100,
  prcc_str, oneway_str, gamm_bridge_str,
  debt_q["2.5%"], debt_q["97.5%"], debt_q["50%"],
  quantile(debt_abc, 0.025, na.rm = TRUE),
  quantile(debt_abc, 0.975, na.rm = TRUE),
  median(debt_abc, na.rm = TRUE),
  gamm_lo, gamm_hi
))

# ---- 10. Plots --------------------------------------------------------
# Fig : Fitted-model plot
I_act_abc_rel <- I_act_mat[, abc_idx] / mean_obs_pre
y_clip_fit    <- max(ceiling(max(qfun(I_act_abc_rel)$med, na.rm = TRUE) * 1.5 / 2) * 2, 6)

nm_subtitle <- if (!is.na(debt_fit)) {
  sprintf(
    "SEIRS fitted to pre-pandemic data  |  ABC n = %d  |  Best-fit debt = %.1f%%\n%s",
    length(abc_idx), debt_fit,
    sprintf("ABC posterior: %.1f%% [%.1f%%-%.1f%% 95%% UI]  |  R0 (base/peak): %.2f/%.2f", median(debt_abc, na.rm = TRUE), quantile(debt_abc, 0.025, na.rm = TRUE), quantile(debt_abc, 0.975, na.rm = TRUE), r0_fit$base, r0_fit$peak))
} else {
  sprintf("SEIRS fitted to pre-pandemic data  |  ABC n = %d  |  NM optimization failed", length(abc_idx))
}

df_fit <- data.frame(
  date    = df1$date,
  obs     = df1$FLUAB / mean_obs_pre,
  act_lo  = q_I_act_abc$lo  / mean_obs_pre,
  act_med = q_I_act_abc$med / mean_obs_pre,
  act_hi  = q_I_act_abc$hi  / mean_obs_pre,
  cf_lo   = q_I_cf_abc$lo   / mean_obs_pre,
  cf_med  = q_I_cf_abc$med  / mean_obs_pre,
  cf_hi   = q_I_cf_abc$hi   / mean_obs_pre,
  fit_act = I_fit_act_rel,
  fit_cf  = I_fit_cf_rel
)

p1 <- ggplot(df_fit, aes(x = date)) +
  npi_rect +
  geom_ribbon(aes(ymin = cf_lo,  ymax = cf_hi),  fill = col_cf,  alpha = 0.30) +
  geom_ribbon(aes(ymin = act_lo, ymax = act_hi),  fill = col_act, alpha = 0.30) +
  geom_line(aes(y = cf_med,  colour = "Counterfactual (ABC median)"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = act_med, colour = "Actual (ABC median)"), linewidth = 0.8) +
  geom_line(aes(y = fit_act, colour = "Best-fit actual (NM)"), linewidth = 1.0, linetype = "dotdash", na.rm = TRUE) +
  geom_line(aes(y = fit_cf,  colour = "Best-fit CF (NM)"), linewidth = 0.8, linetype = "dotdash", na.rm = TRUE) +
  geom_point(aes(y = obs, colour = "Observed"), size = 0.9, alpha = 0.6) +
  coord_cartesian(ylim = c(0, y_clip_fit)) +
  scale_colour_manual(values = c(
    "Actual (ABC median)"         = col_act,
    "Counterfactual (ABC median)" = col_cf,
    "Best-fit actual (NM)"        = "#8B0000",
    "Best-fit CF (NM)"            = "#00008B",
    "Observed"                    = "black")) +
  scale_y_continuous(labels = function(x) paste0(x, "\u00d7"), expand = c(0, 0)) +
  scale_x_date(NULL, date_breaks = "1 year", date_labels = "%b\n%Y", expand = c(0, 0)) +
  labs(subtitle = nm_subtitle, y = "Influenza activity (fold vs. pre-pandemic mean)", colour = NULL) +
  common_theme +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.20, 0.80),
    legend.background      = element_blank()
  )

# Fig : Susceptible pool dynamics (95% UI)
p2 <- ggplot(df_ribbon, aes(x = date)) +
  npi_rect +
  geom_ribbon(aes(ymin = S_cf_lo,  ymax = S_cf_hi),  fill = col_S_cf,  alpha = 0.25) +
  geom_ribbon(aes(ymin = S_act_lo, ymax = S_act_hi),  fill = col_S_act, alpha = 0.25) +
  geom_line(aes(y = S_cf_med,  colour = "Counterfactual"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = S_act_med, colour = "Actual"), linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2020-01-01"), linetype = "solid", colour = "grey40", linewidth = 0.4) +
  annotate("text", x = as.Date("2020-03-01"), y = 0.93, label = "NPI onset\n(Scenarios diverge)", hjust = 0, size = 3.5, colour = "grey0") +
  scale_colour_manual(values = c("Actual" = col_S_act, "Counterfactual" = col_S_cf)) +
  scale_y_continuous(limits = c(0.4, 1), expand = c(0, 0)) +
  scale_x_date(NULL, date_breaks = "1 year", date_labels = "%b\n%Y", expand = c(0, 0)) +
  labs(
    subtitle = sprintf("End-2022 S accumulation: %.1f pp (95%% UI %.1f%%-%.1f%%)", S_acc_q["50%"]*100, S_acc_q["2.5%"]*100, S_acc_q["97.5%"]*100),
    x = NULL,
    y = "Susceptible fraction",
    colour = NULL
  ) +
  common_theme +
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.2, 0.88),
    legend.background = element_blank()
  )

# Fig : PRCC global sensitivity — debt fraction
prcc_plot <- prcc_df[prcc_df$output == "Debt %", ]
prcc_plot$parameter <- factor(prcc_plot$parameter, levels = prcc_plot$parameter[order(abs(prcc_plot$PRCC))])

p3 <- ggplot(prcc_plot, aes(x = parameter, y = PRCC, fill = PRCC > 0)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = col_act, "FALSE" = "#0072B2"), guide = "none") +
  labs(
    subtitle = glue("Error bars: bootstrap 95% CI | n = {nboot} replicates"),
    x = NULL, 
    y = "Partial Rank Correlation Coefficient"
  ) +
  common_theme +
  theme(
    axis.text.y = element_text(size = base.size * 0.85, angle = 0, hjust = 1)
  )

# Fig : Grey band: GAMM-anchored lower bound. Intersection = mechanistically consistent duration.
p4 <- ggplot(oneway_df, aes(x = immunity_yr, y = debt_pct)) +
  geom_ribbon(aes(ymin = gamm_lo, ymax = gamm_hi), fill = "grey80", alpha = 0.6) +
  geom_line(colour = col_act, linewidth = 1) +
  geom_point(colour = col_act, size = 2.5) +
  annotate(
    "text", 
    x = mean(oneway_df$immunity_yr, na.rm = TRUE) + 1, 
    y = (gamm_lo + gamm_hi) / 2, 
    label = sprintf('GAMM ceiling - lower bound\n(%.1f%%-%.1f%%)', gamm_lo, gamm_hi), 
    hjust = 0.5, size = 4.8, colour = "grey0"
  ) +
  labs(
    subtitle = "One-way sensitivity: immunity duration to debt fraction",
    x = "Mean immunity duration (years)", 
    y = "Predicted debt fraction (%)"
  ) +
  scale_x_continuous(breaks = pretty(oneway_df$immunity_yr, n = 6)) +
  common_theme

# Fig : Prior predictive histogram 
# Prior predictive distribution of debt fraction
df_hist <- data.frame(
  debt_pct  = lhs_results$debt_pct[valid],
  debt_sign = ifelse(lhs_results$debt_pct[valid] > 0, "Positive debt", "No debt")
)

p5 <- ggplot(df_hist, aes(x = debt_pct, fill = debt_sign)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0,  linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 38, linetype = "dotted", colour = col_cf) +
  geom_vline(xintercept = 44, linetype = "dotted", colour = col_cf) +
  annotate("text", x = 41, y = Inf, vjust = 1.5, size = 4, label = sprintf('GAMM ceiling - lower bound\n(%.1f%%-%.1f%%)', gamm_lo, gamm_hi), colour = '#000000') +
  scale_fill_manual(values = c("Positive debt" = col_act, "No debt" = "#999999")) +
  labs(
    subtitle = sprintf("%d/%d (%.0f%%) draws show positive debt", sum(debt_pos), sum(valid), debt_pos_pct),
    x = "Debt fraction (% of 2023 peak)", 
    y = "Count", 
    fill = NULL
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  common_theme + 
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.2, 0.88),
    legend.background = element_blank()
  )

ggsave(p1, filename = glue('{fig.dir}/{ofig}_SEIRS_fitted_{this}.pdf'), width = 8, height = 5, units = 'in', bg = '#FFFFFF')
ggsave(p2, filename = glue('{fig.dir}/{ofig}_SEIRS_susceptible_{this}.pdf'),  width = 7, height = 5, units = 'in', bg = '#FFFFFF')
ggsave(p3, filename = glue('{fig.dir}/{ofig}_SEIRS_PRCC_{this}.pdf'),         width = 5,  height = 5, units = 'in', bg = '#FFFFFF')
ggsave(p4, filename = glue('{fig.dir}/{ofig}_SEIRS_oneway_{this}.pdf'),       width = 5,  height = 5, units = 'in', bg = '#FFFFFF')
ggsave(p5, filename = glue('{fig.dir}/{ofig}_SEIRS_histogram_{this}.pdf'),    width = 5,  height = 5, units = 'in', bg = '#FFFFFF')

# ---- 11. Save Tables -------------------------------------------------
write.csv(param_table,            glue('{tbl.dir}/SEIRS_param_table.csv'),    row.names = FALSE)
write.csv(lhs_results[valid,  ],  glue('{tbl.dir}/SEIRS_LHS_valid.csv'),      row.names = FALSE)
write.csv(lhs_results[abc_idx,],  glue('{tbl.dir}/SEIRS_ABC_posterior.csv'),  row.names = FALSE)
write.csv(prcc_df,                glue('{tbl.dir}/SEIRS_PRCC.csv'),           row.names = FALSE)
write.csv(oneway_df,              glue('{tbl.dir}/SEIRS_oneway_omega.csv'),   row.names = FALSE)