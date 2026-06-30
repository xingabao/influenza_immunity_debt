# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lubridate)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
this    <- 'Macau'

# Plotting constants
base.size   <- 18
base.family <- 'serif'
base.col    <- '#000000'
cols <- c(
  'Pre-COVID (Training)' = '#888888', 
  'COVID (Suppression)' = '#56B4E9',    
  'Post-COVID (Rebound)' = '#CC79A7'  
)

common_theme <- theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle = element_text(size = base.size * 0.9, color = base.col, margin = margin(r = 5)),
    axis.text.x   = element_text(size = base.size * 0.8, color = base.col),
    axis.text.y   = element_text(size = base.size * 0.8, color = base.col, angle = 90, hjust = 0.5),
    axis.title.y  = element_text(size = base.size * 0.9, face = 'bold', margin = margin(r = 5)),
    plot.margin   = margin(t = 5, b = 5, r = 10, l = 5)
  )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Data Loading and Feature Engineering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("[Step 1/5] Loading and preparing data ...")

dat.raw <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx'))

if (this == 'Macau') {
  dat.raw <- dat.raw %>%
    dplyr::select(date, FLUAB, nFLUAB, nSample) %>%
    filter(!is.na(FLUAB))
} else {
  dat.raw <- dat.raw %>%
    dplyr::select(date, FLUAB)
}

dat.model <- dat.raw %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(date)) %>%
  mutate(
    year      = year(date),
    week_num  = week(date),
    year_fact = factor(year)
  )

if (this == 'Macau') {
  dat.model <- dat.model %>%
    filter(!is.na(FLUAB)) %>%
    mutate(
      n_pos     = nFLUAB,
      n_tot     = nSample,
      obs_rate  = n_pos / n_tot,
      log_tests = log(n_tot + 1)
    )
} else {
  dat.model <- dat.model %>%
    mutate(
      obs_rate  = if (max(FLUAB, na.rm = TRUE) > 1) FLUAB / 100 else FLUAB,
      n_pos     = NA,
      n_tot     = NA,
      log_tests = 0
    )
}

dat.model <- dat.model %>%
  mutate(
    period = case_when(
      date < as.Date('2020-01-01')                                  ~ 'Pre-COVID (Training)',
      date >= as.Date('2020-01-01') & date <= as.Date('2023-01-01') ~ 'COVID (Suppression)',
      TRUE                                                          ~ 'Post-COVID (Rebound)'
    ),
    period = factor(period, levels = c('Pre-COVID (Training)', 'COVID (Suppression)', 'Post-COVID (Rebound)'))
  )

dat.train <- dat.model %>% filter(period == 'Pre-COVID (Training)')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Fitting GAMM model
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("[Step 2/5] Fitting GAMM model ...")

if (this == 'Macau') {
  f_gam <- as.formula(
    "cbind(n_pos, n_tot - n_pos) ~ s(week_num, bs = 'cc', k = 10) + s(log_tests, k = 5) + s(year_fact, bs = 're')"
  )
} else {
  f_gam <- as.formula(
    "obs_rate ~ s(week_num, bs = 'cc', k = 10) + s(year_fact, bs = 're')"
  )
}

gam_model_fixed <- gam(
  f_gam,
  family = quasibinomial(link = 'logit'),
  data   = dat.train,
  method = 'REML'
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Counterfactual simulation
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("[Step 3/5] Running counterfactual simulation (N = 1000 draws) ...")

set.seed(123)
N_sim <- 1000

dat.pred <- dat.model
if (this == 'Macau') {
  median_log_tests   <- median(dat.train$log_tests, na.rm = TRUE)
  dat.pred$log_tests <- median_log_tests
}

Xp      <- predict(gam_model_fixed, newdata = dat.pred, type = 'lpmatrix')
re_cols <- grep("s\\(year_fact\\)", colnames(Xp))
if (length(re_cols) > 0) Xp[, re_cols] <- 0

beta          <- coef(gam_model_fixed)
Vb            <- vcov(gam_model_fixed)
mrand         <- mvrnorm(N_sim, beta, Vb)
pred_link_sim <- Xp %*% t(mrand)
pred_resp_sim <- plogis(pred_link_sim)

dat.result <- dat.model %>%
  mutate(
    pred_rate = rowMeans(pred_resp_sim),
    ci_lower  = apply(pred_resp_sim, 1, quantile, probs = 0.025),
    ci_upper  = apply(pred_resp_sim, 1, quantile, probs = 0.975),
    diff_rate = obs_rate - pred_rate
  )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4. Magnitude Ceiling Test + Drift vs. Immunity Debt Decomposition
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("[Step 4/5] Running Magnitude Ceiling Test + Drift / Immunity Debt Decomposition ...")

# Per-year posterior summaries
# --------------------------------
all_years <- sort(unique(year(dat.result$date)))

annual_list <- lapply(all_years, function(yr) {
  
  idx <- which(year(dat.result$date) == yr)
  if (length(idx) < 20) return(NULL)
  
  obs_vec    <- dat.result$obs_rate[idx]
  sim_mat    <- pred_resp_sim[idx, , drop = FALSE]
  excess_mat <- matrix(obs_vec, nrow = length(idx), ncol = N_sim) - sim_mat
  
  mean_sims <- colMeans(excess_mat, na.rm = TRUE)
  peak_sims <- apply(excess_mat, 2, max, na.rm = TRUE)
  
  breadth_sims <- colMeans(excess_mat > 0, na.rm = TRUE)
  
  period_yr <- dat.result %>%
    filter(year(date) == yr) %>%
    count(period, sort = TRUE) %>%
    slice(1) %>% pull(period) %>% as.character()
  
  list(
    year         = yr,
    period       = period_yr,
    n_weeks      = length(idx),
    mean_exc_est = mean(mean_sims),
    mean_exc_lwr = quantile(mean_sims, 0.025),
    mean_exc_upr = quantile(mean_sims, 0.975),
    peak_exc_est = mean(peak_sims),
    peak_exc_lwr = quantile(peak_sims, 0.025),
    peak_exc_upr = quantile(peak_sims, 0.975),
    breadth_est  = mean(breadth_sims),
    breadth_lwr  = quantile(breadth_sims, 0.025),
    breadth_upr  = quantile(breadth_sims, 0.975),
    mean_sims    = list(mean_sims),
    peak_sims    = list(peak_sims),
    breadth_sims = list(breadth_sims)
  )
})

df_annual <- bind_rows(Filter(Negate(is.null), annual_list)) %>%
  mutate(
    period   = factor(period, levels = names(cols)),
    is_pre   = period == 'Pre-COVID (Training)',
    is_post  = period == 'Post-COVID (Rebound)',
    is_covid = year %in% c(2020, 2021, 2022)
  )

# Historical ceiling: pre-pandemic pooled reference distributions
# --------------------------------
df_pre  <- filter(df_annual, is_pre)
df_post <- filter(df_annual, is_post)

pool_pre_peak    <- unlist(df_pre$peak_sims)
pool_pre_mean    <- unlist(df_pre$mean_sims)
pool_pre_breadth <- unlist(df_pre$breadth_sims)

cln_pk_50   <- quantile(pool_pre_peak, 0.500)
cln_pk_90   <- quantile(pool_pre_peak, 0.900)
cln_pk_95   <- quantile(pool_pre_peak, 0.950)
cln_pk_975  <- quantile(pool_pre_peak, 0.975)
cln_mn_975  <- quantile(pool_pre_mean, 0.975)

# Season breadth ceiling: above 97.5th pct → unusually widespread excess
cln_breadth_975 <- quantile(pool_pre_breadth, 0.975)

mu_pre_pk      <- mean(pool_pre_peak);    sd_pre_pk      <- sd(pool_pre_peak)
mu_pre_mn      <- mean(pool_pre_mean);    sd_pre_mn      <- sd(pool_pre_mean)
mu_pre_breadth <- mean(pool_pre_breadth); sd_pre_breadth <- sd(pool_pre_breadth)

yr2019 <- filter(df_annual, year == 2019)

# Post-pandemic statistical inference per year
# --------------------------------
inf_list <- lapply(seq_len(nrow(df_post)), function(i) {
  
  pk_sims <- df_post$peak_sims[[i]]
  mn_sims <- df_post$mean_sims[[i]]
  br_sims <- df_post$breadth_sims[[i]]
  
  fold_2019 <- if (isTRUE(yr2019$peak_exc_est > 0))
    df_post$peak_exc_est[i] / yr2019$peak_exc_est else NA_real_
  
  data.frame(
    year                    = df_post$year[i],
    peak_exc_est            = df_post$peak_exc_est[i],
    peak_exc_lwr            = df_post$peak_exc_lwr[i],
    peak_exc_upr            = df_post$peak_exc_upr[i],
    mean_exc_est            = df_post$mean_exc_est[i],
    mean_exc_lwr            = df_post$mean_exc_lwr[i],
    mean_exc_upr            = df_post$mean_exc_upr[i],
    breadth_est             = df_post$breadth_est[i],
    breadth_lwr             = df_post$breadth_lwr[i],
    breadth_upr             = df_post$breadth_upr[i],
    p_peak_above_ceiling    = mean(pk_sims > cln_pk_975),
    p_mean_above_ceiling    = mean(mn_sims > cln_mn_975),
    p_breadth_above_ceiling = mean(br_sims > cln_breadth_975),
    z_peak                  = (df_post$peak_exc_est[i] - mu_pre_pk)      / sd_pre_pk,
    z_mean                  = (df_post$mean_exc_est[i] - mu_pre_mn)      / sd_pre_mn,
    z_breadth               = (df_post$breadth_est[i]  - mu_pre_breadth) / sd_pre_breadth,
    pctrank_peak            = mean(pool_pre_peak    < df_post$peak_exc_est[i]),
    pctrank_breadth         = mean(pool_pre_breadth < df_post$breadth_est[i]),
    fold_over_2019          = fold_2019
  )
})
df_inf <- bind_rows(inf_list)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5. Building plots
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("[Step 5/5] Building plots ...")

col_pre    <- 'gray62'
col_post   <- '#CC79A7'
col_covid  <- '#56B4E9'
col_ceil   <- '#1B3A8C'
col_2019   <- '#E07B00'
col_median <- 'gray42'

df_pA <- df_annual %>%
  mutate(
    bar_fill  = case_when(
      year == 2019 ~ col_2019,
      is_covid     ~ 'gray82',
      is_post      ~ col_post,
      TRUE         ~ col_pre
    ),
    bar_alpha = ifelse(is_covid, 0.40, 0.88),
    label_key = year == 2019 | is_post
  )

ymax_A <- max(df_pA$peak_exc_upr, na.rm = TRUE)
ymin_A <- min(df_pA$peak_exc_lwr, na.rm = TRUE)
yr_rng <- ymax_A - ymin_A
yr_min <- min(all_years)

# Panel A: Magnitude Ceiling Test
gg <- ggplot(df_pA, aes(x = year, y = peak_exc_est)) +
  
  annotate("rect", xmin = 2019.55, xmax = 2022.45, ymin = -Inf, ymax = Inf, fill = col_covid, alpha = 0.07) +
  annotate("text", x = 2021, y = ymin_A + yr_rng * 0.02, label = "NPI suppression\n(excluded from ceiling)", size = base.size / 5.0, color = 'gray28', family = base.family, hjust = 0.5, lineheight = 0.80) +
  
  geom_hline(yintercept = cln_pk_975, linetype = 'dashed', color = col_ceil,   linewidth = 0.90) +
  geom_hline(yintercept = cln_pk_50,  linetype = 'dotted', color = col_median, linewidth = 0.65) +
  
  geom_col(aes(fill = I(bar_fill), alpha = I(bar_alpha)), width = 0.70) +
  geom_errorbar(aes(ymin = peak_exc_lwr, ymax = peak_exc_upr), width = 0.18, color = 'gray12', linewidth = 0.40) +
  
  geom_text(
    data = filter(df_pA, label_key),
    aes(y = peak_exc_upr + yr_rng * c(0.020, 0.020, 0.040, 0.020), label = sprintf("%+.1f%%", peak_exc_est * 100)),
    size = base.size / 4.5, family = base.family, color = '#000000', vjust = 0
  ) +
  
  annotate(
    "text", x = yr_min + 0.2, y = cln_pk_975 + yr_rng * 0.050,
    label = "Historical ceiling (97.5th pct of pre-pandemic distribution)",
    hjust = 0, size = base.size / 4.0, color = col_ceil, family = base.family
  ) +
  annotate(
    "text", x = yr_min + 0.2, y = cln_pk_50 - yr_rng * 0.32,
    label = "Historical median",
    hjust = 0, size = base.size / 4.0, color = 'gray28', family = base.family
  ) +
  annotate(
    "text", x = 2019, y = filter(df_pA, year == 2019)$peak_exc_upr + yr_rng * 0.13,
    label = "Strong\ndrift year",
    hjust = 0.5, size = base.size / 4.5, lineheight = 0.82,
    color = col_2019, family = base.family, fontface = 'italic'
  ) +
  
  scale_x_continuous(breaks = seq(yr_min, max(all_years), by = 1), expand = c(0.015, 0.015)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.08, 0.22))) +
  labs(x = NULL, y = "Peak Annual Excess Rate") +
  common_theme +
  theme(panel.grid.minor = element_blank())

ggsave(gg, filename = glue('{fig.dir}/Magnitude_ceiling_test_{this}.pdf'), width = 11, height = 5, units = 'in', bg = 'white')
message(glue("[Saved] MagnitudeCeilingTest_{this}.pdf"))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6. Tables
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Table 1: Annual Summary
# --------------------------------
df_pre_z <- df_pre %>%
  mutate(
    z_peak_score            = (peak_exc_est - mu_pre_pk)      / sd_pre_pk,
    z_breadth_score         = (breadth_est  - mu_pre_breadth) / sd_pre_breadth,
    p_peak_above_ceiling    = NA_real_,
    p_breadth_above_ceiling = NA_real_,
    pctrank_peak            = NA_real_,
    pctrank_breadth         = NA_real_,
    fold_over_2019          = NA_real_
  )

df_post_z <- df_inf %>%
  rename(z_peak_score = z_peak, z_breadth_score = z_breadth)

df_z_all <- bind_rows(
  df_pre_z  %>% dplyr::select(year, z_peak_score, z_breadth_score, p_peak_above_ceiling, p_breadth_above_ceiling, pctrank_peak, pctrank_breadth, fold_over_2019),
  df_post_z %>% dplyr::select(year, z_peak_score, z_breadth_score, p_peak_above_ceiling, p_breadth_above_ceiling, pctrank_peak, pctrank_breadth, fold_over_2019)
)

tbl_annual <- df_annual %>%
  filter(!is_covid) %>%
  left_join(df_z_all, by = "year") %>%
  dplyr::select(
    year, period, n_weeks,
    peak_exc_est, peak_exc_lwr, peak_exc_upr,
    mean_exc_est, mean_exc_lwr, mean_exc_upr,
    breadth_est,  breadth_lwr,  breadth_upr,
    z_peak_score, z_breadth_score,
    p_peak_above_ceiling, p_breadth_above_ceiling,
    pctrank_peak, pctrank_breadth, fold_over_2019
  ) %>%
  mutate(
    across(c(peak_exc_est, peak_exc_lwr, peak_exc_upr,
             mean_exc_est, mean_exc_lwr, mean_exc_upr), ~ round(.x * 100, 3)),
    across(c(breadth_est, breadth_lwr, breadth_upr),    ~ round(.x * 100, 3)),
    across(c(z_peak_score, z_breadth_score, fold_over_2019), ~ round(.x, 3)),
    across(c(p_peak_above_ceiling, p_breadth_above_ceiling,
             pctrank_peak, pctrank_breadth), ~ round(.x, 4))
  ) %>%
  rename(
    Year                    = year,
    Period                  = period,
    N_weeks                 = n_weeks,
    Peak_Excess_pct         = peak_exc_est,
    Peak_Excess_LCL_pct     = peak_exc_lwr,
    Peak_Excess_UCL_pct     = peak_exc_upr,
    Mean_Excess_pct         = mean_exc_est,
    Mean_Excess_LCL_pct     = mean_exc_lwr,
    Mean_Excess_UCL_pct     = mean_exc_upr,
    Season_Breadth_pct      = breadth_est,
    Season_Breadth_LCL_pct  = breadth_lwr,
    Season_Breadth_UCL_pct  = breadth_upr,
    Z_score_peak            = z_peak_score,
    Z_score_breadth         = z_breadth_score,
    P_peak_above_ceiling    = p_peak_above_ceiling,
    P_breadth_above_ceiling = p_breadth_above_ceiling,
    Pctrank_peak            = pctrank_peak,
    Pctrank_breadth         = pctrank_breadth,
    Fold_vs_2019            = fold_over_2019
  ) %>%
  arrange(Year)

# Table 2: Ceiling Reference Values
# --------------------------------
tbl_ceiling_ref <- data.frame(
  Metric = c(
    "Location",
    "N pre-pandemic years (ceiling basis)",
    "N pooled posterior draws",
    "Historical peak excess — Median [%]",
    "Historical peak excess — 90th pct [%]",
    "Historical peak excess — 95th pct [%]",
    "Historical peak excess — 97.5th pct / Ceiling [%]",
    "Historical mean excess — 97.5th pct [%]",
    "Historical peak excess — Mean [%]",
    "Historical peak excess — SD [%]",
    "Historical season breadth — 97.5th pct / Ceiling [%]",
    "Historical season breadth — Mean [%]",
    "Historical season breadth — SD [%]",
    "2019 peak excess — Estimate [%]",
    "2019 peak excess — 95% CI lower [%]",
    "2019 peak excess — 95% CI upper [%]"
  ),
  Value = c(
    this,
    as.character(nrow(df_pre)),
    as.character(length(pool_pre_peak)),
    as.character(round(cln_pk_50       * 100, 3)),
    as.character(round(cln_pk_90       * 100, 3)),
    as.character(round(cln_pk_95       * 100, 3)),
    as.character(round(cln_pk_975      * 100, 3)),
    as.character(round(cln_mn_975      * 100, 3)),
    as.character(round(mu_pre_pk       * 100, 3)),
    as.character(round(sd_pre_pk       * 100, 3)),
    as.character(round(cln_breadth_975 * 100, 3)),
    as.character(round(mu_pre_breadth  * 100, 3)),
    as.character(round(sd_pre_breadth  * 100, 3)),
    as.character(round(yr2019$peak_exc_est * 100, 3)),
    as.character(round(yr2019$peak_exc_lwr * 100, 3)),
    as.character(round(yr2019$peak_exc_upr * 100, 3))
  )
)

writexl::write_xlsx(
  list(
    "Annual_Summary"    = tbl_annual,
    "Ceiling_Reference" = tbl_ceiling_ref
  ),
  path = glue('{tbl.dir}/MagnitudeCeilingTest_{this}.xlsx')
)

message("\n[Magnitude Ceiling Test] All done.\n")
