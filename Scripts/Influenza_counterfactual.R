# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(scales))) 
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(patchwork))) 
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'HK'  # Options: 'HK' or 'Macau'

# Plotting Constants
base.size <- 7
base.family <- 'sans'
base.col <- '#000000'
cols <- c(
  'Pre-COVID (Training)' = '#888888', 
  'COVID (Suppression)' = '#56B4E9',    
  'Post-COVID (Rebound)' = '#CC79A7'    
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Data Loading and Feature Engineering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dat.raw <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx'))
if (this == 'Macau') {
  dat.raw <- dat.raw %>% dplyr::select(date, FLUAB, nFLUAB, nSample) %>% filter(!is.na(dat.raw$FLUAB))
} else {
  dat.raw <- dat.raw %>% dplyr::select(date, FLUAB)
}

# Normalize and Create Features
dat.model <- dat.raw %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(date)) %>%
  mutate(
    year = year(date),
    week_num = week(date),
    # Year as factor for Random Effects
    year_fact = factor(year)
  )

if (this == 'Macau') {
  dat.model <- dat.model %>%
    filter(!is.na(FLUAB)) %>%
    mutate(
      n_pos = nFLUAB,
      n_tot = nSample,
      obs_rate = n_pos / n_tot,
      # Log testing volume for controlling testing behavior
      log_tests = log(n_tot + 1) 
    )
} else {
  # Handling HK data
  dat.model <- dat.model %>%
    mutate(
      obs_rate = if (max(FLUAB, na.rm = TRUE) > 1) FLUAB / 100 else FLUAB,
      n_pos = NA, 
      n_tot = NA,
      log_tests = 0 
    )
}

# Define Periods
dat.model <- dat.model %>%
  mutate(
    period = case_when(
      date < as.Date('2020-01-01') ~ 'Pre-COVID (Training)',
      date >= as.Date('2020-01-01') & date <= as.Date('2023-01-01') ~ 'COVID (Suppression)',
      TRUE ~ 'Post-COVID (Rebound)'
    ),
    period = factor(period, levels = c('Pre-COVID (Training)', 'COVID (Suppression)', 'Post-COVID (Rebound)'))
  )

dat.train <- dat.model %>% filter(period == 'Pre-COVID (Training)')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Main GAMM Modeling
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (this == 'Macau') {
  # Binomial count data with volume control + Random Effects
  f_gam <- as.formula("cbind(n_pos, n_tot - n_pos) ~ s(week_num, bs = 'cc', k = 10) + s(log_tests, k = 5) + s(year_fact, bs = 're')")
} else {
  # Rate data (HK) + Random Effects
  f_gam <- as.formula("obs_rate ~ s(week_num, bs = 'cc', k = 10) + s(year_fact, bs = 're')")
}

message("Fitting Main GAMM Model...")

gam_model_fixed <- gam(
  f_gam,
  family = quasibinomial(link = 'logit'), 
  data = dat.train,
  method = 'REML'
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Counterfactual Simulation (The Baseline)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set.seed(123)
N_sim <- 1000

# Prepare Prediction Data
dat.pred <- dat.model

# 1. Fix testing volume to Pre-COVID median (remove testing behavior bias)
if (this == 'Macau') {
  median_log_tests <- median(dat.train$log_tests, na.rm = TRUE)
  dat.pred$log_tests <- median_log_tests
}

# 2. Generate Prediction Matrix
Xp <- predict(gam_model_fixed, newdata = dat.pred, type = 'lpmatrix')

# 3. Zero out Random Effects (Simulate an "Average" Year)
re_cols <- grep("s\\(year_fact\\)", colnames(Xp))
if (length(re_cols) > 0) {
  Xp[, re_cols] <- 0
}

# Simulation
beta <- coef(gam_model_fixed)
Vb <- vcov(gam_model_fixed)
mrand <- mvrnorm(N_sim, beta, Vb)

pred_link_sim <- Xp %*% t(mrand)
pred_resp_sim <- plogis(pred_link_sim)

# Summarize results
dat.result <- dat.model %>%
  mutate(
    pred_rate = rowMeans(pred_resp_sim),
    ci_lower  = apply(pred_resp_sim, 1, quantile, probs = 0.025),
    ci_upper  = apply(pred_resp_sim, 1, quantile, probs = 0.975),
    diff_rate = obs_rate - pred_rate
  )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4. Sensitivity Analysis / Negative Control
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
message("Running Sensitivity Analysis (Negative Control)...")

valid_start <- as.Date('2018-01-01')
valid_end   <- as.Date('2020-01-01')

dat.sens.train <- dat.train %>% filter(date < valid_start)
dat.sens.valid <- dat.train %>% filter(date >= valid_start & date < valid_end)

has_sensitivity <- nrow(dat.sens.valid) > 10

if (has_sensitivity) {
  gam_sens <- gam(f_gam, family = quasibinomial(link = 'logit'), data = dat.sens.train, method = 'REML')
  
  dat.sens.pred <- dat.sens.valid
  if (this == 'Macau') dat.sens.pred$log_tests <- median(dat.sens.train$log_tests, na.rm = TRUE)
  
  Xp_sens <- predict(gam_sens, newdata = dat.sens.pred, type = 'lpmatrix')
  
  re_cols_sens <- grep("s\\(year_fact\\)", colnames(Xp_sens))
  if (length(re_cols_sens) > 0) Xp_sens[, re_cols_sens] <- 0
  
  beta_sens <- coef(gam_sens)
  mrand_sens <- mvrnorm(1000, beta_sens, vcov(gam_sens))
  pred_sens_mat <- plogis(Xp_sens %*% t(mrand_sens))
  
  dat.sens.valid$pred_mean <- rowMeans(pred_sens_mat)
  dat.sens.valid$ci_lwr <- apply(pred_sens_mat, 1, quantile, probs = 0.025)
  dat.sens.valid$ci_upr <- apply(pred_sens_mat, 1, quantile, probs = 0.975)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5. Statistical Inference (Quarterly Trend)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
target_period <- c('COVID (Suppression)', 'Post-COVID (Rebound)')
dat.post <- dat.result %>% filter(period %in% target_period) %>% arrange(date)
df_trend <- data.frame()

# Exact P value formatting: plain text for source data, plotmath for figure labels
format_p_exact <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < .Machine$double.eps) {
      return('2.220 × 10^-16')
    }
    if (x >= 0.001) {
      return(formatC(x, format = 'f', digits = 3))
    }
    exponent <- floor(log10(x))
    sprintf('%.3f × 10^%d', x / 10^exponent, exponent)
  }, character(1))
}

format_p_plotmath <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < .Machine$double.eps) {
      return('italic(P) == 2.220 %*% 10^-16')
    }
    if (x >= 0.001) {
      return(sprintf('italic(P) == %s', formatC(x, format = 'f', digits = 3)))
    }
    exponent <- floor(log10(x))
    sprintf('italic(P) == %.3f %*% 10^%d', x / 10^exponent, exponent)
  }, character(1))
}

if (nrow(dat.post) > 0) {
  dat.post <- dat.post %>%
    mutate(
      time_block = cut(date, breaks = '3 months', start.on.monday = FALSE),
      block_start = as.Date(as.character(time_block))
    )
  
  unique_blocks <- unique(dat.post$block_start)
  results_list <- list()
  
  for (i in seq_along(unique_blocks)) {
    b_start <- unique_blocks[i]
    sub_dat <- dat.post %>% filter(block_start == b_start)
    original_indices <- which(dat.model$date %in% sub_dat$date)
    
    obs_vec <- sub_dat$obs_rate
    sim_mat <- pred_resp_sim[original_indices, , drop = FALSE]
    diff_mat <- sweep(sim_mat * -1, 1, obs_vec, '+')
    block_mean_diffs <- colMeans(diff_mat, na.rm = TRUE)
    
    est_diff <- mean(block_mean_diffs)
    ci_diff  <- quantile(block_mean_diffs, probs = c(0.025, 0.975))
    p_val    <- sum(block_mean_diffs <= 0) / N_sim 
    
    results_list[[i]] <- data.frame(
      block_start = b_start,
      excess_est  = est_diff,
      excess_lwr  = ci_diff[1],
      excess_upr  = ci_diff[2],
      p_value     = p_val,
      p_value_exact = format_p_exact(p_val),
      is_signif   = p_val < 0.05,
      p_label     = format_p_plotmath(p_val)
    )
  }
  df_trend <- do.call(rbind, results_list)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6. Visualization
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rects <- dat.result %>%
  group_by(period) %>%
  summarise(xmin = min(date), xmax = max(date), .groups = 'drop')

common_theme <- theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle = element_text(size = base.size, color = base.col, margin = margin(r = 5)),
    axis.text.x = element_text(size = base.size, color = base.col),
    axis.text.y = element_text(size = base.size, color = base.col, angle = 90, hjust = 0.5),
    axis.title.y = element_text(size = base.size, face = 'bold', margin = margin(r = 5)),
    plot.margin = margin(t = 5, b = 5, r = 10, l = 5)
  )

# Plot A: Time Series
p_ts <- ggplot() +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period), alpha = 0.2, inherit.aes = FALSE) +
  geom_ribbon(data = dat.result, aes(x = date, ymin = ci_lower, ymax = ci_upper), fill = '#00008B', alpha = 0.2) +
  geom_line(data = dat.result, aes(x = date, y = pred_rate, linetype = 'Counterfactual Baseline'), color = '#00008B', linewidth = 0.4) +
  geom_line(data = dat.result, aes(x = date, y = obs_rate, linetype = 'Observed Data'), color = 'black', linewidth = 0.5, alpha = 0.8) +
  scale_fill_manual(values = cols, guide = 'none') +
  scale_linetype_manual(name = NULL, values = c('Observed Data' = 'solid', 'Counterfactual Baseline' = 'dashed')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(NULL, date_breaks = "1 year", date_labels = "%b\n%Y", expand = c(0, 0)) +
  labs(
    subtitle = 'Observed vs. Expected',
    y = 'Positivity Rate',
    x = NULL
  ) +
  common_theme +
  theme(
    legend.position = c(ifelse(test = this == 'Macau', yes = 0.20, no = 0.98), 0.98), 
    legend.justification = c(1, 1),
    legend.key.spacing.y = unit(-0.20, units = 'cm'),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.text.x = element_blank()
  )

# Plot B: Weekly Difference
p_diff <- ggplot(dat.result, aes(x = date, y = diff_rate)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period), alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = 'gray30', linewidth = 0.5) +
  geom_col(aes(fill = period, color = period), width = 7, alpha = 0.6) + 
  geom_smooth(method = 'loess', span = 0.1, color = 'black', se = FALSE, linewidth = 0.5) +
  scale_fill_manual(values = cols, guide = 'none') +
  scale_color_manual(values = cols, guide = 'none') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(NULL, date_breaks = "1 year", date_labels = "%b\n%Y", expand = c(0, 0)) +
  labs(
    subtitle = 'Weekly Excess Positivity',
    y = 'Difference', 
    x = NULL
  ) +
  common_theme

# Plot C: Quarterly Trend
line.colors <- df_trend$block_start
date_range <- range(line.colors)
calc_breaks <- breaks_width("3 months")
final_breaks <- calc_breaks(date_range)
calc_labels <- label_date(format = "%b\n%Y")
final_labels <- calc_labels(final_breaks)
final_labels[length(final_breaks)] <- ""
coll <- rep(base.col, length(final_labels))
coll[which(final_breaks %in% c('2023-01-01', '2023-04-01', '2023-07-01', '2023-10-01', '2024-01-01'))] <- '#E50914'

if (nrow(df_trend) > 0) {
  trend_breaks <- breaks_width("3 months")(range(df_trend$block_start))
  trend_labels <- label_date(format = "%b\n%Y")(trend_breaks)
  trend_labels[length(trend_labels)] <- "" 
  trend_labels[seq(2, length(trend_labels), by = 2)] <- ""
  
  p_trend <- ggplot(df_trend, aes(x = block_start, y = excess_est)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
    geom_errorbar(aes(ymin = excess_lwr, ymax = excess_upr), width = 20, color = '#E50914') +
    geom_line(color = '#E50914', group = 1, alpha = 0.5) +
    geom_point(aes(size = is_signif, fill = is_signif), shape = 21, color = '#E50914') +
    geom_text(aes(label = scales::percent(excess_est, 0.1), vjust = ifelse(this == 'HK', -2.5, -3.5)), size = base.size / ggplot2::.pt) +
    geom_text(data = df_trend %>% filter(p_value < 0.05), aes(label = p_label, y = excess_lwr, vjust = 1.5), size = base.size / ggplot2::.pt, color = 'gray30', family = base.family, parse = TRUE) +
    geom_rect(data = data.frame(xmin = -Inf, xmax = as.Date('2023-01-01')), aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = '#56B4E9', alpha = 0.2, inherit.aes = FALSE) +
    geom_rect(data = data.frame(xmin = as.Date('2023-01-01'), xmax = as.Date('2024-01-01')), aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = 0), fill = c('#141473'), alpha = 0.2, inherit.aes = FALSE) +
    geom_rect(data = data.frame(xmin = as.Date('2023-01-01'), xmax = Inf), aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = c('#CC79A7'), alpha = 0.2, inherit.aes = FALSE) +
    scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 1.5), guide = 'none') +
    scale_fill_manual(values = c(`FALSE` = 'white', `TRUE` = '#E50914'), guide = 'none') +
    scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = 0.2)) +
    scale_x_date(limits = as.Date(c('2019-10-15', '2025-08-15')), breaks = trend_breaks, labels = trend_labels, expand = c(0, 0)) +
    labs(
      subtitle = 'Trend of immunological debt (3-month intervals)',
      y = 'Mean Excess Rate',
      x = NULL
    ) +
    common_theme + theme(
      axis.text.x = element_text(size = base.size, family = base.family, color = coll)
    )
} else {
  p_trend <- ggplot() + theme_void()
}

if (this == 'HK') {
  p_trend <-  p_trend + 
    annotate("text", x = as.Date("2023-07-01"), y = -0.15, size = base.size / ggplot2::.pt, label = "Immune Debt", family = base.family, fontface  = 'bold') +
    annotate("text", x = as.Date("2021-06-01"), y = 0.04, size = base.size / ggplot2::.pt, label = "Pandemic Period", family = base.family) +
    annotate("text", x = as.Date("2024-12-01"), y = -0.20, size = base.size / ggplot2::.pt, label = "Post-pandemic Period", family = base.family)
} else {
  p_trend <-  p_trend + 
    annotate("text", x = as.Date("2023-07-01"), y = -0.10, size = base.size / ggplot2::.pt, label = "Immune Debt", family = base.family, fontface  = 'bold') +
    annotate("text", x = as.Date("2021-06-01"), y = 0.15, size = base.size / ggplot2::.pt, label = "Pandemic Period", family = base.family) +
    annotate("text", x = as.Date("2024-12-01"), y = 0.15, size = base.size / ggplot2::.pt, label = "Post-pandemic Period", family = base.family)
}

# Plot D: Sensitivity Analysis (Negative Control)
if (has_sensitivity) {
  
  trend_breaks <- as.Date(c('2018-01-01', '2018-07-01', '2019-01-01', '2019-07-01', '2020-01-01'))
  trend_labels <- label_date(format = "%b\n%Y")(trend_breaks)
  trend_labels[length(trend_labels)] <- ""
  
  p_sens <- ggplot() +
    geom_ribbon(data = dat.sens.valid, aes(x = date, ymin = ci_lwr, ymax = ci_upr), fill = 'gray50', alpha = 0.3) +
    geom_line(data = dat.sens.valid, aes(x = date, y = pred_mean, color = 'Predicted (Hold-out)'), linetype = 'dashed') +
    geom_line(data = dat.sens.valid, aes(x = date, y = obs_rate, color = 'Observed'), linewidth = 0.5) +
    scale_color_manual(name = NULL, values = c('Observed' = 'black', 'Predicted (Hold-out)' = 'blue')) +
    scale_x_date(limits = as.Date(c('2018-01-01', '2020-01-01')), breaks = trend_breaks, labels = trend_labels, expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      subtitle = 'Sensitivity Analysis: Negative Control (Predicting 2018-2019)', 
      y = 'Positivity Rate', 
      x = 'Date'
    ) +
    common_theme +
    theme(
      legend.position = c(0.98, 0.98), 
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha('white', 0.5)),
      legend.key.spacing.y = unit(-0.20, units = 'cm'),
      axis.title.x.bottom = element_text(size = base.size, family = base.family, color = base.col, face = 'bold', margin = margin(t = 3))
    )
} else {
  p_sens <- ggplot() + theme_void() + labs(subtitle = "Insufficient data for Sensitivity Analysis")
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 7. Final Output
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Assemble
layout <- "
AAAA
BBBB
CCCC
DDDD
"
gg <- p_ts / p_diff / p_trend / p_sens + 
  plot_layout(heights = c(0.8, 0.8, 1.2, 0.8)) + theme(plot.margin = margin(b = -10))

print(summary(gam_model_fixed))
print(gg)

# Save to file
width = 180; height = 165
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}.pdf'), width = width, height = height, units = 'mm', device = cairo_pdf, bg = '#FFFFFF')

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
if (this == 'HK') {
  fig.label <- 'Fig.S5'
} else {
  fig.label <- 'Fig.2'
}

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  writeLines(c(
    glue('Source data for {fig.label}{panel}: {title}')
  ), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
  close(con)
}

write_source_txt(
  dat.result, panel = 'A',
  title = 'Observed weekly influenza positivity rates versus GAMM counterfactual baseline with 95% simulation-based credible interval'
)

write_source_txt(
  dat.result, panel = 'B',
  title = 'Weekly excess positivity (observed minus counterfactual baseline) with Loess smoothed trend'
)

write_source_txt(
  df_trend, panel = 'C',
  title = 'Quarterly quantification of immunity debt: mean excess rate with 95% simulation-based confidence interval and exact one-sided P values'
)

if (has_sensitivity) {
  write_source_txt(
    dat.sens.valid, panel = 'D',
    title = 'Negative control out-of-sample validation: model trained on 2010-2017 projected onto held-out 2018-2019 seasons'
  )
}