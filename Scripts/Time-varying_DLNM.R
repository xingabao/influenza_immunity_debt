# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(lmtest)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(splines)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'Macau'   # Options: 'HK' or 'Macau'
base.size <- 16
base.family <- 'serif'
base.col <- '#000000'

# Load data
file.path <- glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')

# Temporal definitions for analysis
sdate <- '2020-01-01'
edate <- '2023-01-01'
cut.date <- as.Date('2023-01-01')

# Model hyperparameters
meteo.main <- 'ave.temp' 
max.lag <- 7
df.var  <- 3 
df.lag  <- 3
df.time <- 7 * 11

# Target outcomes
target.list <- c('FLUA', 'FLUB')

# ==============================================================================
# Section 1: Data Pre-processing and Feature Engineering
# ==============================================================================
prep_data_from_excel <- function(
    file.path, 
    cut.date,
    this = 'Macau', 
    sdate = as.Date('2020-01-01'), 
    edate = as.Date('2023-01-01')
) {
  
  if (!file.exists(file.path)) stop(glue("Error: Input file not found at {file.path}"))
  
  dat.raw <- readxl::read_excel(file.path)
  if (!inherits(dat.raw$date, 'Date')) dat.raw$date <- as.Date(dat.raw$date)
  
  # Calculate diurnal temperature range if not present
  if (!"dif.temp" %in% names(dat.raw) && "max.temp" %in% names(dat.raw)) {
    dat.raw$dif.temp <- dat.raw$max.temp - dat.raw$min.temp
  } else if (!"dif.temp" %in% names(dat.raw)) {
    dat.raw$dif.temp <- 0
  }
  
  # Type coercion for consistency
  dat.raw$FLUA <- as.numeric(dat.raw$FLUA)
  dat.raw$FLUB <- as.numeric(dat.raw$FLUB)
  
  # Define epidemiological periods and time indices
  dat.pl <- dat.raw %>%
    arrange(date) %>%
    mutate(
      Period = case_when(
        date < sdate ~ 'Pre-COVID',
        date > edate ~ 'Post-COVID',
        TRUE ~ 'COVID'
      ),
      post2023 = ifelse(date >= cut.date, 1, 0),
      t = as.numeric(date - min(date, na.rm = TRUE)),
      yday = as.numeric(format(date, "%j"))
    )
  
  # Filter for study specific timeframe
  if (this == 'Macau') dat.pl <- dat.pl %>% filter(date > as.Date('2014-12-31'))
  
  # Select relevant covariates
  vars_needed <- c(
    'date', 't', 'yday', 'post2023', 'Period',
    'FLUA', 'FLUB',
    'PM10', 'PM2.5', 'SO2', 'NO2', 'O3', 'CO',
    'dif.temp', 'ave.temp', 'dew.temp', 'pressure', 'humidity', 'sunshine', 'wind.speed', 'precipitation'
  )
  
  vars_exist <- intersect(vars_needed, names(dat.pl))
  dat.m <- dat.pl %>% dplyr::select(all_of(vars_exist)) %>% na.omit() %>% arrange(date)
  
  list(dat.m = dat.m)
}

add_time_factors <- function(df) {
  df %>% mutate(
    dow = factor(
      weekdays(date),
      levels = c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday')
    )
  )
}

cat(glue('\n[INFO] Loading and processing dataset: {file.path} ...\n'))
dat_list <- prep_data_from_excel(file.path, cut.date, this)
dat.m <- dat_list$dat.m

# ==============================================================================
# Section 2: DLNM Construction with Step-Change Interaction
# ==============================================================================
# This function fits a quasi-Poisson regression model incorporating a cross-basis
# function for temperature and an interaction term for the post-2023 period.
fit_dlnm_step_change <- function(df, outcome, meteo.var, max.lag, df.var = 3, df.lag = 3, df.time = 8) {
  
  df <- add_time_factors(df)
  
  # Construct Cross-basis matrix (exposure-lag-response)
  cb <- crossbasis(
    df[[meteo.var]],
    lag = max.lag,
    argvar = list(fun = 'poly', degree = 2), 
    arglag = list(fun = 'ns', df = df.lag)
  )
  colnames(cb) <- paste0('cb_', seq_len(ncol(cb)))
  
  # Construct Interaction term (Cross-basis * Post-Period Indicator)
  cb_post <- cb * df$post2023
  colnames(cb_post) <- paste0('cb_post_', seq_len(ncol(cb_post)))
  
  df2 <- bind_cols(df, as.data.frame(cb), as.data.frame(cb_post))
  has_dow <- length(unique(df2$dow[!is.na(df2$dow)])) >= 2
  
  cb_terms      <- paste(colnames(cb), collapse = ' + ')
  cb_post_terms <- paste(colnames(cb_post), collapse = ' + ')
  
  # Define potential confounders
  covars <- c('PM2.5', 'PM10', 'SO2', 'NO2', 'O3', 'CO', 'humidity', 'pressure', 'dew.temp', 'sunshine', 'wind.speed', 'precipitation', 'dif.temp')
  covars_use <- intersect(covars, names(df))
  covars_str <- paste(covars_use, collapse = ' + ')
  
  # Model A: Full Interaction Model
  rhs_full <- paste0(
    cb_terms, ' + ', cb_post_terms,
    ' + post2023 + ns(t, df=', df.time, ')',
    if (has_dow) ' + dow' else '',
    if (covars_str != "") paste(' + ', covars_str) else ''
  )
  fmla_full <- as.formula(paste0(outcome, ' ~ ', rhs_full))
  m_full <- glm(fmla_full, family = quasipoisson(link = 'log'), data = df2)
  
  # Model B: Null Model (No Interaction) for Likelihood Ratio Test
  rhs_null <- paste0(
    cb_terms, 
    ' + post2023 + ns(t, df=', df.time, ')',
    if (has_dow) ' + dow' else '',
    if (covars_str != "") paste(' + ', covars_str) else ''
  )
  fmla_null <- as.formula(paste0(outcome, ' ~ ', rhs_null))
  m_null <- glm(fmla_null, family = quasipoisson(link = 'log'), data = df2)
  
  # Statistical Testing (Wald Test for interaction terms & LRT)
  wald_res <- wald_test_block(m_full, 'cb_post_')
  lrt_res <- lrtest(m_null, m_full)
  
  list(
    m1 = m_full,
    m_null = m_null,
    cb = cb,          
    cb_post_mat = cb_post, 
    df_used = df2,
    meteo.var = meteo.var,
    max.lag = max.lag,
    wald_test = wald_res,
    lrt_test = lrt_res
  )
}

# Helper: Extract coefficients and covariance matrix for specific blocks
get_block <- function(model, prefix) {
  cf <- coef(model); V <- vcov(model)
  idx <- grep(paste0('^', prefix), names(cf))
  list(coef = cf[idx], vcov = V[idx, idx, drop = FALSE], names = names(cf)[idx])
}

# Helper: Perform Wald test on a block of coefficients
wald_test_block <- function(model, term_prefix = 'cb_post_') {
  cf <- coef(model); V <- vcov(model)
  idx <- grep(paste0('^', term_prefix), names(cf))
  b <- cf[idx]; VV <- V[idx, idx, drop = FALSE]
  if (length(b) == 0) return(data.frame(chi_sq = NA, df = 0, p_value = NA))
  stat <- as.numeric(t(b) %*% solve(VV, b))
  df <- length(b)
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  data.frame(chi_sq = stat, df = df, p_value = p)
}

# ==============================================================================
# Section 3: Extraction of Cumulative Risk Estimates
# ==============================================================================
extract_cum_effect_pre_post <- function(fit, df_original, contrast = c(0.10, 0.90), cen_quantile = 0.50) {
  
  # Define reference and exposure temperatures based on quantiles
  x <- df_original[[fit$meteo.var]]
  q_cold <- as.numeric(quantile(x, contrast[1], na.rm = TRUE))  # Cold exposure (Risk)
  q_warm <- as.numeric(quantile(x, contrast[2], na.rm = TRUE))  # Warm reference
  cen    <- as.numeric(quantile(x, cen_quantile, na.rm = TRUE)) # Centering value
  
  model <- fit$m1
  cf <- coef(model)
  vc <- vcov(model)
  
  # Map coefficients to Main (Pre) and Interaction (Difference) terms
  idx_main <- match(colnames(fit$cb), names(cf))
  idx_int  <- match(colnames(fit$cb_post_mat), names(cf))
  
  # 1. Pre-Pandemic Parameters (Main Effect only)
  coef_pre <- cf[idx_main]
  vcov_pre <- vc[idx_main, idx_main, drop = FALSE]
  
  # 2. Post-Pandemic Parameters (Main + Interaction)
  # Beta_post = Beta_main + Beta_interaction
  coef_post <- cf[idx_main] + cf[idx_int]
  
  # Var(Post) = Var(Main) + Var(Int) + 2*Cov(Main, Int)
  v_aa <- vc[idx_main, idx_main]
  v_bb <- vc[idx_int, idx_int]
  v_ab <- vc[idx_main, idx_int]
  vcov_post <- v_aa + v_bb + v_ab + t(v_ab)
  
  # 3. Compute Predictions (Relative Risks)
  # Note: Centering at q_warm ensures RR is relative to the warm reference
  
  # A. Pre-Pandemic
  pred_pre <- crosspred(fit$cb, coef = coef_pre, vcov = vcov_pre, model.link = 'log', cen = q_warm, at = q_cold, cumul = TRUE)
  logRR_pre    <- pred_pre$cumfit[1]
  logRR_pre_se <- pred_pre$cumse[1]
  
  # B. Post-Pandemic
  pred_post <- crosspred(fit$cb, coef = coef_post, vcov = vcov_post, model.link = 'log', cen = q_warm, at = q_cold, cumul = TRUE)
  logRR_post    <- pred_post$cumfit[1]
  logRR_post_se <- pred_post$cumse[1]
  
  # C. Interaction Effect (Ratio of RRs)
  coef_int <- cf[idx_int]
  vcov_int <- vc[idx_int, idx_int]
  pred_int <- crosspred(fit$cb, coef = coef_int, vcov = vcov_int, model.link = 'log', cen = q_warm, at = q_cold, cumul = TRUE)
  logRR_diff    <- pred_int$cumfit[1]
  logRR_diff_se <- pred_int$cumse[1]
  
  # Compile results
  data.frame(
    q_cold = q_cold, 
    q_warm = q_warm, 
    cen = cen,
    
    # Pre-Pandemic Estimates
    RR_pre = exp(logRR_pre),
    RR_pre_low = exp(logRR_pre - 1.96 * logRR_pre_se),
    RR_pre_high = exp(logRR_pre + 1.96 * logRR_pre_se),
    
    # Post-Pandemic Estimates
    RR_post = exp(logRR_post),
    RR_post_low = exp(logRR_post - 1.96 * logRR_post_se),
    RR_post_high = exp(logRR_post + 1.96 * logRR_post_se),
    
    # Comparison (Ratio of RRs)
    RR_ratio_post_vs_pre = exp(logRR_diff),
    RR_ratio_post_vs_pre_low = exp(logRR_diff - 1.96 * logRR_diff_se),
    RR_ratio_post_vs_pre_high = exp(logRR_diff + 1.96 * logRR_diff_se),
    
    # Significance Tests
    Wald_P_Value = fit$wald_test$p_value,
    LRT_P_Value = fit$lrt_test$`Pr(>Chisq)`[2]
  )
}

# ==============================================================================
# Section 4: Counterfactual Out-of-Sample Prediction
# ==============================================================================
perform_oos_prediction_fixed <- function(df, cut.date, meteo.var, max.lag, df.var, df.lag, target_col) {
  
  df <- add_time_factors(df)
  if (all(is.na(df$dow))) df$dow <- as.factor(format(df$date, '%u'))
  
  # Generate basis for the full dataset
  cb_all <- crossbasis(
    df[[meteo.var]], lag = max.lag,
    argvar = list(fun = 'poly', degree = 2), 
    arglag = list(fun = 'ns', df = df.lag)
  )
  colnames(cb_all) <- paste0('cb_', seq_len(ncol(cb_all)))
  
  df_mod <- bind_cols(df, as.data.frame(cb_all))
  
  # Define Training Set (Pre-COVID only)
  train_idx <- which(df_mod$Period == 'Pre-COVID')
  
  cb_terms <- paste(colnames(cb_all), collapse = ' + ')
  
  # Dynamic covariate selection based on data availability
  potential_covars <- c('PM2.5', 'humidity', 'dow')
  covars_use <- intersect(potential_covars, names(df_mod))
  final_covars <- c()
  for (cv in covars_use) {
    vec <- df_mod[[cv]][train_idx]
    vec_clean <- na.omit(vec)
    if (length(vec_clean) == 0) next 
    if (is.factor(vec_clean) || is.character(vec_clean)) {
      if (length(unique(vec_clean)) >= 2) final_covars <- c(final_covars, cv)
    } else {
      final_covars <- c(final_covars, cv)
    }
  }
  covars_str <- paste(final_covars, collapse = ' + ')
  
  # Model specification: Cross-basis + Seasonality (yday) + Confounders
  rhs <- paste0(cb_terms, ' + ns(yday, df=3)') 
  if (covars_str != "") rhs <- paste0(rhs, ' + ', covars_str)
  
  fmla_pred <- as.formula(paste0(target_col, ' ~ ', rhs))
  
  # Fit on training data
  m_pred <- glm(fmla_pred, family = quasipoisson(link = 'log'), data = df_mod, subset = train_idx, na.action = na.exclude)
  
  # Predict on test data (Post-2023)
  test_idx <- which(df_mod$date >= cut.date)
  pred_vals <- predict(m_pred, newdata = df_mod[test_idx, ], type = 'response')
  
  res_oos <- df_mod[test_idx, c('date', target_col)]
  colnames(res_oos)[2] <- 'Observed'
  res_oos$Predicted <- pred_vals
  res_oos$Virus <- target_col
  
  return(res_oos)
}

# ==============================================================================
# Section 5: Main Analysis Pipeline (Influenza A & B)
# ==============================================================================
all_effects_list <- list()
all_curves_list <- list()
all_oos_list <- list()

cat('\n[START] Initiating analysis loop for targets: ', paste(target.list, collapse = ', '), '\n')

for (target in target.list) {
  
  lbl <- ifelse(target == 'FLUA', yes = 'Influneza A', 'Influenza B')
  cat(glue('\n>>> Processing Target: {target} <<<\n\n'))
  
  # 1. Model Fitting
  fit_obj <- fit_dlnm_step_change(dat.m, target, meteo.main, max.lag, df.var, df.lag, df.time)
  
  # 2. Effect Extraction (Point estimates at cold quantile)
  eff_df <- extract_cum_effect_pre_post(fit_obj, dat.m) %>% mutate(Virus = target)
  all_effects_list[[target]] <- eff_df
  
  # 3. Counterfactual Prediction (OOS)
  oos_res <- perform_oos_prediction_fixed(dat.m, cut.date, meteo.main, max.lag, df.var, df.lag, target)
  # Scale to percentage for visualization
  oos_res$Observed <- oos_res$Observed * 100
  oos_res$Predicted <- oos_res$Predicted * 100
  all_oos_list[[target]] <- oos_res
  
  # Visualization: Observed vs Counterfactual
  p_oos <- ggplot(oos_res, aes(x = date)) +
    geom_line(aes(y = Observed, color = 'Observed'), linewidth = 0.8, alpha = 0.8) +
    geom_line(aes(y = Predicted, color = 'Counterfactual'), linewidth = 0.8, linetype = 'dashed') +
    scale_color_manual(values = c('Observed' = '#E50914', 'Counterfactual' = '#00A087')) +
    scale_y_continuous(limits = c(0, max(oos_res$Predicted, oos_res$Observed) * 1.2), expand = c(0, 0)) +
    scale_x_date(
      date_breaks = "3 months",
      date_labels = "%b\n%Y",
      expand = c(0, 0)
    ) +
    labs(
      subtitle = glue('Observed vs Predicted ({lbl})'),
      y = 'Positive Rate (%)', x = NULL, color = NULL
    ) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      legend.position = 'inside',
      legend.position.inside = c(0.8, 0.9),
      legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
      plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = 'black', size = base.size * 0.9, margin = margin(t = 5), family = base.family),
      axis.text.y = element_text(color = 'black', size = base.size * 0.9, family = base.family),
      axis.title.y = element_text(size = base.size * 0.9, face = 'bold', margin = margin(r = 10), family = base.family),
      axis.title.x = element_blank(),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6, color = 'black')
    ); p_oos
  
  ggsave(glue('{fig.dir}/{ofig}_{meteo.main}_oos_{target}_{this}.pdf'), p_oos, width = 6, height = 3)
}

# ==============================================================================
# Section 6: Visualization of Risk Shifts (Dumbbell Plot)
# ==============================================================================
final_eff_df <- bind_rows(all_effects_list)

# Export summary statistics
print(final_eff_df)
write.csv(final_eff_df, glue('{tbl.dir}/{ofig}_{meteo.main}_effects_{this}.csv'), row.names = FALSE)

# Generate Dumbbell Plot: Comparison of RR at cold temperatures
p_dumbbell_all <- ggplot(final_eff_df) +
  geom_segment(aes(x = RR_pre, xend = RR_post, y = Virus, yend = Virus), color = 'grey70', size = 1.5) +
  geom_point(aes(x = RR_pre, y = Virus, color = 'Pre-Pandemic'), size = 5) +
  geom_point(aes(x = RR_post, y = Virus, color = 'Post-Pandemic'), size = 5) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'black', alpha = 0.3) +
  geom_text(aes(x = RR_pre, y = Virus, label = sprintf('%.3f', RR_pre)), vjust = 2.5, color = '#00A087', fontface = 'bold', size = base.size / 3.88,  family = base.family) +
  geom_text(aes(x = RR_post, y = Virus, label = sprintf('%.3f', RR_post)), vjust = 2.5, color = '#E50914', fontface = 'bold', size = base.size / 3.88, family = base.family) +
  scale_color_manual(values = c('Pre-Pandemic' = '#00A087', 'Post-Pandemic' = '#E50914')) +
  scale_y_discrete(labels = c('FLUA' = 'Influenza A', 'FLUB' = 'Influenza B')) +
  scale_x_continuous(expand = c(0.05, 0.02)) +
  labs(
    subtitle = 'Comparison of risk sensitivity at cold temperatures',
    x = 'Relative Risk (RR)', y = NULL, color = 'Period'
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(ifelse(this == 'HK', 0.80, 0.25), 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = 'black', size = base.size * 0.9, margin = margin(t = 5), family = base.family),
    axis.text.y = element_text(color = 'black', size = base.size * 0.9, family = base.family, angle = 90, hjust = 0.5),
    axis.title.y = element_text(size = base.size * 0.9, face = 'bold', margin = margin(r = 10), family = base.family),
    axis.title.x = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6, color = 'black')
  ); p_dumbbell_all

ggsave(glue('{fig.dir}/{ofig}_{meteo.main}_dumbbell_{this}.pdf'), p_dumbbell_all, width = 4.5, height = 4)

# ==============================================================================
# Section 7: Sensitivity Analysis (Lag Structure)
# ==============================================================================
# Function to calculate QAIC for model selection
get_qaic <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, lambda = fitted(model), log = TRUE))
  k <- length(coef(model))
  return(-2 * (loglik / phi) + 2 * k)
}

lags_to_test <- c(3, 4, 5, 6, 7, 8, 9, 10, 11)
sens_res_list <- list()
qaic_res_list <- list()

cat('\n[SENSITIVITY] Running lag sensitivity analysis for FLUA and FLUB...\n')

for (target in target.list) {
  cat(glue('  > Target: {target} ...\n'))
  
  for (lg in lags_to_test) {
    # Fit model with varying max lags
    fit_sens <- fit_dlnm_step_change(dat.m, target, meteo.main, max.lag = lg, df.var, df.lag, df.time)
    
    # Store QAIC
    qaic_res_list[[paste(target, lg, sep = '_')]] <- data.frame(
      Virus = target,
      Lag = lg,
      QAIC = get_qaic(fit_sens$m1)
    )
    
    # Store Effect Estimates
    sens_res_list[[paste(target, lg, sep = '_')]] <- extract_cum_effect_pre_post(fit_sens, dat.m) %>% 
      mutate(Virus = target, Lag = lg)
  }
}

df_qaic_all <- bind_rows(qaic_res_list)
df_sens_all <- bind_rows(sens_res_list)
if (FALSE) write.csv(df_qaic_all, glue('{tbl.dir}/{ofig}_{meteo.main}_QAIC_sensitivity_{this}.csv'), row.names = FALSE)

# Visualization: Robustness of RR Ratio across lags
p_sens_all <- ggplot(df_sens_all, aes(x = factor(Lag), y = RR_ratio_post_vs_pre, group = Virus, color = Virus)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c('FLUA' = '#E50914', 'FLUB' = '#00A087'), labels = c('FLUA' = 'Influenza A', 'FLUB' = 'Influenza B')) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2), labels = c('0', '0.5', '1.0', '1.5', '2.0'), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  labs(
    title = ,
    subtitle = 'Robustness of RR ratio',
    x = 'Max Lag (Weeks)',
    y = 'RR Ratio (Post / Pre)'
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.8),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    legend.title = element_blank(),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = 'black', size = base.size * 0.9, margin = margin(t = 5), family = base.family),
    axis.text.y = element_text(color = 'black', size = base.size * 0.9, family = base.family, angle = 90, hjust = 0.5),
    axis.title.y = element_text(size = base.size * 0.9, face = 'bold', margin = margin(r = 10), family = base.family),
    axis.title.x = element_blank(),
    axis.line = element_line(size = 0.6),
    axis.ticks = element_line(size = 0.6, color = 'black')
  ); p_sens_all

ggsave(glue('{fig.dir}/{ofig}_{meteo.main}_sensitivity_{this}.pdf'), p_sens_all, width = 4, height = 4)