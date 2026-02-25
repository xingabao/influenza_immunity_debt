# Load R packages
suppressMessages(suppressWarnings(library(zoo)))
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(rEDM)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(future)))
suppressMessages(suppressWarnings(library(lmtest)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(ISOweek)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(RTransferEntropy)))
Sys.setlocale('LC_TIME', 'C')
start.time <- Sys.time()

# Define a function
process_gisaid_metadata <- function(metafile) {
  
  # Check if file exists
  if (!file.exists(metafile)) stop("File not found: ", metafile)
  
  # Read Data
  raw_meta <- suppressWarnings(readxl::read_excel(metafile))
  
  # Ensure required columns exist (flexible matching)
  req_cols <- c("Collection_Date", "Subtype", "Lineage")
  if (!all(req_cols %in% colnames(raw_meta))) {
    stop("Missing required columns: Collection_Date, Subtype, or Lineage")
  }
  
  # Process Data
  final_dat <- raw_meta %>%
    dplyr::select(`Collection date` = Collection_Date, Type = Subtype, Lineage = Lineage) %>%
    filter(str_length(`Collection date`) == 10) %>% 
    mutate(date = as.Date(`Collection date`)) %>%
    filter(!is.na(date)) %>%
    mutate(
      Virus_Type = case_when(
        str_detect(Type, "A") | str_detect(Lineage, "H1N1|H3N2|pdm09") ~ "FLUA",
        str_detect(Type, "B") | str_detect(Lineage, "Victoria|Yamagata") ~ "FLUB",
        TRUE ~ "Other"
      )
    ) %>%
    filter(Virus_Type %in% c("FLUA", "FLUB")) %>%
    mutate(WeekDate = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(WeekDate, Virus_Type) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(
      names_from = Virus_Type, 
      values_from = Count, 
      values_fill = list(Count = 0)
    ) %>%
    rename(date = WeekDate) %>%
    arrange(date)
  
  final_dat$date <- final_dat$date - 1
  
  # Complete time series
  if (nrow(final_dat) > 0) {
    full_grid <- data.frame(date = seq(min(final_dat$date), max(final_dat$date), by = "week"))
    
    final_dat <- full_grid %>%
      left_join(final_dat, by = "date") %>%
      replace_na(list(FLUA = 0, FLUB = 0)) %>%
      filter(date >= as.Date('2010-01-01'))
  }
  
  return(final_dat)
}

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
seed <- 100
this <- 'HK'  # Options: 'HK' or 'Macau'

# Plotting Constants
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Load data
set.seed(seed)
dat.cli <- readxl::read_excel(glue('{dat.dir}/{this}/Climate.xlsx'))

# Prepare plot data
if (this == 'Macau') {
  dat.fac <- readxl::read_excel(glue('{dat.dir}/{this}/influenza.xlsx'))
  dat.pl <- dat.fac %>% dplyr::select(date, FLUA, FLUB) %>%
    left_join(dat.cli %>% dplyr::select(date, ave.temp, humidity)) %>%
    mutate(Year = year(date), Week = week(date)) %>%
    filter(Year > 2014) %>%
    arrange(date) %>%
    mutate(TimeIndex = row_number())
} else {
  dat.fac <- process_gisaid_metadata(glue('{dat.dir}/{this}/gisaid_epiflu_isolates.xls')) %>%
    filter(FLUA > 0 | FLUB > 0)
  dat.pl <- dat.fac %>%
    left_join(dat.cli %>% dplyr::select(date, ave.temp, humidity)) %>%
    mutate(Year = year(date), Week = week(date)) %>%
    arrange(date) %>%
    mutate(TimeIndex = row_number())
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove seasonality (cyclic trends) from time series data and standardize residuals
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
remove_seasonality <- function(x, period = 52) {
  
  if (any(is.na(x))) x <- na.approx(x, rule = 2) 
  
  t <- 1:length(x)
  fit <- lm(x ~ sin(2 * pi * t / period) + cos(2 * pi * t / period))
  res <- residuals(fit)
  
  if (sd(res) < 1e-6) return(rep(0, length(res)))
  
  return(as.numeric(scale(res)))
}

dat.pl$resFLUA <- remove_seasonality(dat.pl$FLUA)
dat.pl$resFLUB <- remove_seasonality(dat.pl$FLUB)

dat.ccm <- dat.pl %>% 
  dplyr::select(TimeIndex, resFLUA, resFLUB) %>%
  filter(!is.na(resFLUA) & !is.na(resFLUB) & !is.infinite(resFLUA) & !is.infinite(resFLUB))

dat.ccm <- as.data.frame(dat.ccm)

cat(glue('Data rows after cleaning: {nrow(dat.ccm)}\n'))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Robust Optimal E Search
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
find_optimal_E <- function(df, col_name, max_E = 10) {
  
  rhos <- numeric(max_E)
  n <- nrow(df)
  
  lib_str <- paste('1', n)
  
  for (e in 1:max_E) {
    out <- tryCatch({
      Simplex(dataFrame = df, lib = lib_str, pred = lib_str, E = e, columns = col_name, target = col_name, showPlot = FALSE)
    }, error = function(err) { return(NULL) })
    
    if (!is.null(out) && nrow(out) > 0) {
      valid_idx <- complete.cases(out$Observations, out$Predictions)
      if (sum(valid_idx) > 10) {
        rhos[e] <- cor(out$Observations[valid_idx], out$Predictions[valid_idx])
      } else {
        rhos[e] <- NA
      }
    } else {
      rhos[e] <- NA
    }
  }
  
  best_e <- which.max(rhos)
  if (length(best_e) == 0) return(2)
  return(best_e)
}

E_A <- find_optimal_E(dat.ccm, 'resFLUA')
E_B <- find_optimal_E(dat.ccm, 'resFLUB')
best_E <- max(E_A, E_B, na.rm = TRUE)
cat(glue('Optimal E selected: {best_E}\n'))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sensitivity Analysis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('\nRunning Sensitivity Analysis for E...\n')
E_range <- (best_E - 2):(best_E + 2)
E_range <- E_range[E_range > 1] 
sens_results <- data.frame()
max_lib <- paste('1', nrow(dat.ccm)) 

for (e_val in E_range) {
  # A -> B
  ccm_AB <- tryCatch(
    CCM(
      dataFrame = dat.ccm, 
      E = e_val, 
      columns = 'resFLUB',
      target = 'resFLUA',
      libSizes = max_lib,
      sample = 1,
      random = FALSE, 
      showPlot = FALSE
    ), 
    error = function(e) NULL
  )
  rho_AB <- if (!is.null(ccm_AB)) tail(ccm_AB$`resFLUA:resFLUB`, 1) else NA
  
  # B -> A
  ccm_BA <- tryCatch(
    CCM(
      dataFrame = dat.ccm, 
      E = e_val, 
      columns = 'resFLUA',
      target = 'resFLUB',
      libSizes = max_lib,
      sample = 1, 
      random = FALSE,
      showPlot = FALSE
    ), error = function(e) NULL
  )
  rho_BA <- if (!is.null(ccm_BA)) tail(ccm_BA$`resFLUB:resFLUA`, 1) else NA
  
  sens_results <- rbind(sens_results, data.frame(E = e_val, Direction = 'Flu A -> Flu B', Rho = rho_AB))
  sens_results <- rbind(sens_results, data.frame(E = e_val, Direction = 'Flu B -> Flu A', Rho = rho_BA))
}

p_sens <- ggplot(sens_results, aes(x = factor(E), y = Rho, fill = Direction)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
  scale_fill_manual(values = c('Flu A -> Flu B' = '#E50914', 'Flu B -> Flu A' = '#00A087')) +
  scale_y_continuous(limits = c(0, round(max(sens_results$Rho) * 1.2, 1)), expand = c(0, 0)) +
  labs(
    subtitle = 'Sensitivity to E', 
    x = 'Embedding Dimension (E)', 
    y = 'Cross Map Skill (Rho)',
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5))
  ); p_sens

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CCM
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
max_len <- nrow(dat.ccm)
lib_str <- paste(seq(10, max_len - 5, by = 25), collapse = ' ')
nsample <- 100
n_surrogates <- 500 

run_ccm_with_surrogates <- function(df, lib_col, tar_col, E, lib_str, n_samp, n_surr) {
  
  col_name_match <- paste0(tar_col, ':', lib_col)
  
  # 1.Real Data
  ccm_real <- tryCatch({
    CCM(
      dataFrame = df, 
      E = E, 
      columns = lib_col, 
      target = tar_col,
      libSizes = lib_str, 
      sample = n_samp, 
      random = TRUE,
      showPlot = FALSE
    )
  }, error = function(e) {
    message('Error in Real Data CCM: ', e$message)
    return(NULL)
  })
  
  if (is.null(ccm_real)) return(NULL)
  
  res_real <- ccm_real %>%
    dplyr::select(LibSize, Rho = all_of(col_name_match)) %>%
    group_by(LibSize) %>%
    summarise(Rho = mean(Rho, na.rm = TRUE), .groups = 'drop') %>%
    mutate(Type = 'Real Data')
  
  # 2.Surrogates (Parallel)
  surr_mat <- SurrogateData(df[[tar_col]], method = 'ebisuzaki', num_surr = n_surr)
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {library(rEDM); library(dplyr)})
  clusterExport(cl, varlist = c('df', 'tar_col', 'lib_col', 'E', 'lib_str', 'col_name_match', 'surr_mat'), envir = environment())
  
  surr_rhos_list <- parLapply(cl, 1:n_surr, function(i) {
    df_surr <- df
    df_surr[[tar_col]] <- surr_mat[, i]
    tryCatch({
      ccm_out <- CCM(
        dataFrame = df_surr,
        E = E, 
        columns = lib_col, 
        target = tar_col,
        libSizes = lib_str,
        sample = 10,
        random = TRUE, 
        showPlot = FALSE
      )
      ccm_out %>%
        dplyr::select(LibSize, Rho = all_of(col_name_match)) %>%
        group_by(LibSize) %>%
        summarise(Rho = mean(Rho, na.rm = TRUE), .groups = 'drop') %>%
        pull(Rho)
    }, error = function(e) {
      return(rep(NA, length(strsplit(lib_str, ' '))))
    })
  })
  stopCluster(cl)
  
  surr_combined <- do.call(cbind, surr_rhos_list)
  
  upper_95 <- apply(surr_combined, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  res_surr <- data.frame(LibSize = res_real$LibSize, Rho = upper_95, Type = '95% Null Limit')
  
  # 3.Convergence Test
  tau_test <- cor.test(res_real$LibSize, res_real$Rho, method = 'kendall', alternative = 'greater')
  
  return(list(real = res_real, limit = res_surr, p_val = tau_test$p.value))
}

cat('\nRunning Main CCM Analysis...\n')
res_A_cause_B <- suppressWarnings(run_ccm_with_surrogates(dat.ccm, 'resFLUB', 'resFLUA', best_E, lib_str, nsample, n_surrogates))
res_B_cause_A <- suppressWarnings(run_ccm_with_surrogates(dat.ccm, 'resFLUA', 'resFLUB', best_E, lib_str, nsample, n_surrogates))

if (!is.null(res_A_cause_B) && !is.null(res_B_cause_A)) {
  dat.plot.ccm <- bind_rows(
    res_A_cause_B$real %>% mutate(Direction = 'Flu A causes Flu B', Class = 'Observed'),
    res_A_cause_B$limit %>% mutate(Direction = 'Flu A causes Flu B', Class = 'Null 95%'),
    res_B_cause_A$real %>% mutate(Direction = 'Flu B causes Flu A', Class = 'Observed'),
    res_B_cause_A$limit %>% mutate(Direction = 'Flu B causes Flu A', Class = 'Null 95%')
  )
  
  # Plot CCM
  p_ccm <- ggplot(dat.plot.ccm, aes(x = LibSize, y = Rho)) +
    geom_line(data = filter(dat.plot.ccm, Class == 'Null 95%'), aes(linetype = Class), color = 'grey50', linewidth = 0.8) +
    geom_line(data = filter(dat.plot.ccm, Class == 'Observed'), aes(color = Direction), linewidth = 1.2) +
    facet_wrap(~Direction) +
    scale_color_manual(values = c('Flu A causes Flu B' = '#E50914', 'Flu B causes Flu A' = '#00A087')) +
    scale_linetype_manual(values = c('Null 95%' = 'dashed'), name = '') +
    scale_y_continuous(limits = c(0, round(max(dat.plot.ccm$Rho) * 1.2, 1)), expand = c(0, 0)) +
    labs(subtitle = 'CCM Causality Test (Deseasonalized)', y = 'Cross Map Skill (Rho)', x = 'Library Size (L)') +
    guides(color = 'none') +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      legend.position = 'inside', 
      legend.position.inside = c(0.80, 0.25),
      strip.background = element_blank(),
      strip.text = element_text(face = 'bold', size = base.size * 0.9, family = base.family, color = base.col),
      axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
      axis.title.y = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
      axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5)),
      plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col)
    ); p_ccm
} else {
  cat('CCM Analysis failed due to data issues.\n')
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Linear Granger (Linear Baseline)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('\nRunning Linear Granger Test (lmtest)...\n')

# Test 1: A -> B
g_A_to_B <- grangertest(resFLUB ~ resFLUA, order = best_E, data = dat.ccm)

# Test 2: B -> A
g_B_to_A <- grangertest(resFLUA ~ resFLUB, order = best_E, data = dat.ccm)

granger_res <- data.frame(
  Direction = c('Flu A -> Flu B', 'Flu B -> Flu A'),
  Stat = c(g_A_to_B$F[2], g_B_to_A$F[2]),
  P_val = c(g_A_to_B$`Pr(>F)`[2], g_B_to_A$`Pr(>F)`[2]),
  Type = 'Linear (F-test)'
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting Granger Causality Test
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
add_sig <- function(df) {
  df %>% mutate(Sig = case_when(
    P_val < 0.001 ~ '***',
    P_val < 0.01 ~ '**',
    P_val < 0.05 ~ '*',
    TRUE ~ 'NS'
  ))
}

granger_res <- add_sig(granger_res)
mm <- round(max(granger_res$Stat) * 1.2, 1)

# Linear Granger
p_granger <- ggplot(granger_res, aes(x = Direction, y = Stat, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.6, alpha = 0.8) +
  geom_text(aes(label = Sig, y = Stat), vjust = -0.5, size = base.size / 4.5, fontface = 'bold', family = base.family) +
  scale_fill_manual(values = c('Flu A -> Flu B' = '#E50914', 'Flu B -> Flu A' = '#00A087')) +
  scale_x_discrete(labels = c('A to B', 'B to A')) +
  scale_y_continuous(limits = c(0, mm), expand = c(0, 0)) +
  labs(subtitle = 'Linear Granger\n  Causality Test', y = 'F-Score', x = 'Causality Direction') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'none',
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5))
  ); p_granger

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Non-linear Granger
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('\nRunning Non-linear Granger Test (lmtest)...\n')

model_data <- dat.pl %>%
  arrange(date) %>%
  mutate(
    FLUA_lag1 = lag(FLUA, 1),
    FLUA_lag2 = lag(FLUA, 2),
    FLUB_lag1 = lag(FLUB, 1),
    FLUB_lag2 = lag(FLUB, 2),
    temp_lag1 = lag(ave.temp, 1),
    humid_lag1 = lag(humidity, 1)
  ) %>%
  drop_na()

# Granger causality test
run_gam_granger <- function(data, target, cause, lags = 1:2) {
  
  base_terms <- paste0(
    paste0("s(", target, "_lag", lags, ")", collapse = " + "), 
    " + s(ave.temp) + s(humidity) + s(Week, bs = 'cc', k = 52) + s(TimeIndex)"
  )
  
  cause_terms <- paste0(paste0("s(", cause, "_lag", lags, ")", collapse = " + "))
  
  f_null <- as.formula(paste(target, "~", base_terms))
  f_full <- as.formula(paste(target, "~", base_terms, "+", cause_terms))
  
  model_null <- gam(f_null, data = data, method = "REML", family = tw(link = "log"))
  model_full <- gam(f_full, data = data, method = "REML", family = tw(link = "log"))
  
  test_res <- anova(model_null, model_full, test = "F")

  p_value <- test_res$`Pr(>Chi)`[2]
  dev_diff <- test_res$Deviance[2]
  is_sig <- ifelse(p_value < 0.05, "Significant", "Not Significant")
  
  print(test_res)
  message(sprintf("\nConclusion: Does %s cause %s? -> %s", cause, target, is_sig))
  message(sprintf(
    "P-value: %.6f | Deviance Improvement: %.4f", 
    ifelse(is.na(p_value), NA, p_value), 
    ifelse(is.na(dev_diff), NA, dev_diff)
  ))
  
  return(list(null = model_null, full = model_full, anova = test_res))
}

# Perform bidirectional testing
res_A_to_B <- run_gam_granger(model_data, target = "FLUB", cause = "FLUA")
res_B_to_A <- run_gam_granger(model_data, target = "FLUA", cause = "FLUB")

ngranger_res <- data.frame(
  Direction = c('Flu A -> Flu B', 'Flu B -> Flu A'),
  Stat = c(res_A_to_B$anova$Deviance[2], res_B_to_A$anova$Deviance[2]),
  P_val = c(res_A_to_B$anova$`Pr(>Chi)`[2], res_B_to_A$anova$`Pr(>Chi)`[2]),
  Type = 'Non-linear (Deviance)'
)

ngranger_res <- add_sig(ngranger_res)
mm <- round(max(ngranger_res$Stat) * 1.2, 1)

# Non-linear Granger
p_ngranger <- ggplot(ngranger_res, aes(x = Direction, y = Stat, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.6, alpha = 0.8) +
  geom_text(aes(label = Sig, y = Stat), vjust = -0.5, size = base.size / 4.5, fontface = 'bold', family = base.family) +
  scale_fill_manual(values = c('Flu A -> Flu B' = '#E50914', 'Flu B -> Flu A' = '#00A087')) +
  scale_x_discrete(labels = c('A to B', 'B to A')) +
  scale_y_continuous(limits = c(0, mm), expand = c(0, 0)) +
  labs(subtitle = 'Non-linear Granger\n  Causality Test', y = 'Deviance', x = 'Causality Direction') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'none',
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5))
  ); p_ngranger

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DLNM Analysis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('\nRunning Refined DLNM Analysis...\n')

lag.w <- 10

cb <- crossbasis(
  dat.pl$FLUA, 
  lag = lag.w, 
  argvar = list(fun = 'ns', df = 3),
  arglag = list(fun = 'ns', df = 3)
)

model.dlnm <- glm(
  FLUB ~ cb + 
    ns(ave.temp, df = 3) + 
    ns(humidity, df = 3) + 
    ns(TimeIndex, df = 3) + 
    ns(Week, df = 3), 
  family = quasipoisson(), 
  data = dat.pl
)

exposure_val <- quantile(dat.pl$FLUA, 0.75, na.rm = TRUE)

pred <- crosspred(cb, model.dlnm, at = exposure_val, cen = 0) 

dat.dlnm <- data.frame(
  Lag = 0:lag.w, 
  RR = as.numeric(pred$matRRfit[1, ]), 
  Lower = as.numeric(pred$matRRlow[1, ]), 
  Upper = as.numeric(pred$matRRhigh[1, ])
)

# Plot
p_dlnm <- ggplot(dat.dlnm, aes(x = Lag, y = RR)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40') +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = '#E74C3C', alpha = 0.2) +
  geom_line(color = '#C0392B', linewidth = 1) +
  scale_x_continuous(breaks = 0:lag.w) +
  labs(
    subtitle = glue('Lagged Effect of Flu A on Flu B'), 
    x = 'Lag (Weeks)', 
    y = 'Relative Risk (RR)'
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5))
  )

print(p_dlnm)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (exists('p_ccm') && exists('p_dlnm') && exists('p_sens') && exists('p_granger')) {
  
  # Layout Design:
  design <- '
  AABBB
  CCEEE
  '
  gg <- p_sens + p_dlnm + (p_granger + p_ngranger) + p_ccm + 
    plot_layout(design = design, heights = c(0.8, 1))
  
  print(gg)
  
} else {
  cat('Some plots failed to generate.\n')
}

# Save to file
width = 12; height = 9
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')

end.time <- Sys.time()
run.time <- end.time - start.time
print(run.time)