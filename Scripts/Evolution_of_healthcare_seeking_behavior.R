# Load R packages
suppressMessages(suppressWarnings(library(zoo)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(rEDM)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggh4x)))
suppressMessages(suppressWarnings(library(lmtest)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(lubridate)))
Sys.setlocale('LC_TIME', 'C')
start.time <- Sys.time()

# -----------------------------
# Reproducibility
# -----------------------------
set.seed(20240101)
RNGkind("L'Ecuyer-CMRG")

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'
this <- 'HK'  # Options: 'HK' or 'Macau'
control <- TRUE

# Define color scheme
col_flu <- '#3C5488FF'
col_search <- '#E64B35FF'
col_cause <- '#00A087FF'
col_reverse <- '#8491B4FF'

# Define dates
sdate <- '2020-01-01'
edate <- '2023-01-01'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Deseasonalization：sin/cos + residual + scale
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
remove_seasonality <- function(x, period = 52) {
  if (any(is.na(x))) x <- na.approx(x, rule = 2)
  t <- 1:length(x)
  fit <- lm(x ~ sin(2 * pi * t / period) + cos(2 * pi * t / period))
  res <- residuals(fit)
  if (sd(res, na.rm = TRUE) < 1e-6) return(rep(0, length(res)))
  return(as.numeric(scale(res)))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Robust optimal E search (Simplex)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
find_optimal_E <- function(df, col_name, max_E = 10) {
  rhos <- numeric(max_E)
  n <- nrow(df)
  lib_str <- paste('1', n)
  
  for (e in 1:max_E) {
    out <- tryCatch({
      Simplex(
        dataFrame = df,
        lib = lib_str,
        pred = NULL,
        E = e,
        columns = col_name,
        target = col_name,
        showPlot = FALSE
      )
    }, error = function(err) NULL)
    
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
  if (length(best_e) == 0 || is.na(best_e) || best_e < 2) return(2)
  return(best_e)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CCM + Surrogates + Kendall convergence test
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
run_ccm_with_surrogates <- function(df, lib_col, tar_col, E, lib_str, n_samp, n_surr, smooth_null = TRUE, smooth_k = 3, cl = NULL) {
  
  col_name_match <- paste0(tar_col, ':', lib_col)
  
  # 1.Real data CCM
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
  
  # 2.Surrogates: ebisuzaki (phase-randomized) for target series
  surr_mat <- SurrogateData(df[[tar_col]], method = 'ebisuzaki', num_surr = n_surr)
  
  n_lib <- length(strsplit(lib_str, ' ')[[1]])
  
  # Surrogate loop: parallel if cluster provided, else serial
  if (!is.null(cl)) {
    # Ensure packages loaded on workers (safe to call repeatedly)
    clusterEvalQ(cl, {
      suppressMessages(suppressWarnings(library(rEDM)))
      suppressMessages(suppressWarnings(library(dplyr)))
    })
    
    # Export needed objects
    clusterExport(
      cl,
      varlist = c('df', 'tar_col', 'lib_col', 'E', 'lib_str', 'col_name_match', 'surr_mat', 'n_samp', 'n_lib'),
      envir = environment()
    )
    
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
          sample = n_samp,
          random = TRUE,
          showPlot = FALSE
        )
        ccm_out %>%
          dplyr::select(LibSize, Rho = all_of(col_name_match)) %>%
          group_by(LibSize) %>%
          summarise(Rho = mean(Rho, na.rm = TRUE), .groups = 'drop') %>%
          pull(Rho)
      }, error = function(e) {
        rep(NA_real_, n_lib)
      })
    })
    
  } else {
    # Serial
    surr_rhos_list <- lapply(1:n_surr, function(i) {
      df_surr <- df
      df_surr[[tar_col]] <- surr_mat[, i]
      
      tryCatch({
        ccm_out <- CCM(
          dataFrame = df_surr,
          E = E,
          columns = lib_col,
          target = tar_col,
          libSizes = lib_str,
          sample = n_samp,
          random = TRUE,
          showPlot = FALSE
        )
        ccm_out %>%
          dplyr::select(LibSize, Rho = all_of(col_name_match)) %>%
          group_by(LibSize) %>%
          summarise(Rho = mean(Rho, na.rm = TRUE), .groups = 'drop') %>%
          pull(Rho)
      }, error = function(e) {
        rep(NA_real_, n_lib)
      })
    })
  }
  
  surr_combined <- do.call(cbind, surr_rhos_list)
  upper_95 <- apply(surr_combined, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  
  if (smooth_null && length(upper_95) >= smooth_k && smooth_k >= 3) {
    upper_95 <- zoo::rollmean(upper_95, k = smooth_k, fill = NA, align = 'center')
    upper_95 <- zoo::na.approx(upper_95, rule = 2)
  }
  
  res_surr <- data.frame(LibSize = res_real$LibSize, Rho = upper_95, Type = '95% Null Limit')
  
  # 3.Convergence test: Kendall tau (increasing rho with LibSize)
  tau_test <- suppressWarnings(
    cor.test(res_real$LibSize, res_real$Rho, method = 'kendall', alternative = 'greater')
  )
  
  return(list(real = res_real, limit = res_surr, p_val = tau_test$p.value))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Non-linear Granger (GAM) with deviance test
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
run_gam_granger <- function(data, target, cause, lags = 1:2) {
  
  get_lag_prefix <- function(x) {
    x <- as.character(x)
    x <- sub("\\.DS$", "", x)
    x <- sub("\\.Z$",  "", x)
    x <- sub("\\.MA$", "", x)
    strsplit(x, "\\.")[[1]][1]
  }
  
  target_prefix <- get_lag_prefix(target)  # e.g. "Flu"
  cause_prefix  <- get_lag_prefix(cause)   # e.g. "Search"
  
  target_lag_cols <- paste0(target_prefix, "_lag", lags)
  cause_lag_cols  <- paste0(cause_prefix,  "_lag", lags)
  
  miss1 <- setdiff(c(target, cause), colnames(data))
  miss2 <- setdiff(c(target_lag_cols, cause_lag_cols), colnames(data))
  if (length(miss1) > 0) stop('Missing columns in data: ', paste(miss1, collapse = ', '))
  if (length(miss2) > 0) stop('Missing lag columns in data: ', paste(miss2, collapse = ', '))
  
  base_terms <- paste0(paste0("s(", target_lag_cols, ")", collapse = " + "), " + s(TimeIndex)")
  cause_terms <- paste0(paste0("s(", cause_lag_cols, ")", collapse = " + "))
  
  f_null <- as.formula(paste0("`", target, "` ~ ", base_terms))
  f_full <- as.formula(paste0("`", target, "` ~ ", base_terms, " + ", cause_terms))
  
  message(paste0("\n================ Detecting: ", cause, " -> ", target, " ================\n"))
  
  model_null <- mgcv::gam(f_null, data = data, method = "REML", family = gaussian())
  model_full <- mgcv::gam(f_full, data = data, method = "REML", family = gaussian())
  
  test_res <- anova(model_null, model_full, test = "F")
  
  p_value <- suppressWarnings(test_res$`Pr(>F)`[2])
  dev_diff <- suppressWarnings(test_res$Deviance[2])
  is_sig <- ifelse(!is.na(p_value) && p_value < 0.05, "Significant", "Not Significant")
  
  print(test_res)
  message(sprintf("\nConclusion: Does %s cause %s? -> %s", cause, target, is_sig))
  message(sprintf(
    "P-value: %s | Deviance Improvement: %s",
    ifelse(is.na(p_value), "NA", sprintf("%.6f", p_value)),
    ifelse(is.na(dev_diff), "NA", sprintf("%.4f", dev_diff))
  ))
  
  return(list(null = model_null, full = model_full, anova = test_res))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dat.raw <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx'))

if (this == 'Macau' && control) {
  dat.search <- data.frame()
  csvs <- list.files(glue('{rt.dir}/data/{this}/Search/Control'))
  for (csv in csvs) {
    csv.file <- glue('{rt.dir}/data/{this}/Search/Control/{csv}')
    dat.tmp <- readr::read_csv(csv.file, skip = 1, show_col_types = FALSE)
    colnames(dat.tmp) <- c('date', 'Search.index')
    dat.search <- rbind(dat.search, dat.tmp)
  }
} else {
  dat.search <- data.frame()
  csvs <- list.files(glue('{rt.dir}/data/{this}/Search'))
  for (csv in csvs) {
    csv.file <- glue('{rt.dir}/data/{this}/Search/{csv}')
    dat.tmp <- readr::read_csv(csv.file, skip = 1, show_col_types = FALSE)
    colnames(dat.tmp) <- c('date', 'Search.index')
    dat.search <- rbind(dat.search, dat.tmp)
  }
}

# Aggregation and merging
dat.search <- dat.search %>%
  mutate(Week = cut(date, 'week')) %>%
  group_by(Week) %>%
  summarise(Search.index = sum(Search.index, na.rm = TRUE), date = min(date), .groups = 'drop') %>%
  dplyr::select(-Week)

dat.merged <- dat.raw %>%
  dplyr::select(date, FLUAB) %>%
  arrange(date) %>%
  mutate(Week = cut(date, 'week')) %>%
  group_by(Week) %>%
  summarise(FLUAB = mean(FLUAB, na.rm = TRUE), date = min(date), .groups = 'drop') %>%
  dplyr::select(-Week) %>%
  mutate(Period = case_when(
    date < sdate ~ 'Pre-Pandemic',
    date > edate ~ 'Post-Pandemic',
    TRUE ~ 'Pandemic'
  )) %>%
  filter(Period %in% c('Pre-Pandemic', 'Post-Pandemic')) %>%
  left_join(dat.search, by = 'date') %>%
  na.omit()

dat.ana <- dat.merged %>%
  group_by(Period) %>%
  mutate(
    Search.MA = rollmean(Search.index, k = 3, fill = NA, align = 'center'),
    Flu.MA = rollmean(FLUAB, k = 3, fill = NA, align = 'center')
  ) %>%
  na.omit() %>%
  mutate(
    Search.Z = as.vector(scale(Search.MA)),
    Flu.Z = as.vector(scale(Flu.MA))
  ) %>%
  ungroup() %>%
  mutate(date = as.Date(date))

dat.ana$Period <- factor(dat.ana$Period, levels = c('Pre-Pandemic', 'Post-Pandemic'))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CCF -> Granger -> Non-linear Granger -> CCM (deseason + sens + surrogates)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
analyze_period_stats <- function(df_subset, period_name, ccm_sample = 100, lag_max_ccf = 10, max_E = 10, sens_span = 2, n_surrogates = 500, lib_step = 25, cl = NULL) {
  
  df_subset <- df_subset %>% arrange(date)
  
  ccf_res <- ccf(df_subset$Flu.MA, df_subset$Search.MA, lag.max = lag_max_ccf, plot = FALSE)
  ccf_df <- data.frame(Lag = ccf_res$lag, ACF = as.numeric(ccf_res$acf), Period = period_name)
  
  g_s_to_f <- grangertest(Flu.MA ~ Search.MA, order = 2, data = df_subset)
  g_f_to_s <- grangertest(Search.MA ~ Flu.MA, order = 2, data = df_subset)
  
  granger_df <- data.frame(
    Period = period_name,
    Direction = c('Search\nto\nInfluenza', 'Influenza\nto\nSearch'),
    F_stat = c(g_s_to_f$F[2], g_f_to_s$F[2]),
    P_val = c(g_s_to_f$`Pr(>F)`[2], g_f_to_s$`Pr(>F)`[2]),
    Type = 'Linear'
  )
  
  tmp <- df_subset %>%
    mutate(
      Flu.DS = remove_seasonality(Flu.MA, period = 52),
      Search.DS = remove_seasonality(Search.MA, period = 52),
      TimeIndex = row_number()
    )
  
  gam_data <- tmp %>%
    mutate(
      Flu_lag1 = lag(Flu.DS, 1),
      Flu_lag2 = lag(Flu.DS, 2),
      Search_lag1 = lag(Search.DS, 1),
      Search_lag2 = lag(Search.DS, 2)
    ) %>%
    drop_na()
  
  # GAM Granger: Search -> Flu
  res_ng_s_to_f <- tryCatch(run_gam_granger(gam_data, target = "Flu.DS", cause = "Search.DS", lags = 1:2), error = function(e) { message(e$message); NULL })
  
  # GAM Granger: Flu -> Search
  res_ng_f_to_s <- tryCatch(run_gam_granger(gam_data, target = "Search.DS", cause = "Flu.DS", lags = 1:2), error = function(e) { message(e$message); NULL })
  
  ngranger_df <- data.frame(
    Period = period_name,
    Direction = c('Search\nto\nInfluenza', 'Influenza\nto\nSearch'),
    Stat = c(
      if (!is.null(res_ng_s_to_f)) res_ng_s_to_f$anova$Deviance[2] else NA,
      if (!is.null(res_ng_f_to_s)) res_ng_f_to_s$anova$Deviance[2] else NA
    ),
    P_val = c(
      if (!is.null(res_ng_s_to_f)) res_ng_s_to_f$anova$`Pr(>F)`[2] else NA,
      if (!is.null(res_ng_f_to_s)) res_ng_f_to_s$anova$`Pr(>F)`[2] else NA
    ),
    Type = 'Non-linear (GAM)'
  )
  
  ccm_data <- tmp %>%
    dplyr::select(date, Flu.DS, Search.DS) %>%
    filter(!is.na(Flu.DS) & !is.na(Search.DS) & is.finite(Flu.DS) & is.finite(Search.DS)) %>%
    as.data.frame()
  
  if (nrow(ccm_data) < 30) {
    ccm_df <- data.frame()
    ccm_plot_df <- data.frame()
    sens_df <- data.frame()
    return(list(ccf = ccf_df, granger = granger_df, ngranger = ngranger_df, ccm = ccm_df, ccm_plot = ccm_plot_df, sens = sens_df))
  }
  
  # Find the optimal E
  E_flu <- find_optimal_E(ccm_data, 'Flu.DS', max_E = max_E)
  E_search <- find_optimal_E(ccm_data, 'Search.DS', max_E = max_E)
  if (is.na(E_flu) || E_flu < 2) E_flu <- 3
  if (is.na(E_search) || E_search < 2) E_search <- 3
  
  # Sensitivity analysis
  sens_results <- data.frame()
  max_len <- nrow(ccm_data)
  lib_vec <- seq(10, max_len - 5, by = lib_step)
  if (length(lib_vec) < 5) lib_vec <- unique(floor(seq(10, max_len - 5, length.out = 8)))
  lib_str <- paste(lib_vec, collapse = ' ')
  lib_str <- glue('{10} {max_len - 5} {lib_step}')
  E_range_flu <- (E_flu - sens_span):(E_flu + sens_span)
  E_range_flu <- E_range_flu[E_range_flu > 1]
  
  E_range_search <- (E_search - sens_span):(E_search + sens_span)
  E_range_search <- E_range_search[E_range_search > 1]
  
  # Search -> Flu
  for (e_val in E_range_flu) {
    print(glue('{e_val}, {period_name}, E_range_flu, {lib_str}'))
    ccm_s_to_f <- tryCatch(
      CCM(
        dataFrame = ccm_data,
        E = e_val,
        columns = 'Flu.DS',
        target = 'Search.DS',
        libSizes = lib_str,
        sample = 1,
        random = FALSE,
        showPlot = FALSE
      ),
      error = function(e) NULL
    )
    rho_s_to_f <- if (!is.null(ccm_s_to_f)) tail(ccm_s_to_f$`Search.DS:Flu.DS`, 1) else NA
    
    sens_results <- rbind(
      sens_results,
      data.frame(Period = period_name, E = e_val, Direction = 'Search to Influenza', Rho = rho_s_to_f)
    )
  }
  
  # Flu -> Search
  for (e_val in E_range_search) {
    print(glue('{e_val}, {period_name}, E_range_search, {lib_str}'))
    ccm_f_to_s <- tryCatch(
      CCM(
        dataFrame = ccm_data,
        E = e_val,
        columns = 'Search.DS',
        target = 'Flu.DS',
        libSizes = lib_str,
        sample = 1,
        random = FALSE,
        showPlot = FALSE
      ),
      error = function(e) NULL
    )
    rho_f_to_s <- if (!is.null(ccm_f_to_s)) tail(ccm_f_to_s$`Flu.DS:Search.DS`, 1) else NA
    
    sens_results <- rbind(
      sens_results,
      data.frame(Period = period_name, E = e_val, Direction = 'Influenza to Search', Rho = rho_f_to_s)
    )
  }
  
  # Main CCM with surrogates
  res_s_to_f <- suppressWarnings(run_ccm_with_surrogates(
    df = ccm_data,
    lib_col = 'Flu.DS',
    tar_col = 'Search.DS',
    E = E_flu,
    lib_str = lib_str,
    n_samp = ccm_sample,
    n_surr = n_surrogates,
    smooth_null = TRUE,
    smooth_k = 3,
    cl = cl
  ))
  
  # Flu -> Search (Flu causes Search) => Search predicts Flu, use E_search
  res_f_to_s <- suppressWarnings(run_ccm_with_surrogates(
    df = ccm_data,
    lib_col = 'Search.DS',
    tar_col = 'Flu.DS',
    E = E_search,
    lib_str = lib_str,
    n_samp = ccm_sample,
    n_surr = n_surrogates,
    smooth_null = TRUE,
    smooth_k = 3,
    cl = cl
  ))
  
  ccm_plot_df <- data.frame()
  ccm_df <- data.frame()
  
  if (!is.null(res_s_to_f) && !is.null(res_f_to_s)) {
    
    ccm_plot_df <- bind_rows(
      res_s_to_f$real  %>% mutate(Direction = 'Search to Influenza', Class = 'Observed', Period = period_name, E = E_flu,    P_kendall = res_s_to_f$p_val),
      res_s_to_f$limit %>% mutate(Direction = 'Search to Influenza', Class = 'Null 95%', Period = period_name, E = E_flu,    P_kendall = res_s_to_f$p_val),
      res_f_to_s$real  %>% mutate(Direction = 'Influenza to Search', Class = 'Observed', Period = period_name, E = E_search, P_kendall = res_f_to_s$p_val),
      res_f_to_s$limit %>% mutate(Direction = 'Influenza to Search', Class = 'Null 95%', Period = period_name, E = E_search, P_kendall = res_f_to_s$p_val)
    )
    
    last_rho_s_to_f <- res_s_to_f$real %>% arrange(LibSize) %>% tail(1) %>% pull(Rho)
    last_rho_f_to_s <- res_f_to_s$real %>% arrange(LibSize) %>% tail(1) %>% pull(Rho)
    
    ccm_df <- data.frame(
      Period = period_name,
      Direction = c('Search\nto\nInfluenza', 'Influenza\nto\nSearch'),
      Rho_last = c(last_rho_s_to_f, last_rho_f_to_s),
      E = c(E_flu, E_search),
      P_kendall = c(res_s_to_f$p_val, res_f_to_s$p_val)
    )
  }
  
  # -----------------------------
  # Ensure function returns results
  # -----------------------------
  return(list(
    ccf = ccf_df,
    granger = granger_df,
    ngranger = ngranger_df,
    ccm = ccm_df,
    ccm_plot = ccm_plot_df,
    sens = sens_results
  ))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Execute a loop
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
periods <- unique(dat.ana$Period)
ccf_combined <- data.frame()
granger_combined <- data.frame()
ngranger_combined <- data.frame()
ccm_summary <- data.frame()
ccm_plot_combined <- data.frame()
sens_combined <- data.frame()

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

for (p in periods) {
  sub_dat <- dat.ana %>% filter(Period == p)
  
  res <- analyze_period_stats(
    df_subset = sub_dat,
    period_name = p,
    ccm_sample = 100,
    lag_max_ccf = 10,
    max_E = 10,
    sens_span = 2,
    n_surrogates = 500,
    lib_step = 25,
    cl = cl
  )
  
  ccf_combined <- rbind(ccf_combined, res$ccf)
  granger_combined <- rbind(granger_combined, res$granger)
  ngranger_combined <- rbind(ngranger_combined, res$ngranger)
  ccm_summary <- rbind(ccm_summary, res$ccm)
  ccm_plot_combined <- rbind(ccm_plot_combined, res$ccm_plot)
  sens_combined <- rbind(sens_combined, res$sens)
}

level. <- c('Pre-Pandemic', 'Post-Pandemic')
ccf_combined$Period <- factor(ccf_combined$Period, levels = level.)
granger_combined$Period <- factor(granger_combined$Period, levels = level.)
ngranger_combined$Period <- factor(ngranger_combined$Period, levels = level.)
ccm_summary$Period <- factor(ccm_summary$Period, levels = level.)
ccm_plot_combined$Period <- factor(ccm_plot_combined$Period, levels = level.)
sens_combined$Period <- factor(sens_combined$Period, levels = level.)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my.theme <- theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold', size = base.size * 0.8, family = base.family, color = base.col),
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    axis.text.y = element_text(size = base.size * 0.8, family = base.family, color = base.col, angle = 90, hjust = 0.5),
    axis.title.y.left = element_text(size = base.size * 0.8, family = base.family, color = base.col, face = 'bold', margin = margin(r = 3)),
    axis.title.x.bottom = element_text(size = base.size * 0.8, family = base.family, color = base.col, face = 'bold', margin = margin(t = 8))
  )

# A: Time Series
dat.long <- dat.ana %>%
  dplyr::select(date, Period, Search.Z, Flu.Z) %>%
  pivot_longer(cols = c('Search.Z', 'Flu.Z'), names_to = 'Type', values_to = 'Value') %>%
  mutate(Type = factor(Type, labels = c('Influenza rate', 'Search index')))

p1 <- ggplot(dat.long, aes(x = date, y = Value, color = Type)) +
  geom_line(size = 0.6, alpha = 0.9) +
  facet_grid(~Period, scales = 'free_x', space = 'free_x') +
  labs(y = 'Intensity (Z-Score)', x = NULL) +
  scale_color_manual(values = c(col_flu, col_search)) +
  scale_y_continuous(expand = c(0, 0)) +
  my.theme +
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, -5, -5, -5),
    legend.position = c(0.10, 0.90)
  ) +
  facetted_pos_scales(
    x = list(
      Period == "Pre-Pandemic" ~ scale_x_date(
        date_breaks = "1 year",
        date_labels = "%b\n%Y",
        expand = c(0, 0)
      ),
      Period == "Post-Pandemic" ~ scale_x_date(
        date_breaks = "6 months",
        date_labels = "%b\n%Y",
        expand = c(0, 0)
      )
    )
  )

# B: CCF
max_ccf <- ccf_combined %>% group_by(Period) %>% slice(which.max(ACF))
p2 <- ggplot(ccf_combined, aes(x = Lag, y = ACF)) +
  geom_bar(stat = 'identity', fill = 'gray80', width = 0.6) +
  geom_bar(data = max_ccf, stat = 'identity', fill = col_search, width = 0.6) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 'dashed', color = 'black', size = 0.2) +
  scale_y_continuous(limits = c(-0.5, 1.0), expand = c(0, 0)) +
  scale_x_continuous(expand = c(ifelse(control, 0.3, 0.05), ifelse(control, 0.05, 0))) +
  facet_wrap(~Period) +
  labs(y = 'Correlation (r)', x = 'Lag (Weeks)') +
  geom_text(
    data = max_ccf,
    aes(label = paste0('Lag:', Lag, '\nr=', sprintf('%.3f', ACF))),
    vjust = ifelse(max_ccf$ACF > 0, -0.2, 1.2),
    size = base.size / 3.88,
    fontface = 'bold', 
    family = base.family
  ) +
  my.theme +
  theme(
    axis.title.x.bottom = element_text(size = base.size * 0.8, family = base.family, color = base.col, face = 'bold', margin = margin(t = -0.5, unit = 'cm'))
  )

# C: Granger causality test (Linear)
granger_combined <- granger_combined %>%
  mutate(Sig = case_when(
    P_val < 0.001 ~ '***',
    P_val < 0.01 ~ '**',
    P_val < 0.05 ~ '*',
    TRUE ~ 'NS'
  ))

p3 <- ggplot(granger_combined, aes(x = Direction, y = F_stat, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.6, alpha = 0.8) +
  geom_text(
    aes(
      label = Sig,
      y = ifelse(F_stat > ifelse(control, 2, 5), F_stat / 2, F_stat),
      vjust = ifelse(F_stat > 5, 0.5, -0.5)
    ),
    color = base.col,
    family = base.family,
    size = base.size / 3.88
  ) +
  facet_wrap(~Period) +
  scale_fill_manual(values = c(col_cause, col_reverse)) +
  labs(y = 'Granger F-Score', x = NULL) +
  my.theme +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col, lineheight = 0.60)
  )

# C2: Non-linear Granger causality test
ngranger_combined <- ngranger_combined %>%
  mutate(Sig = case_when(
    P_val < 0.001 ~ '***',
    P_val < 0.01 ~ '**',
    P_val < 0.05 ~ '*',
    TRUE ~ 'NS'
  ))

p3b <- ggplot(ngranger_combined, aes(x = Direction, y = Stat, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.6, alpha = 0.8) +
  geom_text(
    aes(
      label = Sig,
      y = ifelse(Stat > ifelse(this == 'Macau', 4, 4), Stat / 2, Stat),
      vjust = ifelse(Stat > 5, 0.5, -0.5)
    ),
    color = base.col,
    family = base.family,
    size = base.size / 3.88
  ) +
  facet_wrap(~Period) +
  scale_fill_manual(values = c(col_cause, col_reverse)) +
  labs(y = 'Non-linear Granger\n(Deviance)', x = NULL) +
  my.theme +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col, lineheight = 0.60)
  )

# D: CCM + Surrogates limit
# ccm_plot_combined contains Observed and Null 95% curves
if (nrow(ccm_plot_combined) > 0) {
  p4 <- ggplot(ccm_plot_combined, aes(x = LibSize, y = Rho)) +
    
    # 1. Plot the Null 95% lines (Dashed)
    # We map color to Direction so they inherit the same colors as the solid lines
    # We set alpha to make them slightly lighter/fainter than the main lines if desired, or keep them solid color
    geom_line(
      data = filter(ccm_plot_combined, Class == 'Null 95%'),
      aes(color = Direction, linetype = Class),
      linewidth = 0.8,
      alpha = 0.7
    ) +
    
    # 2. Plot the Observed lines (Solid)
    geom_line(
      data = filter(ccm_plot_combined, Class == 'Observed'),
      aes(color = Direction),
      linewidth = 1.2
    ) +
    
    facet_grid(~Period, scales = 'free_x', space = 'free_x') +
    
    # Define colors for both directions
    scale_color_manual(values = c('Search to Influenza' = col_cause, 'Influenza to Search' = col_reverse)) +
    
    # Define linetype manually
    scale_linetype_manual(values = c('Null 95%' = 'dashed'), name = '') +
    
    labs(x = 'Library Size (L)', y = 'Prediction Skill (rho)') +
    scale_y_continuous(expand = c(0, 0)) +
    my.theme +
    theme(
      legend.background = element_rect(fill = NA),
      plot.margin = margin(t = -.1, unit = 'cm'),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
  
  if (this == 'Macau') {
    if  (control) {
      p4 <- p4 + theme(legend.position = c(0.55, 0.35))
    } else {
      p4 <- p4 + theme(legend.position = c(0.55, 0.75))
    }
  } else {
    p4 <- p4 + theme(legend.position = c(0.25, 0.38))
  }
} else {
  p4 <- ggplot() + theme_void() + labs(title = 'CCM failed (insufficient data).')
}

# E: CCM E sensitivity
if (nrow(sens_combined) > 0) {
  p5 <- ggplot(sens_combined, aes(x = factor(E), y = Rho, fill = Direction)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
    facet_wrap(~Period) +
    scale_fill_manual(values = c('Search to Influenza' = col_cause, 'Influenza to Search' = col_reverse)) +
    labs(y = 'Cross Map Skill (Rho)', x = 'Embedding Dimension (E)') +
    my.theme +
    theme(
      legend.position = 'none',
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
} else {
  p5 <- ggplot() + theme_void()
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
layout_design <- '
AAAAAA
BBCCDD
EEFFFF
'

gg <- p1 / (p3 + p3b + p2) / (p5 + p4) /
  plot_layout(heights = c(1, 0.8, 1))

print(gg)

# Save to file
width = 14; height = 14
if (this == 'Macau' && control) {
  ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}_control.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')
} else {
  ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')
}

end.time <- Sys.time()
run.time <- end.time - start.time
print(run.time)
