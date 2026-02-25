# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Define a function
is_integer_string_col <- function(x) {
  num_x <- suppressWarnings(as.numeric(x))
  if (sum(is.na(num_x)) > sum(is.na(x))) {
    return(FALSE)
  }
  
  if (all(is.na(num_x))) {
    return(FALSE)
  }
  
  all(num_x %% 1 == 0, na.rm = TRUE)
}

find_cols_with_keyword <- function(df, keyword = '人次') {
  
  check_col <- function(x) {
    x_char <- as.character(x)
    any(str_detect(x_char, keyword), na.rm = TRUE)
  }
  
  has_keyword_vec <- vapply(df, check_col, FUN.VALUE = logical(1))
  
  col_indices <- which(has_keyword_vec)
  
  return(col_indices)
}

# Load data
if (!exists('dat.tou')) {
  xlsxs <- list.files(glue('{dat.dir}/Macau/Tourist'))
  cola <- '總數'
  colb <- c('亞洲', '美洲', '歐洲', '大洋洲', '非洲及其他')
  colb <- c('總數\nTOTAL GERAL', '總數\nTotal')
  dat.tou <- data.frame()
  for (xlsx in xlsxs) {
    dftmp <- suppressMessages(suppressWarnings(readxl::read_excel(glue('{dat.dir}/Macau/Tourist/{xlsx}'), skip = 0, col_names = FALSE, sheet = '1'))) %>% as.data.frame()
    indices <- find_cols_with_keyword(dftmp, '人次')
    if (length(indices) != 0) {
      dftmp <- dftmp[, c('...1', colnames(dftmp)[indices])]
    }
    dftmp <- dftmp %>% filter(`...1` %in% c(cola, colb))
    lindices <- vapply(dftmp, is_integer_string_col, FUN.VALUE = logical(1))
    indices <- which(lindices)
    if (FALSE) print(glue('{xlsx}\t{dftmp[1, 1]}\t{dftmp[1, indices[1]]}, {dftmp[1, indices[2]]}'))
    tmp <- data.frame(year = str_extract(xlsx, '\\d{4}'), month = str_extract(xlsx, '(?<=M)\\d+'), num = dftmp[1, indices[2]])
    dat.tou <- rbind(dat.tou, tmp)
  }
  
  dat.tou <- dat.tou %>% as_tibble()
  dat.fac <- readxl::read_excel(glue('{dat.dir}/Macau/FLU-AQ-CLIMATE.xlsx')) %>% dplyr::select(Year, Week, Flu = FLUAB) %>% filter(Year > 2009)
}

# Arrange data
dat.tou.monthly <- dat.tou %>%
  as_tibble() %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    visitors = as.numeric(num),
    date = make_date(year, month, 15)
  ) %>%
  dplyr::select(date, year, month, visitors) %>%
  arrange(date) %>%
  filter(!is.na(visitors))

dat.flu.weekly <- dat.fac %>%
  mutate(
    Year = as.integer(Year),
    Week = as.integer(Week),
    Flu = as.numeric(Flu)
  ) %>%
  filter(!is.na(Year) & !is.na(Week) & !is.na(Flu)) %>%
  mutate(
    date = suppressWarnings(as.Date(paste(Year, Week, 1, sep = '-'), format = '%Y-%U-%u'))
  ) %>%
  filter(!is.na(date)) %>%
  arrange(date)

#
dat.agg <- dat.flu.weekly %>%
  mutate(month_date = floor_date(date, 'month')) %>%
  group_by(month_date) %>%
  summarise(flu_rate_monthly = mean(Flu, na.rm = TRUE)) %>%
  mutate(year = year(month_date), month = month(month_date))

dat.merged.stat <- inner_join(
  dat.tou.monthly %>% mutate(month_date = floor_date(date, 'month')), 
  dat.agg, 
  by = 'month_date'
) %>%
  mutate(month_fac = factor(month.x))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analysis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cor_res <- cor.test(dat.merged.stat$visitors, dat.merged.stat$flu_rate_monthly)
r_val <- round(cor_res$estimate, 3)
p_val_cor <- format.pval(cor_res$p.value, eps = 0.001, digits = 3)

dat.recent <- dat.merged.stat %>% filter(year.x >= 2023)
has_recent_data <- nrow(dat.recent) > 10

if (has_recent_data) {
  
  # Model: Influenza ~ Tourists + Month Factor
  fit_model <- lm(flu_rate_monthly ~ visitors + month_fac, data = dat.recent)
  sum_model <- summary(fit_model)
  
  coef_vis <- sum_model$coefficients['visitors', 'Estimate']
  pval_vis <- sum_model$coefficients['visitors', 'Pr(>|t|)']
  r2_adj   <- sum_model$adj.r.squared
  
  pval_vis_txt <- format.pval(pval_vis, eps = 0.001, digits = 3)
  if (substr(pval_vis_txt, 1, 1) == '<') {
    pval_vis_str <- paste0('P < ', substr(pval_vis_txt, 2, 100))
  } else {
    pval_vis_str <- paste0('P = ', pval_vis_txt)
  }
  
  stats_label <- paste0(
    'Statistical analysis:\n',
    '• Correlation (Pearson):', '\n    r = ', r_val, ' (P ', ifelse(cor_res$p.value < 0.001, '< 0.001', paste0('= ', p_val_cor)), ')\n',
    '• Regression (Adj. for seasonality):\n',
    '    Visitor Coef. = ', formatC(coef_vis, format = 'e', digits = 3), ', ', pval_vis_str, ', Adj. R2 = ', sprintf('%.3f', r2_adj)
  )
} else {
  stats_label <- 'Insufficient data for regression analysis'
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Dual-Frequency Plot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
max_flu <- max(dat.flu.weekly$Flu, na.rm = TRUE)
max_tourist <- max(dat.tou.monthly$visitors, na.rm = TRUE)

coeff <- max_tourist / max_flu * 1.25

# Colors
col_tourist <- '#00A087'
col_flu <- '#E50914'

# Plot
p_dual <- ggplot() +
  geom_col(data = dat.tou.monthly, aes(x = date, y = visitors / coeff, fill = 'Monthly tourists'), width = 25, alpha = 0.5) +
  geom_line(data = dat.flu.weekly, aes(x = date, y = Flu, color = 'Weekly influenza rate'), linewidth = 0.5, alpha = 0.9) +
  annotate(
    'label', x = as.Date('2010-11-10'), y = max_flu * 1.15, label = stats_label, 
    hjust = 0, vjust = 1, size = base.size / 3.88, family = base.family,
    fill = '#FFFFFF', alpha = 0.8, lineheight = 1.5
  ) +
  scale_x_date(NULL, date_breaks = '1 year', date_labels = '%b\n%Y', expand = c(0, 0), limits = as.Date(c('2010-01-01', '2026-01-01'))) +
  scale_y_continuous(
    name = 'Weekly Influenza Rate',
    expand = c(0, 0),
    limits = c(0, max_flu * 1.20), 
    labels = label_percent(accuracy = 0.1), 
    sec.axis = sec_axis(~ . * coeff, name = 'Monthly Visitor Arrivals', labels = label_number(accuracy = 0.1, scale_cut = cut_short_scale()))
  ) +
  scale_fill_manual(name = NULL, values = c('Monthly tourists' = col_tourist)) +
  scale_color_manual(name = NULL, values = c('Weekly influenza rate' = col_flu)) +
  labs(subtitle = 'Overlay of high-frequency surveillance (Weekly) and demographic flow (Monthly)') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family,color = col_flu, face = 'bold', margin = margin(r = 10)),
    axis.text.y.left = element_text(size = base.size * 0.8, family = base.family, color = col_flu),
    axis.line.y.left = element_line(color = col_flu),
    axis.title.y.right = element_text(size = base.size * 0.9, family = base.family,color = col_tourist, face = 'bold', angle = 90, margin = margin(l = 10)),
    axis.text.y.right = element_text(size = base.size * 0.8, family = base.family, color = col_tourist),
    axis.line.y.right = element_line(color = col_tourist),
    legend.position = c(0.65, 0.85),
    legend.direction = 'horizontal',
    legend.spacing.y = unit(-0.4, 'cm'),
    legend.background = element_blank()
  ); p_dual

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Seasonal Pattern
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dat.season <- dat.merged.stat %>%
  group_by(month.x) %>%
  summarise(
    mean_vis = mean(visitors, na.rm=TRUE),
    mean_flu = mean(flu_rate_monthly, na.rm=TRUE),
    se_vis = sd(visitors, na.rm=TRUE)/sqrt(n()),
    se_flu = sd(flu_rate_monthly, na.rm=TRUE)/sqrt(n())
  )

scale_season <- max(dat.season$mean_vis) / max(dat.season$mean_flu)

p_season <- ggplot(dat.season, aes(x = month.x)) +
  geom_col(aes(y = mean_vis / scale_season, fill = 'Avg. Tourists'), alpha = 0.6, width = 0.7) +
  geom_errorbar(aes(ymin = (mean_vis - se_vis) / scale_season, ymax = (mean_vis + se_vis) / scale_season), width = 0.1, color = col_tourist, linewidth = 0.5, position = position_nudge(x = -0.05)) +
  geom_line(aes(y = mean_flu, group = 1, color = 'Avg. Flu Rate'), linewidth = 0.5, position = position_nudge(x = 0.15)) +
  geom_point(aes(y = mean_flu, color = 'Avg. Flu Rate'), size = 1.5, position = position_nudge(x = 0.15)) +
  geom_errorbar(aes(ymin = mean_flu - se_flu, ymax = mean_flu + se_flu), width = 0.1, color = col_flu, linewidth = 0.5, position = position_nudge(x = 0.15)) +
  scale_y_continuous(
    name = 'Average Influenza Rate',
    expand = c(0, 0),
    limits = c(0, 0.11),
    labels = label_percent(accuracy = 0.1), 
    sec.axis = sec_axis(~ . * scale_season, name = 'Average Visitors', labels = label_number(accuracy = 0.1, scale_cut = cut_short_scale())) 
  ) +
  scale_x_continuous(breaks = 1:12, labels = glue('{month.abb}.'), name = NULL) +
  scale_fill_manual(values = c('Avg. Tourists' = col_tourist)) +
  scale_color_manual(values = c('Avg. Flu Rate' = col_flu)) +
  labs(subtitle = 'Seasonal aggregation (monthly average)', x = NULL) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'none', 
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family,color = col_flu, face = 'bold', margin = margin(r = 10)),
    axis.text.y.left = element_text(size = base.size * 0.8, family = base.family, color = col_flu),
    axis.line.y.left = element_line(color = col_flu),
    axis.title.y.right = element_text(size = base.size * 0.9, family = base.family,color = col_tourist, face = 'bold', angle = 90, margin = margin(l = 10)),
    axis.text.y.right = element_text(size = base.size * 0.8, family = base.family, color = col_tourist),
    axis.line.y.right = element_line(color = col_tourist),
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col)
  ); p_season

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gg <- p_dual / p_season + plot_layout(heights = c(2, 1))

# Show plot
print(gg)

# Save to file
width = 12; height = 9
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')