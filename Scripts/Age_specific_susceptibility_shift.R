# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 16
base.family <- 'serif'
base.col <- '#000000'
flu <- 'FLUA' # Options: 'FLUA' or 'FLUB'
if (flu == 'FLUA') {
  s.week <- 18
} else {
  s.week <- 12
}
age.pal <- c('#404040', '#E50611', '#F781BF', '#E67E24', '#19A955', '#3D99DA', '#831763')

# Load data
dat.raw <- readxl::read_excel(glue('{dat.dir}/Macau/FLU-CL-AQ-age.xlsx'))
dat.fac <- dat.raw %>% dplyr::select(date, Age, all_of(flu))

# Prepare data
dat.sur <- dat.fac %>%
  filter(!is.na(.data[[flu]])) %>%
  mutate(
    date = as.POSIXct(date),
    Age = factor(Age, levels = c('[0, 5)', '[5, 10)', '[10, 20)', '[20, 35)', '[35, 50)', '[50, 65)', '65+')),
    year = year(date),
    week = isoweek(date),
    period = case_when(
      year <= 2019 ~ 'Pre_Pandemic',
      year >= 2023 ~ 'Post_Pandemic',
      TRUE ~ 'Excluded'
    )
  ) %>%
  filter(period != 'Excluded')

# Influenza season alignment - Dynamic calculation
dat.aln <- dat.sur %>%
  mutate(
    season_week = ifelse(week >= s.week, week - (s.week - 1), week + (52 - (s.week - 1)))
  ) %>%
  filter(season_week > 0 & season_week <= 53)

# ----------------------------------------------
# 1. Create "survival curves" (Figures A & B)
# To illustrate "seasonal epidemic time shifts"
# ----------------------------------------------
# Calculate the weekly average positivity rate and accumulate it.
dat.survival <- dat.aln %>%
  group_by(period, Age, season_week) %>%
  summarise(mean_rate = mean(.data[[flu]], na.rm = TRUE), .groups = 'drop') %>%
  arrange(period, Age, season_week) %>%
  group_by(period, Age) %>%
  mutate(
    cum_rate = cumsum(mean_rate),
    cdf = cum_rate / max(cum_rate),
    inv_cdf = 1 - (cum_rate / max(cum_rate) * 0.95)
  ) %>%
  ungroup()

# K-S Test
get_p_value_labels <- function(df, target_period, ref_age = '[0, 5)') {
  
  ref_data <- df %>% 
    filter(period == target_period, Age == ref_age) %>% 
    pull(cdf)
  
  p_values <- df %>%
    filter(period == target_period) %>%
    group_by(Age) %>%
    summarise(
      test_data = list(cdf),
      .groups = 'drop'
    ) %>%
    rowwise() %>%
    mutate(
      ks_res = list(ks.test(unlist(test_data), ref_data)),
      p_val = ks_res$p.value,
      p_label = case_when(
        Age == ref_age ~ glue('{Age} (Ref.)'),
        p_val < 0.001  ~ glue("{Age} (*P* < 0.001)"),
        p_val < 0.05   ~ glue("{Age} (*P* = {format(p_val, digits=2)})"),
        p_val < 0.01   ~ glue("{Age} (*P* = {format(p_val, digits=2)})"),
        TRUE           ~ glue('{Age} ~ (NS)')
      )
    ) %>%
    mutate(p_label = fct_inorder(p_label)) %>%
    dplyr::select(Age, p_label)
  
  return(p_values)
}

# Calculate Pre and Post separately.
labels_pre <- get_p_value_labels(dat.survival, 'Pre_Pandemic', ref_age = '[0, 5)')
labels_post <- get_p_value_labels(dat.survival, 'Post_Pandemic', ref_age = '[0, 5)')

# Plot function
plot.survival.with.p <- function(data_sub, title_str, label_df) {
  
  plot_data <- data_sub %>%
    left_join(label_df, by = 'Age')
  
  # K-S Test vs Reference Group [0, 5)
  ggplot(plot_data, aes(x = season_week, y = inv_cdf, color = p_label, group = Age)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    scale_color_manual(values = age.pal) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      subtitle = title_str,
      x = NULL,
      y = 'Proportion of Total Cases Remaining',
      color = 'Age Group'
    ) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text = element_text(size = base.size, family = base.family, color = base.col),
      axis.title.y = element_text(size = base.size, family = base.family, color = base.col, face = 'bold', margin = margin(l = 0.5, unit = 'cm')),
      legend.position = c(0.05, 0.60), 
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = 'transparent', color = NA),
      legend.key.height = unit(0.6, 'cm'),
      legend.text = element_markdown(size = base.size * 0.9, family = base.family, color = base.col),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    guides(color = guide_legend(ncol = 1))
}

p_a <- plot.survival.with.p(dat.survival %>% filter(period == 'Pre_Pandemic'), 'Pre-Pandemic (2010-2019)', labels_pre)
p_b <- plot.survival.with.p(dat.survival %>% filter(period == 'Post_Pandemic'), 'Post-Pandemic (2023-2025)', labels_post) + labs(y = NULL)

# ----------------------------------------------
# 2.Calculate cumulative incidence curves (for Figures C & D)
# a.Determine the average incidence rate for each period (Pre/Post), each age group, and each week (1-40) during the flu season.
# b. Calculate the cumulative sum -> This represents the cumulative incidence curve.
# ----------------------------------------------
dat.curve <- dat.aln %>%
  group_by(period, Age, season_week) %>%
  summarise(mean_rate = mean(.data[[flu]], na.rm = TRUE), .groups = 'drop') %>%
  arrange(period, Age, season_week) %>%
  group_by(period, Age) %>%
  mutate(cumulative_rate = cumsum(mean_rate)) 

get_magnitude_p_values <- function(df, target_period, ref_age = '[0, 5)') {
  ref_data <- df %>% 
    filter(period == target_period, Age == ref_age) %>% 
    arrange(season_week) %>%
    pull(mean_rate)
  
  p_values <- df %>%
    filter(period == target_period) %>%
    group_by(Age) %>%
    summarise(
      test_data = list(mean_rate),
      .groups = 'drop'
    ) %>%
    rowwise() %>%
    mutate(
      wilcox_res = list(tryCatch(
        wilcox.test(unlist(test_data), ref_data, paired = TRUE),
        error = function(e) list(p.value = 1)
      )),
      p_val = wilcox_res$p.value,
      p_label = case_when(
        Age == ref_age ~ glue('{Age} (Ref.)'),
        p_val < 0.001  ~ glue("{Age} (*P* < 0.001)"),
        p_val < 0.05   ~ glue("{Age} (*P* = {format(p_val, digits=2)})"),
        p_val < 0.01   ~ glue("{Age} (*P* = {format(p_val, digits=2)})"),
        TRUE           ~ glue('{Age} ~ (NS)')
      )
    ) %>%
    mutate(p_label = fct_inorder(p_label)) %>%
    dplyr::select(Age, p_label)
  
  return(p_values)
}

labels_mag_pre <- suppressWarnings(get_magnitude_p_values(dat.curve, 'Pre_Pandemic', ref_age = '[0, 5)'))
labels_mag_post <- suppressWarnings(get_magnitude_p_values(dat.curve, 'Post_Pandemic', ref_age = '[0, 5)'))

my.pal <- brewer.pal(7, 'Dark2') 

plot.curve <- function(data_sub, title_str, label_df) {
  plot_data <- data_sub %>%
    left_join(label_df, by = 'Age')
  
  ggplot(plot_data, aes(x = season_week, y = cumulative_rate, color = p_label, group = Age)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    scale_color_manual(values = age.pal) +
    labs(
      subtitle = title_str,
      x = glue('Weeks into Flu Season'),
      y = 'Cumulative % of Seasonal Total',
      color = 'Age Group' 
    ) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text = element_text(size = base.size, family = base.family, color = base.col),
      axis.title.x = element_text(size = base.size, family = base.family, face = 'bold', color = base.col, margin = margin(t = 0.2, unit = 'cm')),
      axis.title.y = element_text(size = base.size, family = base.family, face = 'bold', color = base.col, margin = margin(l = 0.5, unit = 'cm')),
      legend.position = c(0.05, 0.95), 
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = 'transparent', color = NA),
      legend.key.height = unit(0.6, 'cm'),
      legend.text = element_markdown(size = base.size * 0.9, family = base.family, color = base.col),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    guides(color = guide_legend(ncol = 1))
}

p_pre <- plot.curve(dat.curve %>% filter(period == 'Pre_Pandemic'), 'Pre-Pandemic: Baseline Accumulation', labels_mag_pre)
p_post <- plot.curve(dat.curve %>% filter(period == 'Post_Pandemic'), 'Post-Pandemic: Accelerated Accumulation', labels_mag_post) + labs(y = NULL)

# ----------------------------------------------
# 3.Calculate the Risk Ratio (RR) and plot a Forest Plot (Figure F).
# ----------------------------------------------
dat.rr <- dat.aln %>%
  group_by(period, Age) %>%
  summarise(
    mean_val = mean(.data[[flu]], na.rm = TRUE),
    sd_val = sd(.data[[flu]], na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = period, values_from = c(mean_val, sd_val, n)) %>%
  mutate(
    RR = mean_val_Post_Pandemic / mean_val_Pre_Pandemic,
    se_log_rr = sqrt((sd_val_Post_Pandemic^2 / (n_Post_Pandemic * mean_val_Post_Pandemic^2)) + (sd_val_Pre_Pandemic^2 / (n_Pre_Pandemic * mean_val_Pre_Pandemic^2)) ),
    lower_ci = exp(log(RR) - 1.96 * se_log_rr),
    upper_ci = exp(log(RR) + 1.96 * se_log_rr),
    significance = case_when(
      lower_ci > 1 ~ 'Significant Increase (Debt)',
      upper_ci < 1 ~ 'Significant Decrease',
      TRUE ~ 'No Significant Change'
    )
  )

# Ratio > 1 indicates higher susceptibility intensity
nudge_amount <- -0.25 
p_f <- ggplot(dat.rr, aes(x = RR, y = fct_rev(Age))) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50') +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = significance), height = 0.3, linewidth = 0.8, position = position_nudge(y = nudge_amount)) +
  geom_point(aes(color = significance), size = 4, position = position_nudge(y = nudge_amount)) +
  geom_text(
    aes(label = sprintf('%.3f [%.3f - %.3f]', RR, lower_ci, upper_ci)), 
    position = position_nudge(y = nudge_amount),
    vjust = -1.25,
    size = base.size / 3.88,
    color = '#000000'
  ) +
  scale_x_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_color_manual(
    values = c(
      'Significant Increase (Debt)' = '#F34149', 
      'No Significant Change' = '#95A5A6', 
      'Significant Decrease' = '#2ECC71'
    )
  ) +
  labs(
    subtitle = 'Immunity debt evidence: Rate Ratio (Post-Pandemic / Pre-Pandemic)',
    x = 'Rate Ratio [95% CI]',
    y = ''
  ) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = base.size, family = base.family, face = 'bold', color = base.col, margin = margin(t = 0.2, unit = 'cm')),
    axis.text = element_text(size = base.size, family = base.family, color = base.col),
    legend.position = c(0.70, ifelse(flu == 'FLUA', yes = 0.30, no = 0.30)), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = 'transparent', color = NA),
    legend.key.height = unit(0.7, 'cm'),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_blank()
  )

# ----------------------------------------------
# 4.Combine Plots
# ----------------------------------------------
gg <- (p_a + p_b) / (p_pre + p_post) / p_f +
  plot_layout(heights = c(1, 1, 0.8))

print(dat.rr %>% dplyr::select(Age, RR, lower_ci, upper_ci, significance))

# Save to file
width = 12; height = 13
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{flu}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')