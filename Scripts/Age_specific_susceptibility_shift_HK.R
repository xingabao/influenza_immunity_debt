# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
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
flu <- 'FLUAB'
s.week <- 20
age.pal <- c('#404040', '#E50611', '#F781BF', '#E67E24', '#19A955', '#3D99DA', '#831763')

# Load data
dat.raw <- readxl::read_excel(glue('{dat.dir}/HK/FLU-CL-AQ-age.xlsx'))
dat.fac <- dat.raw %>% dplyr::select(date, Age, all_of(flu))

# Prepare data
dat.sur <- dat.fac %>%
  filter(!is.na(.data[[flu]])) %>%
  mutate(
    date = as.POSIXct(date),
    Age = factor(Age, levels = c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 49)', '[50, 64)', '65+')),
    year = year(date),
    week = isoweek(date),
    period = case_when(
      year <= 2019 ~ 'Pre_COVID',
      year >= 2023 ~ 'Post_COVID',
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
# Calculate the Risk Ratio (RR) and plot a Forest Plot
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
    RR = mean_val_Post_COVID / mean_val_Pre_COVID,
    se_log_rr = sqrt((sd_val_Post_COVID^2 / (n_Post_COVID * mean_val_Post_COVID^2)) + (sd_val_Pre_COVID^2 / (n_Pre_COVID * mean_val_Pre_COVID^2)) ),
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
    legend.position = c(0.70, ifelse(flu == 'FLUA', yes = 0.25, no = 0.30)), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = 'transparent', color = NA),
    legend.key.height = unit(0.7, 'cm'),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_blank()
  )

# ----------------------------------------------
# Combine Plots
# ----------------------------------------------
gg <- p_f

print(dat.rr %>% dplyr::select(Age, RR, lower_ci, upper_ci, significance))

# Save to file
width = 12; height = 5
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{flu}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')