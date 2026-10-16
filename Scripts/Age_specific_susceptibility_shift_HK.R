# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(ggtext)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 7
base.family <- 'sans'
base.col <- '#000000'
flu <- 'FLUAB'
s.week <- 20
age.pal <- c('#404040', '#E50611', '#F781BF', '#E67E24', '#19A955', '#3D99DA', '#831763')

# Vectorized formatter, works on scalars and vectors
format_p_exact <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < .Machine$double.eps) {
      return(sprintf('&lt; 2.220 × 10<sup>%d</sup>', floor(log10(.Machine$double.eps))))
    }
    if (x >= 0.001) {
      return(formatC(x, format = 'f', digits = 3))
    }
    exponent <- floor(log10(x))
    sprintf('%.3f × 10<sup>%d</sup>', x / 10^exponent, exponent)
  }, character(1))
}

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
    # z-test P value, consistent with the delta-method CI
    p_exact = 2 * pnorm(abs(log(RR) / se_log_rr), lower.tail = FALSE),
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
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = significance), height = 0.3, linewidth = 0.5, position = position_nudge(y = nudge_amount)) +
  geom_point(aes(color = significance), size = 2.5, position = position_nudge(y = nudge_amount)) +
  geom_richtext(
    aes(label = sprintf('%.3f [%.3f - %.3f], *P* = %s', RR, lower_ci, upper_ci, format_p_exact(p_exact))), 
    position = position_nudge(y = nudge_amount),
    vjust = -1.25,
    size = base.size / ggplot2::.pt,
    color = '#000000',
    fill = NA,
    label.color = NA,
    label.padding = unit(0, 'lines')
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
  guides(
    color = guide_legend(
      override.aes = list(
        size = 1.5,
        linewidth = 0.5,
        linetype = 1 
      )
    )
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
    legend.key.spacing.y = unit(-0.20, units = 'cm'),
    legend.text = element_text(size = base.size, family = base.family, color = base.col),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_blank()
  )

# ----------------------------------------------
# Combine Plots
# ----------------------------------------------
gg <- p_f

print(dat.rr %>% dplyr::select(Age, RR, lower_ci, upper_ci, significance))

# Save to file
width = 180; height = 60
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{flu}.pdf'), width = width, height = height, units = 'mm', device = cairo_pdf, bg = '#FFFFFF')


# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
fig.label <- 'Fig.S13'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  writeLines(c(
    glue('Source data for {fig.label}{panel}: {title}')
  ), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
  close(con)
}

write_source_txt(
  dat.rr,
  panel = '',
  title = 'Rate ratios (RR, 95% CI) of post- vs pre-pandemic influenza susceptibility by age group, Hong Kong surveillance data'
)