# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(broom)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(lubridate)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'Macau'
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'
flu <- 'FLUA'  # Options: 'FLUA' or 'FLUB'

# Load data
dat.raw <- readxl::read_excel(glue('{dat.dir}/Macau/FLU-CL-AQ-day.oe.xlsx'))

# Prepare plot data
sdate <- '2020-01-01'
edate <- '2023-01-01'
envs <- c('ave.temp')


if (flu == 'FLUA') { noTest = 'noFLUA'; neTest = 'neFLUA'} else { noTest = 'noFLUB'; neTest = 'neFLUB' }

dat.fac <- dat.raw %>% dplyr::select(date, holiday, nSample, all_of(noTest), all_of(neTest)) %>%
  arrange(date) %>%
  mutate(TimeIndex = row_number()) %>%
  mutate(
    Period = case_when(
      date < sdate ~ 'Pre-Pandemic',
      date > edate ~ 'Post-Pandemic',
      TRUE ~ 'Pandemic'
    )
  ) %>%
  mutate(
    ER.ratio = .data[[neTest]] / (.data[[neTest]] + .data[[noTest]]),
  ) %>% filter(Period %in% c('Pre-Pandemic', 'Post-Pandemic'))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Verify whether the emergency burden has significantly increased in the post-pandemic era
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# We only compare Pre-Pandemic and Post-Pandemic periods, excluding the interference from the Pandemic period.
dat.ana <- dat.fac %>% 
  filter(Period != 'Pandemic') %>%
  mutate(Period = factor(Period, levels = c('Pre-Pandemic', 'Post-Pandemic'))) %>%
  dplyr::select(all_of(neTest), all_of(noTest), ER.ratio, Period) %>%
  na.omit()

model.burden <- glm(cbind(get(neTest), get(noTest)) ~ Period, data = dat.ana, family = binomial(link = 'logit'))
model.result <- tidy(model.burden, exponentiate = TRUE, conf.int = TRUE) 
print(model.result)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Visualization
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Colors
color.opd <- '#00A087'
color.er <- '#E50914'

# Severity Ratio
# Structural Shift in Healthcare Burden
p2 <- ggplot(dat.ana, aes(x = Period, y = ER.ratio, fill = Period)) +
  geom_violin(alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.2, color = 'black', alpha = 0.8,) +
  scale_fill_manual(values = c('Pre-Pandemic' = color.opd, 'Post-Pandemic' = color.er)) +
  scale_y_continuous(limits = c(-0.1, 1.4), breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', '1.25'), expand = c(0, 0)) +
  labs(
    subtitle = 'Comparison of ED admission ratios',
    y = 'Proportion of cases in ED'
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  stat_summary(fun = mean, geom = 'point', shape = 23, size = 3, fill = 'white') +
  stat_compare_means(
    comparisons = list(c('Pre-Pandemic', 'Post-Pandemic')),
    method = 'wilcox.test',
    label = 'p.signif',
    vjust = 0.5
  ) + theme(
    legend.position = 'none',
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(r = 1)),
    axis.title.x.bottom = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5)),
    axis.text = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.margin = margin(t = 0.05, b = 0, unit = 'in')
  ); p2

# Forest Plot
dat.plot <- model.result %>% 
  filter(term != '(Intercept)') %>%
  mutate(
    label_text = sprintf("%.3f [%.3f - %.3f]", estimate, conf.low, conf.high)
  )
p3 <- ggplot(dat.plot %>% filter(term != '(Intercept)'), aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), color = color.er, width = 0.03, linewidth = 1.0) + 
  geom_point(color = color.er, size = 3, shape = 18) +
  geom_text(
    aes(label = label_text), 
    vjust = -1.5,
    size = base.size * 0.25,
    family = base.family,
    color = base.col
  ) +
  labs(
    subtitle = 'Odds Ratio (Post-Pandemic vs. Pre-Pandemic baseline)',
    x = 'Odds Ratio (95% CI)',
    y = ''
  ) +
  scale_x_continuous(limits = c(0, 3.5), breaks = c(0, 1, 2, 3), expand = c(0, 0)) +
  scale_y_discrete(labels = c('PeriodPost-Pandemic' = 'Post-Pandemic Effect')) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.margin = margin(t = 0.05, b = 0, unit = 'in'),
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    axis.text.y = element_text(size = base.size * 0.8, family = base.family, color = base.col, angle = 90, hjust = 0.5),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(r = 1)),
    axis.title.x.bottom = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5))
  ); p3

# Combination
gg <- p2 | p3

# Show plot
print(gg)

# Save to file
width = 12; height = 4
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{flu}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')