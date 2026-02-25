# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ISOweek)))
suppressMessages(suppressWarnings(library(stringr)))
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
this <- 'Macau'  # Options: 'HK' or 'Macau'

# Plotting Constants
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Load data
dat.fac <- readxl::read_excel(glue('{dat.dir}/{this}/influenza.xlsx'))

# Prepare plot data
sdate <- '2020-01-01'
edate <- '2023-01-01'

dat.pl <- dat.fac %>% dplyr::select(date, FLUA, FLUB) %>%
  arrange(date) %>%
  mutate(TimeIndex = row_number()) %>%
  mutate(
    Period = case_when(
      date < sdate ~ 'Pre-COVID',
      date > edate ~ 'Post-COVID',
      TRUE ~ 'COVID'
    )
  )

dat.long <- dat.pl %>%
  pivot_longer(cols = c('FLUA', 'FLUB'), names_to = 'Subtype', values_to = 'Cases')

df.prop <- dat.long %>%
  group_by(date) %>%
  mutate(Total = sum(Cases)) %>%
  ungroup() %>%
  mutate(Proportion = ifelse(Total == 0, 0, Cases / Total))

# Colors
cols <- c('FLUA' = '#E50914', 'FLUB' = '#00A087')

# ++++++++++++++++++++++++++++++++++++++++++
# Plot A: Absolute Epidemic Intensity (Showing Rebound Height)
# ++++++++++++++++++++++++++++++++++++++++++
p1 <- ggplot(dat.long, aes(x = as.Date(date), y = Cases, fill = Subtype)) +
  geom_area(alpha = 0.6, position = 'identity') + 
  geom_smooth(aes(color = Subtype), se = FALSE, span = 0.1, size = 0.8) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  scale_x_date(
    date_breaks = '2 years',
    date_labels = '%Y',
    limits = as.Date(c(ifelse(this == 'Macau', '2010-01-01', '2014-01-01'), '2026-01-01')), 
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  annotate('rect', xmin = as.Date(sdate), xmax = as.Date(edate), ymin = -Inf, ymax = Inf, alpha = 0.1, fill = 'gray30') +
  annotate('text', x = as.Date('2021-07-31'), y = 0.25, label = 'NPIs / COVID-19 Pandemic', color = base.col, size = base.size / 3.88, family = base.family) +
  labs(subtitle = NULL, y = 'Weekly Cases', x = NULL) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'none',
    text = element_text(size = base.size, family = base.family, color = base.col),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.margin = margin(r = 20, b = 20)
  )
    
# ++++++++++++++++++++++++++++++++++++++++++
# Plot B: Relative Competitive Landscape (Showing Strain Replacement)
# ++++++++++++++++++++++++++++++++++++++++++
p2 <- ggplot(df.prop, aes(x = as.Date(date), y = Proportion, fill = Subtype)) +
  geom_area(alpha = 0.8, color = 'white', size = 0.1) +
  scale_fill_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    name = 'Influenza type',
    labels = c('FLUA' = 'Influenza A', 'FLUB' = 'Influenza B'),
    guide = guide_legend(title.position = 'top', title.hjust = 0.5) 
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_date(
    date_breaks = '2 years',
    date_labels = '%Y',
    limits = as.Date(c(ifelse(this == 'Macau', '2010-01-01', '2014-01-01'), '2026-01-01')), 
    expand = c(0, 0)
  ) +
  annotate('rect', xmin = as.Date(sdate), xmax = as.Date(edate), ymin = -Inf, ymax = Inf, alpha = 0.1, fill = 'gray30') +
  labs(subtitle = 'Strain competition: relative share (influenza A vs B)', y = 'Proportion of Cases', x = 'Year') +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(size = base.size, family = base.family, color = base.col),
    legend.position = 'inside',
    legend.position.inside = c(0.80, 0.72),
    legend.background = element_rect(fill = '#FFFFFF66'),
    legend.direction = 'horizontal',
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold'),
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.margin = margin(r = 20)
  )

# Combination
gg <- p1 / p2 + plot_layout(heights = c(1, 1))

# Show plot
print(gg)

# Save to file
width = 12; height = 5
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')