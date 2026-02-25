# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ISOweek)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(cowplot)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'
flu <- 'FLUA'  # Options: 'FLUA' or 'FLUB'

# Load data
if (file.exists(glue('{dat.dir}/Macau/INFLU-U2.csv.gz'))) {
  # The primary clinical dataset comprising individual-level electronic health records 
  # from Kiang Wu Hospital is not publicly available due to patient privacy regulations
  # and ethical restrictions regarding the protection of personal health information;
  # however, anonymized data supporting the findings of this study may be made available 
  # to qualified researchers upon reasonable request to the corresponding authors,
  # subject to approval by the institutional review board and the execution of a data sharing agreement.
  dat.raw <- data.table::fread(glue('{dat.dir}/Macau/INFLU-U2.csv.gz'))
  
  # Prepare plot data
  dat.pl <- dat.raw %>% dplyr::select(Year, Week, Local, FLUA, FLUB) %>% filter(Local != -1)
  saveRDS(dat.pl, file = glue('{dat.dir}/Macau/Tourist_vs_Local_dual_axis.rds'))
} else {
  dat.pl <- readRDS(glue('{dat.dir}/Macau/Tourist_vs_Local_dual_axis.rds'))
}

process_virus_counts <- function(dt, virus_name) {
  
  if (virus_name == 'FLUA') {
    dt <- dt %>% filter(FLUA == 1)
  } else {
    dt <- dt %>% filter(FLUB == 1)
  }
  
  counts <- dt[get(virus_name) > 0, .(Count = .N), by = .(Year, Week, Local)]
  
  counts[, Week_Date := ISOweek2date(sprintf('%04d-W%02d-1', Year, Week))]
  
  res <- dcast(counts, Week_Date ~ Local, value.var = 'Count', fun.aggregate = sum, fill = 0)
  
  setnames(res, 'Week_Date', 'Week')
  if ('0' %in% names(res)) setnames(res, '0', 'Tourist') else res[, Tourist := 0]
  if ('1' %in% names(res)) setnames(res, '1', 'Local')   else res[, Local := 0]
  
  setcolorder(res, c('Week', 'Tourist', 'Local'))
  setorder(res, Week)
  
  return(res)
}

df_ts_A <- process_virus_counts(dat.pl, 'FLUA')
df_ts_B <- process_virus_counts(dat.pl, 'FLUB')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A: Dual-Axis Time Series
# +++++++++++++++++++++++++++++++++++++++++++++++++++
if (flu == 'FLUA') {
  df_ts <- df_ts_A
} else {
  df_ts <- df_ts_B
}

scale_factor <- max(df_ts$Local) / max(df_ts$Tourist)

plot_A <- ggplot(df_ts, aes(x = Week)) +
  
  geom_line(aes(y = Local, color = 'Local Residents'), linewidth = 1.0, linetype = 'dashed') +
  geom_area(aes(y = Local), fill = '#00A087', alpha = 0.1) +

  geom_line(aes(y = Tourist * scale_factor, color = 'Tourists'), linewidth = 1.0) +
  
  scale_y_continuous(
    name = 'Local Cases (Weekly)',
    sec.axis = sec_axis(~ . / scale_factor, name = 'Tourist Cases (Weekly)'),
    expand = c(0, 0)
  ) +
  scale_x_date(
    date_breaks = '1 year',
    date_labels = '%Y',
    limits = as.Date(c('2016-01-01', '2026-01-01')), 
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c('Local Residents' = '#00A087', 'Tourists' = '#E50914')) +
  
  labs(x = 'Date', color = '') +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.05, 0.85),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = 'transparent', color = NA),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'in'),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = '#00A087', face = 'bold', margin = margin(r = 1)),
    axis.title.y.right = element_text(size = base.size * 0.9, family = base.family, color = '#E50914', face = 'bold', margin = margin(l = 1)),
    axis.text = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.4, l = 0.1, unit = 'in'),
  )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B: Cross-Correlation Function, CCF
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ts_tourist_diff <- diff(df_ts$Tourist)
ts_local_diff   <- diff(df_ts$Local)

ccf_res <- ccf(ts_tourist_diff, ts_local_diff, plot = FALSE, lag.max = 10, na.action = na.pass)

df_ccf <- data.frame(
  Lag = ccf_res$lag,
  ACF = ccf_res$acf
)

ci_threshold <- 1.96 / sqrt(length(ts_tourist_diff))

# Plot
plot_B <- ggplot(df_ccf, aes(x = Lag, y = ACF)) +
  geom_bar(
    stat = 'identity', 
    fill = ifelse(df_ccf$ACF > ci_threshold & df_ccf$Lag < 0, '#E50914', 'gray70'), 
    width = 0.6
  ) +
  geom_hline(yintercept = c(ci_threshold, -ci_threshold), linetype = 'dashed', color = 'blue') +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(
    x = 'Lag (Weeks)', 
    y = 'Cross-Correlation (Differenced)'
  ) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    panel.grid = element_blank(),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(r = 1)),
    axis.title.x.bottom = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5)),
    axis.text = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    axis.title.x.top = element_blank(),
    plot.margin = margin(t = 0, r = -1, b = 0.1, l = 0.1, unit = 'in')
  )
  
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C: Spatial Heatmap
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dat.pl <- dat.raw %>% dplyr::select(KID, Local, FLUA, FLUB, WGS84LNG, WGS84LAT) %>% 
  filter(Local != -1) %>% 
  filter(WGS84LNG != 'NONE' | WGS84LAT != 'NONE')

if (flu == 'FLUA') {
  PDF <- readxl::read_excel(glue('{tbl.dir}/Table.Influenza.Peaks.A.xlsx'))
} else {
  PDF <- readxl::read_excel(glue('{tbl.dir}/Table.Influenza.Peaks.B.xlsx'))
}
MACAU <- sf::st_read(glue('{dat.dir}/Macau/820000.geojson'))[c('name','geometry')]

dat.spatial <- data.frame()
for (index in 1:nrow(PDF)) {
  sdate = PDF[index, 'Start']$Start
  datea = sdate + 7 * 24 * 60 * 60
  dateb = sdate + 35 * 24 * 60 * 60
  
  tmpa <- dat.pl[KID > as.Date(sdate) & KID <= as.Date(datea) & get(flu) == 1, .(WGS84LNG, WGS84LAT)]
  tmpb <- dat.pl[KID > as.Date(datea) & KID <= as.Date(dateb) & get(flu) == 1, .(WGS84LNG, WGS84LAT)]
  tmpa$Phase <- 'Early Phase'
  tmpb$Phase <- 'Spread Phase'
  
  colnames(tmpa) <- c('lon', 'lat', 'Phase')
  colnames(tmpb) <- c('lon', 'lat', 'Phase')
  
  dat.spatial <- rbind(dat.spatial, tmpa, tmpb)
}
dat.spatial <- as.data.frame(dat.spatial)
dat.spatial$lon <- as.numeric(dat.spatial$lon)
dat.spatial$lat <- as.numeric(dat.spatial$lat)
dat.spatial <- dat.spatial[dat.spatial$lon < 114, ]

plot_C <- ggplot(dat.spatial, aes(x = lon, y = lat)) +
  ggplot2::geom_sf(data = MACAU, fill = 'NA', colour = '#444444', linewidth = 0.10, inherit.aes = FALSE) +
  stat_density_2d(aes(fill = ..level..), geom = 'polygon', color = 'white', alpha = 0.8, linewidth = 0.40) +
  scale_fill_viridis_c(option = 'magma') +
  facet_wrap(~Phase, strip.position = 'bottom') +
  annotate('text', x = 113.553889, y = 22.223556, label = 'Tourist Hub', color = '#E50914', size = base.size / 4.88) +
  annotate(
    'segment', 
    x = 113.553889, y = 22.220556,
    xend = 113.550889, yend = 22.210556,
    color = '#E50914',
    linewidth = 0.75,
    arrow = arrow(length = unit(0.2, 'cm'))
  ) +
  annotate('text', x = 113.583889, y = 22.190556, label = 'High-Density\nResidential Area', color = '#A1F7F3', size = base.size / 4.88, lineheight = 0.85) +
  annotate(
    'segment', 
    x = 113.583889, y = 22.195556,
    xend = 113.555889, yend = 22.208556,
    color = '#A1F7F3',
    linewidth = 0.75,
    arrow = arrow(length = unit(0.2, 'cm'))
  ) +
  annotate(
    'segment', 
    x = 113.583889, y = 22.195556,
    xend = 113.543889, yend = 22.203556,
    color = '#A1F7F3',
    linewidth = 0.75,
    arrow = arrow(length = unit(0.2, 'cm'))
  ) +
  annotate(
    'segment', 
    x = 113.583889, y = 22.185356,
    xend = 113.556889, yend = 22.159256,
    color = '#A1F7F3',
    linewidth = 0.75,
    arrow = arrow(length = unit(0.2, 'cm'))
  ) +
  labs(x = 'Spatial Diffusion Pattern', y = '') +
  theme_dark(base_size = base.size, base_family = base.family) +
  theme(
    axis.title.x = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 0.20, unit = 'in')),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_line(color = NA),
    legend.position = 'none',
    plot.margin = margin(t = 0, r = -1, b = 0.2, l = 0, unit = 'in')
  )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bottom_row <- plot_grid(
  plot_B, plot_C,
  ncol = 2,
  rel_widths = c(1, 1.4) 
)

gg <- plot_grid(
  plot_A, bottom_row,
  ncol = 1,
  rel_heights = c(1, 1.4) 
)

# Show plot
print(gg)

# Save to file
width = 12; height = 10
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{flu}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')