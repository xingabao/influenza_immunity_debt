# Load R packages
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(sp)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggspatial)))
suppressMessages(suppressWarnings(library(patchwork)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

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
  UDF <- dat.raw[dat.raw$WGS84LNG != 'NONE', ]

  IDF <- UDF %>% 
    dplyr::filter(FLUA != -1 & FLUB != -1) %>% 
    dplyr::filter(KID >= as.Date('2010-01-01')) %>% 
    dplyr::select(KID, FLUA, FLUB, WGS84LNG, WGS84LAT)
  
  saveRDS(IDF, file = glue('{dat.dir}/Macau/Kriging_based_spatial_interpolation.rds'))
} else {
  IDF <- readRDS(glue('{dat.dir}/Macau/Kriging_based_spatial_interpolation.rds'))
}

# 
MACAU <- sf::st_read(glue('{dat.dir}/Macau/820000.geojson'))[c('name','geometry')]

# 
MACAU.Z <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/澳门特别行政区.json'))
MACAU.A <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/嘉模堂区.json'))
MACAU.B <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/风顺堂区.json'))
MACAU.C <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/花地玛堂区.json'))
MACAU.D <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/大堂区.json'))
MACAU.E <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/花王堂区.json'))
MACAU.F <- sf::st_read(glue('{dat.dir}/Macau/澳门行政区/望德堂区.json'))
MACAU.G <- MACAU.Z %>% filter(name == '圣方济各堂区'); MACAU.G$shi <- '圣方济各堂区'; MACAU.G$name <- '路环区'
MACAU.H <- MACAU.Z %>% filter(name == '新城A区'); MACAU.H$shi <- '新城A区'; MACAU.H$name <- '新城A区及港珠澳大桥澳门口岸管理区'
MACAU.I <- MACAU.Z %>% filter(name == '路氹填海区'); MACAU.I$shi <- '路氹填海区'; MACAU.I$name <- '路氹填海区'
MACAU.J <- MACAU.Z %>% filter(name == '澳门大学'); MACAU.J$shi <- '澳门大学'; MACAU.J$name <- '澳门大学'
distrincts <- rbind(MACAU.A, MACAU.B, MACAU.C, MACAU.D, MACAU.E, MACAU.F, MACAU.G, MACAU.H, MACAU.I, MACAU.J) %>% dplyr::select(shi, name, lng, lat, geometry)
MACAU <- rbind(MACAU, MACAU.J %>% dplyr::select(name, geometry))

# 
flu.types = c("AB", "A", "B")

#
IDF_sf <- st_as_sf(
  IDF,
  coords = c("WGS84LNG", "WGS84LAT"),
  crs = 4326
)

#
IDF_with_region <- st_join(
  IDF_sf,
  distrincts[, c("name")],
  join = st_within
)

IDF_with_region <- IDF_with_region %>% filter(!is.na(name))

a.date <- '2020-04-01'
b.date <- '2021-03-31'
c.date <- '2023-04-30'

pre.pan <- IDF_with_region %>% dplyr::filter(KID < as.Date(b.date))
acu.pan <- IDF_with_region %>% dplyr::filter(KID >= as.Date(a.date) & KID < as.Date(b.date))
tra.pan <- IDF_with_region %>% dplyr::filter(KID >= as.Date(b.date) & KID < as.Date(c.date))
pos.pan <- IDF_with_region %>% dplyr::filter(KID >= as.Date(c.date))

# Define a general summarizing function
summarize.flu <- function(data, flu.type = flu.types) {
  flu.type <- match.arg(flu.type)
  data <- data %>%
    mutate(
      FLU = case_when(
        flu.type == "AB" ~ ifelse(FLUA == 1 | FLUB == 1, 1, 0),
        flu.type == "A"  ~ ifelse(FLUA == 1, 1, 0),
        flu.type == "B"  ~ ifelse(FLUB == 1, 1, 0)
      )
    ) %>%
    group_by(name) %>%
    summarise(
      Total_samples = n(),
      Positive_cases = sum(FLU),
      positive.rate = Positive_cases / nrow(data) * 100,
      .groups = "drop"
    ) %>% 
    as.data.frame() %>% 
    dplyr::select(name, positive.rate)
  return(data)
}

# Apply function for each type
cc <- distrincts %>% as.data.frame() %>% dplyr::select(name, lng, lat)
epi.dat.pre.A <- summarize.flu(pre.pan, "A") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.acu.A <- summarize.flu(acu.pan, "A") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.tra.A <- summarize.flu(tra.pan, "A") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.pos.A <- summarize.flu(pos.pan, "A") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.pre.B <- summarize.flu(pre.pan, "B") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.acu.B <- summarize.flu(acu.pan, "B") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.tra.B <- summarize.flu(tra.pan, "B") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.pos.B <- summarize.flu(pos.pan, "B") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.pre.AB <- summarize.flu(pre.pan, "AB") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.acu.AB <- summarize.flu(acu.pan, "AB") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.tra.AB <- summarize.flu(tra.pan, "AB") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)
epi.dat.pos.AB <- summarize.flu(pos.pan, "AB") %>% left_join(cc) %>% dplyr::select(lng, lat, positive.rate) %>% rename('longitude' = lng, 'latitude' = lat, 'value' = positive.rate)

if (TRUE) {
  epi.dat.pre.A$G1 <- 'A'; epi.dat.pre.A$G2 <- 'pre'
  epi.dat.acu.A$G1 <- 'A'; epi.dat.acu.A$G2 <- 'acu'
  epi.dat.tra.A$G1 <- 'A'; epi.dat.tra.A$G2 <- 'tra'
  epi.dat.pos.A$G1 <- 'A'; epi.dat.pos.A$G2 <- 'pos'
  epi.dat.pre.B$G1 <- 'B'; epi.dat.pre.B$G2 <- 'pre'
  epi.dat.acu.B$G1 <- 'B'; epi.dat.acu.B$G2 <- 'acu'
  epi.dat.tra.B$G1 <- 'B'; epi.dat.tra.B$G2 <- 'tra'
  epi.dat.pos.B$G1 <- 'B'; epi.dat.pos.B$G2 <- 'pos'
  epi.dat.pre.AB$G1 <- 'AB'; epi.dat.pre.AB$G2 <- 'pre'
  epi.dat.acu.AB$G1 <- 'AB'; epi.dat.acu.AB$G2 <- 'acu'
  epi.dat.tra.AB$G1 <- 'AB'; epi.dat.tra.AB$G2 <- 'tra'
  epi.dat.pos.AB$G1 <- 'AB'; epi.dat.pos.AB$G2 <- 'pos'
  
  epi.dat.merge <- rbind(
    epi.dat.pre.A, epi.dat.acu.A, epi.dat.tra.A, epi.dat.pos.A,
    epi.dat.pre.B, epi.dat.acu.B, epi.dat.tra.B, epi.dat.pos.B,
    epi.dat.pre.AB, epi.dat.acu.AB, epi.dat.tra.AB, epi.dat.pos.AB
  )
  
  range(epi.dat.merge$value)
  write.csv(epi.dat.merge, glue('{dat.dir}/Macau/epi.dat.merge.csv'), row.names = FALSE)
  openxlsx::write.xlsx(x = epi.dat.merge, file = glue('{tbl.dir}/Table.epi.dat.merge.xlsx'))
}

# Colors
pals <- c("#FFFFCC", "#FFFFCC", "#FFF4B2", "#FEEA9A", "#FEDE82", "#FECD6A", "#FEB752", "#FDA245", "#FD8D3C", "#FC6931", "#F84628", "#EA2820", "#D8121E", "#C20324", "#A20026", "#800026")

plot.gg <- function(dat, flu.type, time.type) {
  
  dat <- dat %>% filter(G1 == flu.type, G2 == time.type)
  
  print(glue('{flu.type}, {time.type}, {range(dat$value)[1]}, {range(dat$value)[2]}'))
  
  dat.sp <- dat
  
  sp::coordinates(dat.sp) <- ~longitude+latitude
  
  dat.sf <- sf::st_as_sf(dat, coords = c('longitude', 'latitude'), crs = 4326)
  
  # Generate a map grid for interpolation calculations
  ex. <- .00005
  map.box <- sf::st_bbox(MACAU)
  map.grid <- base::expand.grid(
    longitude = seq(
      from = map.box[1] * (1 - ex.),
      to = map.box[3] * (1 + ex.), length.out = 500),
    latitude = seq(
      from = map.box[2] * (1 - ex.), 
      to = map.box[4] * (1 + ex.), length.out = 500)
  )
  
  # Convert to spatial points grid in sp format
  coordinates(map.grid) <- ~longitude+latitude
  
  # Set to grid structure for spatial interpolation
  gridded(map.grid) <- TRUE
  
  # Spatial Variogram analysis to examine spatial correlation of observations
  variogram. <- gstat::variogram(object = value ~ 1, locations = dat.sf, data = dat.sf)
  if (FALSE) plot(variogram.)
  
  # Fit a Semivariogram model to determine spatial interpolation weights
  vgm. <- gstat::vgm(psill = 600, model = 'Lin', nugget = 0, range = 20)
  fit.variogram. <- suppressWarnings({ gstat::fit.variogram(object = variogram., model = vgm.) })
  if (FALSE) plot(variogram., fit.variogram.)
  
  # Kriging interpolation calculation
  krige.res <- gstat::krige(formula = value ~ 1, locations = dat.sp, newdata = map.grid, model = vgm.)
  
  # Convert Kriging interpolation results to data frame
  krige.dat <- as.data.frame(krige.res)
  
  colnames(krige.dat)[1:3] <- c('Longitude', 'Latitude', 'value')
  
  if (TRUE) {
    krige.dat.sf <- st_as_sf(x = krige.dat, coords = c('Longitude', 'Latitude'), crs = 4326)
    krige.dat.sf <- krige.dat.sf[st_make_valid(distrincts), ]
  }
  
  # Plot Kriging interpolation heatmap
  gg <- ggplot() +
    ggplot2::geom_sf(data = krige.dat.sf, mapping = aes(color = value)) +
    ggplot2::geom_sf(data = MACAU, fill = 'NA', colour = '#666666', linewidth = 0.10) +
    ggplot2::geom_sf(data = dat.sf, shape = 1, size = 3, color = '#666666') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(min(krige.dat$Longitude), max(krige.dat$Longitude), length.out = 5)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_gradientn(
      colors = pals,
      na.value = pals[1],
      breaks = seq(0, 1.5, 0.25),
      limits = c(0, 1.5)
    ) +
    coord_sf(clip = 'off') + 
    guides(color = guide_colorbar(title = NULL, title.position = "top", label.position = "bottom", title.hjust = 0.5, barwidth = 30)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'NA'),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = 'top',
      legend.direction = "horizontal",
      legend.title = element_text(family = base.family, color = base.col, size = base.size, margin = margin(b = 10)),
      legend.text = element_text(family = base.family, color = base.col, size = base.size, margin = margin(t = 10)),
      legend.background = element_blank(),
      legend.key.height = unit(.2, units = 'in'),
      legend.key.width  = unit(.6, units = 'in'),
      legend.key.size = unit(.6, units = 'in'),
      legend.frame = element_rect(color = '#FFFFFF'),
      legend.ticks = element_line(color = '#FFFFFF', linewidth = 2),
    )
  
  if (flu.type == 'AB' & time.type == 'pos') {
    gg  +
      annotation_north_arrow(style = north_arrow_nautical(text_size = base.size, text_family = base.family, text_col = base.col), location = 'tr', pad_x = unit(-0.50, 'cm')) +
      annotation_scale(location = 'br', text_cex = base.size / 12, text_family = base.family, line_width = 1.5, pad_x = unit(-0.50, 'cm'))
  } else {
    gg
  }
}

ALL <- list(); ind = 1
epi.dat.merge <- data.table::fread(glue('{dat.dir}/Macau/epi.dat.merge.csv'))
gg.pre.A <- plot.gg(epi.dat.merge, flu.type = 'A', time.type = 'pre'); ALL[[ind]] <- gg.pre.A; ind <- ind + 1
gg.acu.A <- plot.gg(epi.dat.merge, flu.type = 'A', time.type = 'acu'); ALL[[ind]] <- gg.acu.A; ind <- ind + 1
gg.tra.A <- plot.gg(epi.dat.merge, flu.type = 'A', time.type = 'tra'); ALL[[ind]] <- gg.tra.A; ind <- ind + 1
gg.pos.A <- plot.gg(epi.dat.merge, flu.type = 'A', time.type = 'pos'); ALL[[ind]] <- gg.pos.A; ind <- ind + 1
gg.pre.B <- plot.gg(epi.dat.merge, flu.type = 'B', time.type = 'pre'); ALL[[ind]] <- gg.pre.B; ind <- ind + 1
gg.acu.B <- plot.gg(epi.dat.merge, flu.type = 'B', time.type = 'acu'); ALL[[ind]] <- gg.acu.B; ind <- ind + 1
gg.tra.B <- plot.gg(epi.dat.merge, flu.type = 'B', time.type = 'tra'); ALL[[ind]] <- gg.tra.B; ind <- ind + 1
gg.pos.B <- plot.gg(epi.dat.merge, flu.type = 'B', time.type = 'pos'); ALL[[ind]] <- gg.pos.B; ind <- ind + 1
gg.pre.AB <- plot.gg(epi.dat.merge, flu.type = 'AB', time.type = 'pre'); ALL[[ind]] <- gg.pre.AB; ind <- ind + 1
gg.acu.AB <- plot.gg(epi.dat.merge, flu.type = 'AB', time.type = 'acu'); ALL[[ind]] <- gg.acu.AB; ind <- ind + 1
gg.tra.AB <- plot.gg(epi.dat.merge, flu.type = 'AB', time.type = 'tra'); ALL[[ind]] <- gg.tra.AB; ind <- ind + 1
gg.pos.AB <- plot.gg(epi.dat.merge, flu.type = 'AB', time.type = 'pos'); ALL[[ind]] <- gg.pos.AB

# Combination
gg <- wrap_plots(ALL, ncol = 4, guides = "collect") +
  plot_annotation(tag_levels = list(c('', '', '', '', '', '', '', '', '', '', '', ''))) &
  theme(legend.position = "bottom", plot.tag = element_text(family = base.family, color = base.col, size = base.size))

# Save to file
width = 11; height = 14
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, bg = '#FFFFFF')