# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale("LC_TIME", "C")

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
chart.plot <- function(cat){
  data <- IFL.DAT.M %>%
    filter(CAT == cat)
  data <- data %>% dplyr::select(CAT, date, positive.rate)
  min_axis <- ifelse(min(data$positive.rate) > 0, 0, min(data$positive.rate)) - 0
  max_axis <- ifelse(max(data$positive.rate) < 0, 0, max(data$positive.rate)) + 0.15
  min_date <- min(data$date)
  data$date <- as.Date(data$date)
  
  if (cat == 'UK') {
    cat. = 'United Kingdom'
  } else if (cat == 'Russia') {
    cat. = 'Russian Federation'
  } else {
    cat. = cat
  }
  
  gg <- data %>%
    ggplot(aes(date, positive.rate)) +
    geom_point(size = 2) +
    geom_line(color = "grey50") +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    annotate("rect", xmin = as.Date("2020-04-01"), xmax = as.Date("2021-03-31"), ymin = 0, ymax = max_axis, alpha = 0.2, fill = "#EE000055") +
    annotate("rect", xmin = as.Date("2021-04-01"), xmax = as.Date("2023-04-30"), ymin = 0, ymax = max_axis, alpha = 0.2, fill = "#41B34955") +
    scale_y_continuous(limits = c(min_axis, max_axis), labels = scales::percent) +
    scale_x_date(NULL, date_breaks = "1 year", date_labels = "%b\n%Y", expand = c(0.01, 0)) +
    guides(color = 'none') +
    labs(y = NULL, x = NULL, subtitle = cat.) +
    theme_bw(base_size = 13, base_family = base.family) +
    theme(
      axis.text = element_text(color = base.col, family = base.family),
      panel.grid.major.x = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA)
    )
  
  if (cat == 'Russia') {
    font.size <- 3
    lineheight <- 0.8
    lwd.size <- 0.38
    gg + annotate("text", x = c(as.Date("2018-08-01")), y = 0.65, size = font.size, label = c("Pre-pandemic\nperiod"), lineheight = lineheight) +
      geom_segment(aes(xend = as.Date("2020-04-01"), y = 0.65, x = as.Date("2021-03-01"), yend = 0.65), arrow = arrow(length = unit(0.1, 'cm')), lwd = lwd.size) +
      geom_segment(aes(xend = as.Date("2023-04-01"), y = 0.65, x = as.Date("2022-06-01"), yend = 0.65 ), arrow = arrow(length = unit(0.1, 'cm')), lwd = lwd.size) +
      annotate("text", x = c(as.Date("2021-10-01")), y = 0.65, size = font.size, label = c("Pandemic\nperiod"), lineheight = lineheight) +
      annotate("text", x = c(as.Date("2020-10-01")), y = 0.50, size = font.size, label = c("Acute\nphase"), color = "#EE0000", lineheight = lineheight) +
      annotate("text", x = c(as.Date("2022-04-01")), y = 0.50, size = font.size, label = c("Transition\nphase"), color = "#3B4992", lineheight = lineheight) +
      annotate("text", x = c(as.Date("2024-02-01")), y = 0.65, size = font.size, label = c("Post-\npandemic\nperiod"), lineheight = lineheight)
  } else {
    gg
  }
}

# Load data
if (!exists('dat.a')) { dat.a <- readxl::read_excel(glue('{dat.dir}/Global.Influenza-A.xlsx')) }
if (!exists('dat.b')) { dat.b <- readxl::read_excel(glue('{dat.dir}/Global.Influenza-B.xlsx')) }

# 'Azerbaijan', 'Poland', 'Germany',
# 'Russian Federation', 'France', 
# 'Kazakhstan', 'Lebanon', 'Lithuania', 
# 'Sweden', 'Slovakia', 'Turkey', 'Ukraine', 'Spain'
dat.raw <- rbind(dat.a, dat.b)
IFL.DAT <- dat.raw %>% 
  filter(ORIGIN_SOURCE == 'SENTINEL', ISO_YEAR > 2012, ISO_YEAR < 2025) %>% 
  distinct() %>% 
  group_by(`COUNTRY/AREA/TERRITORY`) %>%
  mutate(N = n()) %>%
  ungroup() %>%
  dplyr::select(date = ISO_SDATE, year = ISO_YEAR, week = ISO_WEEK, CAT = `COUNTRY/AREA/TERRITORY`, nFLUAB = INF_ALL, nSample = SPEC_PROCESSED_NB, nn = INF_NEGATIVE) %>%
  mutate(CAT = recode(
    CAT, 
    "Bolivia (Plurinational State of)" = "Bolivia",
    "Czechia" = "Czech Republic",
    "Lao People's Democratic Republic" = "Laos",
    "Netherlands (Kingdom of the)" = "Netherlands",
    "Republic of Korea" = "South Korea",
    "Republic of Moldova" = "Moldova",
    "Russian Federation" = "Russia",
    "Türkiye" = "Turkey",
    "C&#244;te d&#8217;Ivoire" = "Ivory Coast",
    "Iran (Islamic Republic of)" = "Iran",
    "United Republic of Tanzania" = "Tanzania",
    "Venezuela (Bolivarian Republic of)" = "Venezuela",
    "Viet Nam" = "Vietnam",
    "United States of America" = "USA"
  )) %>%
  mutate(
    nFLUAB = as.numeric(nFLUAB),
    nSample = as.numeric(nSample),
    nn = as.numeric(nn)
  ) %>%
  mutate(
    nFLUAB_imputed = case_when(
      is.na(nFLUAB) & !is.na(nSample) & !is.na(nn) ~ nSample - nn,
      TRUE ~ nFLUAB
    ),
    nSample_imputed = case_when(
      is.na(nSample) & !is.na(nFLUAB) & !is.na(nn) ~ nFLUAB + nn,
      TRUE ~ nSample
    )
  ) %>%
  mutate(
    nFLUAB = nFLUAB_imputed,
    nSample = nSample_imputed
  ) %>%
  filter(
    !is.na(nFLUAB), !is.na(nSample),
    nSample > 0,
    nFLUAB >= 0,
    nFLUAB <= nSample
  ) %>%
  mutate(
    nFLUAB = as.integer(nFLUAB),
    nSample = as.integer(nSample)
  ) %>% 
  dplyr::select(-nFLUAB, -nSample) %>%
  dplyr::select(date, year, week, CAT, nFLUAB = nFLUAB_imputed, nSample = nSample_imputed) %>%
  mutate(month = substr(date, 6, 7))

# Map
uk.regions <- c(
  "United Kingdom, England",
  "United Kingdom, Northern Ireland",
  "United Kingdom, Scotland",
  "United Kingdom, Wales"
)
IFL.DAT.H <- IFL.DAT %>%
  mutate(
    CAT = case_when(
      CAT %in% uk.regions ~ "UK",
      TRUE ~ CAT
    )
  ) %>%
  group_by(CAT) %>%
  summarise(
    Total_samples = sum(nSample, na.rm = TRUE),
    Positive_cases = sum(nFLUAB, na.rm = TRUE),
    positive.rate = Positive_cases / Total_samples * 100
  ) %>%
  filter(Total_samples > 1000) %>%
  ungroup()

# Calculate rate
IFL.DAT.M <- IFL.DAT %>%
  mutate(
    CAT = case_when(
      CAT %in% uk.regions ~ "UK",
      TRUE ~ CAT
    )
  ) %>%
  group_by(CAT, year, week) %>%
  summarise(
    Total_samples = sum(nSample, na.rm = TRUE),
    Positive_cases = sum(nFLUAB, na.rm = TRUE),
    positive.rate = Positive_cases / Total_samples,
    date = as.Date(min(date, na.rm = TRUE))
  ) %>%
  ungroup()

if (TRUE) {
  world_map <-  map_data("world") %>%
    filter(region != "Antarctica") %>%
    as_tibble() %>%
    fuzzyjoin::regex_left_join(maps::iso3166, c(region = "mapname"))
  C = 0; D = 0
  for (region in IFL.DAT.H$CAT) {
    if (region %in% unique(world_map$region)) {
      C = C + 1
    } else {
      print(region)
      D = D + 1
    }
  }
}

# 960 670 1433
map.dat <- map_data("world") %>%
  filter(region != "Antarctica") %>%
  as_tibble() %>%
  fuzzyjoin::regex_left_join(maps::iso3166, c(region = "mapname"))

world.map <- map.dat %>%
  left_join(IFL.DAT.H, by = c(region = "CAT")) %>%
  dplyr::select(long, lat, group, positive.rate)

map.lab <- map.dat %>% dplyr::select(long, lat, region) %>% group_by(region) %>% summarise(
  long = mean(long),
  lat = mean(lat)
)

# Draw plot
countries.lines <- data.frame()
lithuania.line <- tibble(x = 50.000000, xend = 25.279800, y = 85.000000, yend = 54.689160); countries.lines <- rbind(countries.lines, lithuania.line)
russia.line <- tibble(x = 178.000000, xend = 120.618423, y = 75.000000, yend = 65.751244); countries.lines <- rbind(countries.lines, russia.line)
ukraine.line <- tibble(x = 170.000000, xend = 170.000000, y = 140.000000, yend = 80.000000); countries.lines <- rbind(countries.lines, ukraine.line)
ukraine.line <- tibble(x = 170.000000, xend = 30.523333, y = 80.000000, yend = 50.450001); countries.lines <- rbind(countries.lines, ukraine.line)
kazakhstan.line <- tibble(x = 160.000000, xend = 71.449074, y = 30.000000, yend = 51.169392); countries.lines <- rbind(countries.lines, kazakhstan.line)
poland.line <- tibble(x = -18.000000, xend = -18.000000, y = 150.000000, yend = 100.000000); countries.lines <- rbind(countries.lines, poland.line)
poland.line <- tibble(x = -18.000000, xend = 21.01182, y = 100.000000, yend = 52.237049); countries.lines <- rbind(countries.lines, poland.line)
france.line <- tibble(x = -183.000000, xend = 2.348800, y = 55.000000, yend = 48.853410); countries.lines <- rbind(countries.lines, france.line)
uk.line <- tibble(x = -228.000000, xend = -0.09184, y = 108.000000, yend = 51.51279); countries.lines <- rbind(countries.lines, uk.line)
slovakia.line <- tibble(x = -180.000000, xend = 17.106740, y = -55.000000, yend = 48.148598); countries.lines <- rbind(countries.lines, slovakia.line)
spain.line <- tibble(x = -169.000000, xend = -3.703790, y = 0.000000, yend = 40.416775); countries.lines <- rbind(countries.lines, spain.line)
azerbaijan.line <- tibble(x = 180.000000, xend = 49.867092, y = -60.000000, yend = 40.409264); countries.lines <- rbind(countries.lines, azerbaijan.line)
germany.line <- tibble(x = -100.000000, xend = 9.735603, y = 90.000000, yend = 52.373920); countries.lines <- rbind(countries.lines, germany.line)
turkey.line <- tibble(x = -30.000000, xend = 32.866287, y = -55.000000, yend = 39.925533); countries.lines <- rbind(countries.lines, turkey.line)
lebanon.line <- tibble(x = 60.000000, xend = 33.893791, y = -55.000000, yend = 35.501777); countries.lines <- rbind(countries.lines, lebanon.line)

ins <- c('China')
map <- ggplot() +
  geom_polygon(data = world.map, aes(long, lat, group = group, fill = positive.rate), color = "#000000", linewidth = 0.3) +
  geom_segment(data = countries.lines, aes(x = x, xend = xend, y = y , yend = yend), color = "grey50", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-350, 350), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-120, 190)) +
  labs(fill = "Infection Rate") +
  geom_text(data = map.lab %>% filter(region %in% ins), aes(x = long, y = lat, label = region), size = 3, color = "black", check_overlap = TRUE) +
  guides(fill = guide_colorbar(title.position = "top", label.position = "bottom", title.hjust = 0.5, barwidth = 30)) +
  scale_fill_gradientn(
    colors = c("#F9E9CD", "#04A3FF", "#FF349C", "#EE3F4D"),
    na.value = "#FFFFFF",
    breaks = seq(0, 50, 5),
    limits = c(0, 50),
    labels = sprintf("%.0f%%", round(seq(0, 50, 5), 1))
  ) +
  theme_void() +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.17, 0.955),
    legend.direction = "horizontal",
    legend.title = element_text(size = 14, margin = margin(b = 10)),
    legend.text = element_text(size = 10),
    legend.background = element_blank(),
    legend.key.height = unit(.2, units = 'in'),
    legend.key.width  = unit(.7, units = 'in'),
    legend.key.size = unit(.7, units = 'in'),
    legend.frame = element_rect(colour = '#FFFFFF'),
    legend.ticks = element_line(colour = '#FFFFFF', size = 2),
  )

germany <- chart.plot('Germany')
russia <- chart.plot('Russia')
turkey <- chart.plot('Turkey')
uk <- chart.plot('UK')
spain <- chart.plot('Spain')
france <- chart.plot('France')
poland <- chart.plot('Poland')
azerbaijan <- chart.plot('Azerbaijan')
kazakhstan <- chart.plot('Kazakhstan')
lebanon <- chart.plot('Lebanon')
lithuania <- chart.plot('Lithuania')
sweden <- chart.plot('Sweden')
slovakia <- chart.plot('Slovakia')
ukraine <- chart.plot('Ukraine')

a <- function(x) { return(x + 0.20) }
aA <- function(x) { return(x + 0.22) }
aa <- function(x) { return(x + 0.24) }
A <- function(x) { return(x + 0.25) }
Aa <- function(x) { return(x + 0.27) }
b <- function(x) {  return(x + 0.17) }
B <- function(x) {  return(x + 0.20) }
Bb <- function(x) {  return(x + 0.22) }
BB <- function(x) {  return(x + 0.25) }

gg <- map +
  inset_element(uk, -0.003, 0.72, aA(-0.003), B(0.72)) +
  inset_element(lithuania, 0.48, 0.645, A(0.48), B(0.645)) +
  inset_element(turkey, 0.24, 0.00, A(0.24), BB(0.00)) +
  inset_element(ukraine, 0.70, 0.81, A(0.70), B(0.81)) +
  inset_element(russia, 0.75, 0.6, A(0.75), B(0.6)) +
  inset_element(kazakhstan, 0.73, 0.34, Aa(0.73), Bb(0.34)) +
  inset_element(france, -0.007, 0.50, aa(-0.007), B(0.50)) +
  inset_element(spain, 0.01, 0.28, A(0.01), Bb(0.28)) +
  inset_element(germany, 0.223, 0.67, aa(0.223), B(0.67)) +
  inset_element(slovakia, 0.00, 0.01, aa(0.00), BB(0.01)) +
  inset_element(azerbaijan, 0.75, 0.06, A(0.75), BB(0.06)) +
  inset_element(poland, 0.36, 0.84, A(0.36), b(0.84)) +
  inset_element(lebanon, 0.50, 0.00, A(0.50), BB(0.00))

# Save to file
width = 20; height = 12
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')