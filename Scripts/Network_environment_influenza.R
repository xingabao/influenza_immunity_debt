# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(qgraph)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(bootnet)))
suppressMessages(suppressWarnings(library(networktools)))
suppressMessages(suppressWarnings(library(NetworkComparisonTest)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'Macau'  # Options: 'HK' or 'Macau'
base.size <- 16
base.family <- 'serif'
base.col <- '#000000'

# Load data
dat.raw <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx'))
dat.raw$dif.temp <- dat.raw$max.temp - dat.raw$min.temp

# Prepare plot data
sdate <- '2020-01-01'
edate <- '2023-01-01'

dat.pl <- dat.raw %>%
  arrange(date) %>%
  mutate(
    Period = case_when(
      date < sdate ~ 'Pre-COVID',
      date > edate ~ 'Post-COVID',
      TRUE ~ 'COVID'
    )
  )

if (this == 'Macau') { dat.pl <- dat.pl %>% filter(date > '2014-12-31') }

data_pre  <- dat.pl %>% filter(Period == 'Pre-COVID') %>% dplyr::select(FLUA, FLUB, PM10, PM2.5, SO2, NO2, O3, CO, dif.temp, ave.temp, dew.temp, pressure, humidity, sunshine, wind.speed, precipitation)
data_post <- dat.pl %>% filter(Period == 'Post-COVID') %>% dplyr::select(FLUA, FLUB, PM10, PM2.5, SO2, NO2, O3, CO, dif.temp, ave.temp, dew.temp, pressure, humidity, sunshine, wind.speed, precipitation)

data_pre <- data_pre %>% na.omit()
data_post <- data_post %>% na.omit()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('Estimating network structure ...\n')

# Update group definitions
group.list <- list(       
  'Meteo' = c('dif.temp', 'min.temp', 'pressure', 'ave.temp', 'dew.temp', 'humidity', 'sunshine', 'wind.speed', 'precipitation'),
  'Pollution' = c('PM2.5', 'PM10', 'SO2', 'NO2', 'CO', 'O3'),
  'Outcome' = c('FLUA', 'FLUB')
)

nw.pre <- suppressWarnings(estimateNetwork(data_pre, default = 'EBICglasso', corMethod = 'cor_auto', tuning = 0.5))
nw.post <- suppressWarnings(estimateNetwork(data_post, default = 'EBICglasso', corMethod = 'cor_auto', tuning = 0.5))

# Define custom labels
label_map <- c(
  'FLUA'          = 'Influenza A',
  'FLUB'          = 'Influenza B',
  'PM10'          = 'PM10',
  'PM2.5'         = 'PM2.5',
  'SO2'           = 'SO2',
  'NO2'           = 'NO2',
  'O3'            = 'O3',
  'CO'            = 'CO',
  'dif.temp'      = 'Temp. Diff.',
  'ave.temp'      = 'Ave. Temp.',
  'dew.temp'      = 'Dew Temp.',
  'pressure'      = 'Pressure',
  'humidity'      = 'Humidity',
  'sunshine'      = 'Sunshine',
  'wind.speed'    = 'Wind Speed',
  'precipitation' = 'Precipitation'
)

orig_vars <- nw.pre$labels 

# Generate the final label vector based on the original variable order
final_labels <- orig_vars
match_idx <- match(orig_vars, names(label_map))
final_labels[!is.na(match_idx)] <- label_map[orig_vars[!is.na(match_idx)]]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cat('Generating custom central layout ...\n')

# Recalculate the layout (maintaining the previous central layout logic)
center_nodes <- c('FLUA', 'FLUB')
outer_nodes <- setdiff(orig_vars, center_nodes)
Layout <- matrix(0, nrow = length(orig_vars), ncol = 2)
rownames(Layout) <- orig_vars

# Set the center
Layout['FLUA', ] <- c(-0.2, 0) 
Layout['FLUB', ] <- c(0.2, 0)

n_outer <- length(outer_nodes)
angles <- seq(0, 2 * pi, length.out = n_outer + 1)[1:n_outer]
for (i in 1:n_outer) {
  Layout[outer_nodes[i], ] <- c(cos(angles[i]), sin(angles[i]))
}

real_vars <- nw.pre$labels 

target_groups <- list(
  'Meteo'     = group.list[['Meteo']],
  'Pollution' = group.list[['Pollution']],
  'Outcome'   = group.list[['Outcome']]
)

group.list.index <- lapply(target_groups, function(x) {
  as.vector(na.omit(match(x, real_vars)))
})

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
net.args <- list(
  layout = Layout, 
  groups = group.list.index,
  labels = final_labels,
  color = c('#00A08731', '#FF990033', '#E5091432'), 
  borders = FALSE,
  shape = 'ellipse',
  vsize = 13,
  vsize2 = 6,
  details = FALSE,
  label.cex = 1.1,
  label.scale = FALSE,
  legend = FALSE,
  label.color = 'black',
  posCol = '#E50914',
  negCol = '#00A087',
  mar = c(1, 2.0, 1, 2.0),
  rescale = TRUE,
  normalize = FALSE, 
  title = NULL
)

pdf(glue('{fig.dir}/{ofig}_{this}_pre.pdf'), width = 6, height = 6, family = 'serif', bg = '#FFFFFF')
do.call(qgraph, c(list(input = nw.pre$graph), net.args))
dev.off()

pdf(glue('{fig.dir}/{ofig}_{this}_post.pdf'), width = 6, height = 6, family = 'serif', bg = '#FFFFFF')
do.call(qgraph, c(list(input = nw.post$graph), net.args))
dev.off()

# ==============================================================================
# Bridge Strength
# ==============================================================================
cat('\nCalculating bridge centrality...\n')

W_pre <- nw.pre$graph
W_post <- nw.post$graph
node_names <- colnames(W_pre)

communities_vec <- rep(NA, length(node_names))
names(communities_vec) <- node_names
for (g_name in names(group.list)) {
  vars_in_group <- group.list[[g_name]]
  idx <- which(node_names %in% vars_in_group)
  communities_vec[idx] <- g_name
}
if (any(is.na(communities_vec))) communities_vec[is.na(communities_vec)] <- 'Unknown'

bridge_pre <- bridge(W_pre, communities = communities_vec)
bridge_post <- bridge(W_post, communities = communities_vec)

get_bs <- function(bridge_obj, node) bridge_obj$`Bridge Strength`[node]

res_df <- data.frame(
  Virus = c('Influneza A', 'Influneza A', 'Influneza B', 'Influneza B'),
  Scenario = c('Pre-Pandemic', 'Post-Pandemic', 'Pre-Pandemic', 'Post-Pandemic'),
  BridgeStrength = c(
    get_bs(bridge_pre, 'FLUA'),
    get_bs(bridge_post, 'FLUA'),
    get_bs(bridge_pre, 'FLUB'),
    get_bs(bridge_post, 'FLUB')
  )
)

cat('------------------------------------------------\n')
cat('Bridge Strength:\n')
print(res_df)
cat('------------------------------------------------\n')

# ==============================================================================
# NCT
# ==============================================================================
cat('Performing Network Comparison Test (NCT)...\n')
nct_res <- NCT(nw.pre, nw.post, it = 1000, test.edges = TRUE, progressbar = FALSE, verbose = FALSE)

# Output the statistical test results
cat('\nStatistical Test Results:\n')
cat('Overall Network Structure Difference (M) P-value:', nct_res$nwinv.pval, '\n')
cat('Overall Network Strength Difference (S) P-value:', nct_res$glstrinv.pval, '\n')

# ==============================================================================
# Plot
# ==============================================================================
res_df$Scenario <- factor(res_df$Scenario, levels = c('Pre-Pandemic', 'Post-Pandemic'))

# Bridge Strength (Sensitivity to Environment)
gg <- ggplot(res_df, aes(x = Virus, y = BridgeStrength, fill = Scenario)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf('%.3f', BridgeStrength)), position = position_dodge(width = 0.8), vjust = -0.5) +
  scale_fill_manual(values = c('#E50914', '#00A087'), breaks = c('Post-Pandemic', 'Pre-Pandemic')) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.3), expand = c(0, 0)) + 
  labs(x = NULL, y = 'Bridge Strength', subtitle = 'Sensitivity to environment') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(size = base.size, family = base.family, color = base.col),
    legend.position = 'inside',
    legend.position.inside = c(0.35, 0.88),
    legend.background = element_rect(fill = '#FFFFFF44'),
    legend.title = element_blank(),
    legend.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    axis.title.y.left = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(r = 5)),
    axis.title.x.bottom = element_text(size = base.size * 0.9, family = base.family, color = base.col, face = 'bold', margin = margin(t = 5)),
    axis.text = element_text(size = base.size * 0.9, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.9, family = base.family, color = base.col)
  )

print(gg)

# Save to file
ggsave(filename = glue('{fig.dir}/{ofig}_{this}.pdf'), plot = gg, width = 3.5, height = 4, units = 'in')