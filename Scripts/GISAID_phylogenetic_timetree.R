# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(treeio)))
suppressMessages(suppressWarnings(library(ggtree)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(lubridate)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'H1N1'  # Options: 'H1N1' or 'H3N2' or 'B'
if (this == 'H1N1') {
  this.label = 'A/H1N1pdm09'
} else if (this == 'H3N2') {
  this.label = 'A/H3N2'
} else if (this == 'B') {
  this.label = 'B'
}

# Plotting Constants
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Load data
tree.file <- glue('{dat.dir}/HK/GISAID/{this}/treetime/timetree.nexus') 
meta.file <- glue('{dat.dir}/HK/gisaid_epiflu_isolates.xls')
next.file <- glue('{dat.dir}/HK/GISAID/nextclade.{this}.tsv')

# Read TreeTime Nexus Tree and Output Dates and Node Confidence
tree <- read.beast(tree.file)

# Load meta data
dat.meta <- suppressWarnings(readxl::read_xls(meta.file))
dat.next <- suppressWarnings(readr::read_tsv(next.file, show_col_types = FALSE))

# Arrange data
meta. <- dat.meta %>%
  mutate(
    accession_id = str_split(`HA Segment_Id`, "\\|", simplify = TRUE)[, 1],
    Collection_Date = suppressWarnings(ymd(Collection_Date)), 
    Year = year(Collection_Date),
    Clade = replace_na(Clade, "Unknown"),
    Subtype = replace_na(Subtype, "Unknown")
  ) %>%
  filter(Year >= 2010) %>% 
  dplyr::select(accession_id, Collection_Date, Year, Subtype, Clade)

next. <- dat.next %>%
  mutate(
    accession_id = str_split(seqName, "\\|", simplify = TRUE)[, 1],
    Clade = replace_na(clade, "Unknown"),
  ) %>% 
  dplyr::select(accession_id, Clade)

meta <- next. %>% left_join(meta. %>% dplyr::select(-Clade))
meta <- meta %>%
  mutate(
    # Logical Explanation:
    # ^          : Matches the beginning of the string
    # [A-Za-z]+  : Matches one or more letters (e.g., 'B', 'J', 'unassigned')
    # (          : Starts a capturing group (optional part)
    # \\.        : Matches a literal dot (.)
    # \\d+       : Matches one or more digits
    # )?         : Ends the capturing group, where ? indicates this group appears 0 or 1 time
    Clade = str_extract(Clade, "^[A-Za-z]+(\\.\\d+)?")
  )

# Associate Metadata with Tree Object
tree.tbl <- as_tibble(tree) %>%
  mutate(
    accession_id = str_split(label, "\\|", simplify = TRUE)[, 1]
  )

tree_data_joined <- left_join(tree.tbl, meta, by = "accession_id")
tree.obj <- as.treedata(tree_data_joined)

max.date <- max(meta$Collection_Date, na.rm = TRUE)
x.limits <- c(2009, year(max.date) + 2)

p_base <- ggtree(tree.obj, mrsd = max.date) 
plot_data <- p_base$data
label.df <- plot_data %>%
  filter(isTip) %>%
  filter(Year >= 2023 | (Year == 2019 & Collection_Date > '2019-07-01')) %>%
  filter(Clade != "Unknown") %>%
  group_by(Clade) %>%
  slice_max(x, n = 1, with_ties = FALSE) %>%
  ungroup()

ignore.clades <- c("unassigned", "Unknown", NA)

# Plot
gg <- ggtree(tree.obj, mrsd = max.date) + 
  theme_tree2() + 
  
  geom_tippoint(
    data = . %>% filter(!Clade %in% ignore.clades & !is.na(Clade)),
    aes(color = Clade, shape = Subtype), 
    size = 2.5, 
    alpha = 0.8
  ) +
  
  geom_label_repel(
    data = label.df,
    aes(x = x, y = y, label = Clade, color = Clade),
    size = 4,
    family = base.family,
    nudge_x = 0.5,
    box.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE,
    fill = "white",
    alpha = 0.9
  ) +
  
  scale_x_continuous(
    limits = x.limits, 
    breaks = seq(2009, 2026, 2)
  ) +
  
  ggsci::scale_color_igv(name = "Clade") + 
  
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 8, 4), name = "Subtype") +
  guides(shape = 'none') +
  
  annotate('rect', xmin = 2020, xmax = 2023, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = 'gray30') +
  annotate('text', x = ifelse(this == 'B', 2016, ifelse(this == 'H3N2', 2015.5, 2016)), y = ifelse(this == 'B', 600, ifelse(this == 'H3N2', 950, 600)), label = ifelse(this == 'H3N2', 'Pre-Pandemic\nPeriod', 'Pre-Pandemic Period'), color = '#00A087', size = base.size / 4.00, family = base.family) +
  annotate('text', x = ifelse(this == 'H1N1', 2022, 2021.5), y = 50, label = 'NPIs / COVID-19 Pandemic', color = base.col, size = base.size / 4.00, family = base.family) +
  annotate('text', x = 2025.5, y = ifelse(this == 'B', 200, ifelse(this == 'H3N2', 400, 300)), label = 'Post-Pandemic\nPeriod', color = '#E50914', size = base.size / 4.00, family = base.family, lineheight = 0.85) +

  labs(subtitle = glue('Influenza {this.label}'), x = NULL) +
  
  theme(
    legend.position = "inside",
    legend.position.inside = c(ifelse(this == 'B', 0.05, ifelse(this == 'H3N2', 0.10, 0.05)), c(ifelse(this == 'B', 0.65, ifelse(this == 'H3N2', 0.65, 0.60)))),
    legend.title = element_blank(),
    legend.text = element_text(size = base.size * 0.7, family = base.family, color = base.col),
    legend.background = element_blank(),
    text = element_text(size = base.size, family = base.family, color = base.col),
    axis.title.x = element_text(size = base.size * 0.8, family = base.family, color = base.col, face = 'bold', margin = margin(t = 10)),
    axis.text.x = element_text(size = base.size * 0.8, family = base.family, color = base.col),
    plot.subtitle = element_text(size = base.size * 0.8, family = base.family, color = base.col, hjust = 0.5),
    axis.line.x = element_line(color = "black")
  )

# Show plot
print(gg)

# Save to file
width = 5; height = 4
ggsave(gg, filename = glue('{fig.dir}/{ofig}_{this}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')