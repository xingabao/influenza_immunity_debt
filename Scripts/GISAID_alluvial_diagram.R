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
suppressMessages(suppressWarnings(library(ggalluvial)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
these <- c('B', 'H1N1', 'H3N2')

# Plotting Constants
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Load data
meta.file <- glue('{rt.dir}/data/HK/gisaid_epiflu_isolates.xls')

# Load meta data
dat.meta <- suppressWarnings(readxl::read_xls(meta.file)) %>%
  mutate(
    accession_id = str_split(`HA Segment_Id`, '\\|', simplify = TRUE)[, 1],
    Year = year(as.Date(Collection_Date))
  ) %>% dplyr::select('Year', 'accession_id', 'Collection_Date', 'Subtype', 'Clade', 'Isolate_Submitter')

openxlsx::write.xlsx(dat.meta, glue('{tbl.dir}/Table.Clades.xlsx'))

dat.nexts <- data.frame()
for (this in these) {
  next.file <- glue('{dat.dir}/HK/GISAID/nextclade.{this}.tsv')
  dat.next <- suppressWarnings(readr::read_tsv(next.file, show_col_types = FALSE))
  next. <- dat.next %>%
    mutate(
      accession_id = str_split(seqName, '\\|', simplify = TRUE)[, 1],
      Clade = replace_na(clade, 'Unknown'),
    ) %>% 
    dplyr::select(accession_id, Clade) %>%
    mutate(
      # Logical Explanation:
      # ^          : Matches the beginning of the string
      # [A-Za-z]+  : Matches one or more letters (e.g., 'B', 'J', 'unassigned')
      # (          : Starts a capturing group (optional part)
      # \\.        : Matches a literal dot (.)
      # \\d+       : Matches one or more digits
      # )?         : Ends the capturing group, where ? indicates this group appears 0 or 1 time
      Clade = str_extract(Clade, '^[A-Za-z]+(\\.\\d+)?')
    )
  dat.nexts <- rbind(dat.nexts, next.)
}

dat.nexts <- dat.nexts %>% filter(Clade != 'unassigned')
meta <- dat.nexts %>% 
  dplyr::left_join(dat.meta %>% dplyr::select(-Clade)) %>%
  mutate(
    Subtype = case_when(
      Subtype == 'B' ~ 'B',
      Subtype == 'A / H1N1' ~ 'A/H1N1pdm09',
      Subtype == 'A / H3N2' ~ 'A/H3N2',
      TRUE ~ Subtype
    ),
    Group = glue('{Subtype} ({Clade})')
  )

# Arrange data
dat.pl <- meta %>%
  count(Year, Group, Subtype) %>%
  group_by(Year) %>%
  mutate(Freq = n / sum(n)) %>%
  ungroup()

label.data <- dat.pl %>%
  arrange(Year, desc(Group)) %>%
  group_by(Year) %>%
  mutate(
    y_pos = cumsum(Freq) - Freq / 2 
  ) %>%
  ungroup() %>%
  group_by(Group) %>%
  filter(Year == min(Year)) %>%
  ungroup()

# Plot
gg <- ggplot(dat.pl, aes(x = Year, y = Freq, alluvium = Group, stratum = Group)) +
  geom_flow(aes(fill = Group), width = 1/4, alpha = 0.4, color = 'white', size = 0.1) +
  geom_stratum(aes(fill = Group), width = 1/4, alpha = 1, color = 'white', size = 0.2) +
  geom_text_repel(
    data = label.data,
    aes(
      x = Year, 
      y = y_pos, 
      label = Group
    ),
    size = 4.5, 
    color = '#000000',
    nudge_x = -0.3,
    direction = 'y',
    hjust = 1,
    family = base.family,
    segment.size = 0.2,
    segment.color = 'grey50',
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_continuous(breaks = unique(dat.pl$Year), expand = c(0, 0)) +
  theme_bw(base_size = base.size, base_family = base.family) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme(
    legend.position = 'none', 
    panel.grid = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = base.size * 0.9, family = base.family, color = base.col)
  )

# Show plot
print(gg)

# Save to file
width = 15; height = 4
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, units = 'in', bg = '#FFFFFF')