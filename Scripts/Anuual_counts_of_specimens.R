# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggchicklet)))
suppressMessages(suppressWarnings(library(hrbrthemes)))

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
  IDF <- dat.raw %>% dplyr::filter(FLUA != -1 & FLUB != -1) %>% 
    dplyr::filter(KID >= as.Date('2010-01-01')) %>% 
    dplyr::select(Year, FLUA, FLUB) %>%
    mutate(Group = case_when(
      FLUA == 1 ~ 'Influenza A',
      FLUB == 1 ~ 'Influenza B',
      TRUE ~ 'Non-Influenza'
    )) %>%
    count(Year, Group) %>%
    rename(Count = n)
  
  saveRDS(IDF, file = glue('{dat.dir}/Macau/Anuual_counts_of_specimens.rds'))
} else {
  IDF <- readRDS(glue('{dat.dir}/Macau/Anuual_counts_of_specimens.rds'))
}

dat <- IDF %>% 
  mutate(Group = factor(Group, levels = c('Influenza A', 'Influenza B', 'Non-Influenza'))) %>% 
  mutate(Year = factor(Year))

if (TRUE) {
  dat.wide <- dat %>%
    mutate(Year = Year) %>% 
    tidyr::pivot_wider(
      id_cols = Year,
      names_from = Group,
      values_from = Count,
      values_fill = 0 
    ) %>%
    arrange(Year)
  
  openxlsx::write.xlsx(dat.wide, glue('{tbl.dir}/Table.Sample.Size.xlsx'))
}

color.mapping <- c(
  'Influenza A' = '#E50914',
  'Influenza B' = '#00A087',
  'Non-Influenza' = '#000000'
)

labels <- dat %>%
  arrange(Year, Group) %>%
  group_by(Year) %>%
  mutate(
    TotalCount = sum(Count),
    Percentage = Count / TotalCount,
    LabelPart = sprintf(
      '<span style="color:%s;">%s (%s)</span>',
      color.mapping[Group],
      scales::label_comma()(Count),
      scales::percent(Percentage, accuracy = 0.1)
    )
  ) %>%
  summarise(
    Label = paste(LabelPart, collapse = ' / '),
    Position = first(TotalCount)
  ) %>%
  mutate(
    Position = case_when(
      Position > 130000 ~ Position - 110000,
      TRUE ~ Position
    )
  )

# Draw Plot
gg <- dat %>%
  ggplot(aes(Year, Count, fill = Group)) +
  geom_chicklet(width = 0.75) +
  scale_y_continuous(
    limits = c(0, 220000),
    expand = c(0, 0.0625), position = 'right',
    breaks = seq(0, 210000, 30000), labels = scales::label_comma() 
  ) +
  geom_richtext(
    data = labels,
    aes(x = Year, y = Position, label = Label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.6,
    nudge_y = 1000,
    size = base.size / 3.88,
    label.color = NA, 
    fill = NA
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      'Influenza A' = '#AE4544',
      'Influenza B' = '#8DBCB8',
      'Non-Influenza' = '#CCCCCC'
    ),
    breaks = setdiff(unique(dat$Group), 'Other')
  ) +
  guides(fill = guide_legend(nrow = 3)) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_ipsum_rc(grid = 'X', base_family = base.family) +
  theme_classic(base_family = base.family, base_size = base.size * 0.9) +
  theme(
    axis.text.x = element_text(color = base.col, family = base.family, size = base.size * 0.9),
    axis.text.y = element_text(color = base.col, family = base.family, size = base.size * 0.9),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.2),
    legend.direction = 'horizontal',
    legend.text = element_text(family = base.family, color = base.col, size = base.size * 0.8),
    legend.background = ggfun::element_roundrect(color = '#636363', linetype = 1)
  )

# Save to file
width = 10; height = 6
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, bg = '#FFFFFF')