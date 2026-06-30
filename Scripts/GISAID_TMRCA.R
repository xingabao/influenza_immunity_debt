# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(treeio)))
suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(purrr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtree)))
suppressMessages(suppressWarnings(library(lubridate)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
these <- c('B', 'H1N1', 'H3N2')

# Period definitions
PRE_FROM      <- 2017
PRE_TO        <- 2020
NPI_FROM      <- 2020
NPI_TO        <- 2023
POST_FROM     <- 2023
MAX_SAMPLE    <- 60

# Plot constants
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

col.pre <- '#56B4E9'
col.post <- '#CC79A7'
col.cross <- '#D55E00'

comp_colors <- c(
  'Pre vs Pre' = col.pre,
  'Pre vs Post' = col.cross,
  'Post vs Post' = col.post
)

# Helpers
# TMRCA from tip years and patristic distance: (year_a + year_b - dist) / 2
calc_tmrca <- function(year_a, year_b, patristic) {
  (year_a + year_b - patristic) / 2
}

# Extract pairwise comparisons from patristic matrix
extract_pairs <- function(labels_a, labels_b, comparison, symmetric = FALSE, pat_mat, year_vec, clade_vec) {
  pairs <- expand.grid(
    label_a = labels_a,
    label_b = labels_b,
    stringsAsFactors = FALSE
  )
  if (symmetric) {
    pairs <- pairs %>% filter(label_a < label_b)
  } else {
    pairs <- pairs %>% filter(label_a != label_b)
  }
  
  pairs %>%
    mutate(
      patristic  = pat_mat[cbind(label_a, label_b)],
      year_a     = year_vec[label_a],
      year_b     = year_vec[label_b],
      clade_a    = clade_vec[label_a],
      clade_b    = clade_vec[label_b],
      tmrca      = calc_tmrca(year_a, year_b, patristic),
      comparison = comparison
    ) %>%
    filter(!is.na(patristic), !is.na(tmrca), tmrca >= 2005, tmrca <= 2026)
}

# Per-subtype analysis
analyze_subtype <- function(this) {
  
  tree.file <- glue('{dat.dir}/HK/GISAID/{this}/treetime/timetree.nexus')
  meta.file <- glue('{dat.dir}/HK/gisaid_epiflu_isolates.xls')
  next.file <- glue('{dat.dir}/HK/GISAID/nextclade.{this}.tsv')
  
  tree     <- read.beast(tree.file)
  dat.meta <- suppressWarnings(readxl::read_xls(meta.file))
  dat.next <- suppressWarnings(readr::read_tsv(next.file, show_col_types = FALSE))
  
  # Clean metadata
  meta. <- dat.meta %>%
    mutate(
      accession_id = str_split(`HA Segment_Id`, '\\|', simplify = TRUE)[, 1],
      Collection_Date = suppressWarnings(ymd(Collection_Date)),
      Year = year(Collection_Date),
      Clade = replace_na(Clade, 'Unknown'),
      Subtype = replace_na(Subtype, 'Unknown')
    ) %>%
    filter(Year >= 2010) %>%
    dplyr::select(accession_id, Collection_Date, Year, Subtype, Clade)
  
  # Merge Nextclade clade calls
  next. <- dat.next %>%
    mutate(
      accession_id = str_split(seqName, '\\|', simplify = TRUE)[, 1],
      Clade = replace_na(clade, 'Unknown')
    ) %>%
    dplyr::select(accession_id, Clade) %>%
    mutate(Clade = str_extract(Clade, '^[A-Za-z]+(\\.\\d+)?'))
  
  meta <- next. %>%
    left_join(meta. %>% dplyr::select(-Clade), by = 'accession_id')
  
  max.date <- max(meta$Collection_Date, na.rm = TRUE)
  
  # Attach metadata to tree
  tree.tbl <- as_tibble(tree) %>%
    mutate(accession_id = str_split(label, '\\|', simplify = TRUE)[, 1])
  
  tree_joined <- left_join(tree.tbl, meta, by = 'accession_id')
  tree.obj <- as.treedata(tree_joined)
  
  p_base <- ggtree(tree.obj, mrsd = max.date)
  plot_data <- p_base$data
  
  # Tip labels and periods
  tip_info <- plot_data %>%
    filter(isTip) %>%
    dplyr::select(node, label, x, accession_id, Clade, Year) %>%
    rename(tip_year = x) %>%
    mutate(
      Period = case_when(
        tip_year >= POST_FROM ~ 'Post-Pandemic',
        tip_year >= NPI_FROM  ~ 'Pandemic',
        tip_year >= PRE_FROM  ~ 'Pre-Pandemic',
        TRUE                  ~ 'Early'
      )
    )
  
  phylo_obj <- as.phylo(tree)
  n_tips <- length(phylo_obj$tip.label)
  
  pre_labels_all <- tip_info %>%
    filter(Period == 'Pre-Pandemic') %>%
    pull(label) %>%
    intersect(phylo_obj$tip.label)
  
  post_labels_all <- tip_info %>%
    filter(Period == 'Post-Pandemic') %>%
    pull(label) %>%
    intersect(phylo_obj$tip.label)
  
  if (length(pre_labels_all) < 2 || length(post_labels_all) < 2) {
    warning(sprintf('[%s] Insufficient tips, skipping ...', this))
    return(NULL)
  }
  
  # Downsample to MAX_SAMPLE
  set.seed(2024)
  pre_samp  <- if (length(pre_labels_all)  > MAX_SAMPLE) 
    sample(pre_labels_all,  MAX_SAMPLE) else pre_labels_all
  post_samp <- if (length(post_labels_all) > MAX_SAMPLE)
    sample(post_labels_all, MAX_SAMPLE) else post_labels_all
  
  all_samp <- unique(c(pre_samp, post_samp))
  
  sub_phylo     <- keep.tip(phylo_obj, all_samp)
  patristic_mat <- cophenetic.phylo(sub_phylo)
  
  valid_tips <- tip_info %>%
    filter(label %in% rownames(patristic_mat)) %>%
    distinct(label, .keep_all = TRUE)
  
  year_vec  <- setNames(valid_tips$tip_year, valid_tips$label)
  clade_vec <- setNames(valid_tips$Clade,    valid_tips$label)
  
  pre_samp  <- intersect(pre_samp,  rownames(patristic_mat))
  post_samp <- intersect(post_samp, rownames(patristic_mat))
  
  # Pairwise comparisons
  df_pre_pre <- extract_pairs(
    pre_samp, pre_samp, 'Pre vs Pre',
    symmetric = TRUE,
    pat_mat = patristic_mat, year_vec = year_vec, clade_vec = clade_vec
  )
  
  df_post_post <- extract_pairs(
    post_samp, post_samp, 'Post vs Post',
    symmetric = TRUE,
    pat_mat = patristic_mat, year_vec = year_vec, clade_vec = clade_vec
  )
  
  df_pre_post <- extract_pairs(
    pre_samp, post_samp, 'Pre vs Post',
    symmetric = FALSE,
    pat_mat = patristic_mat, year_vec = year_vec, clade_vec = clade_vec
  )
  
  all_pairs_df <- bind_rows(df_pre_pre, df_post_post, df_pre_post) %>%
    mutate(subtype = this)
  
  # Clade-level summary (Pre vs Post only)
  ignore_clades <- c('Unknown', 'unassigned', NA)
  
  clade_summary <- df_pre_post %>%
    filter(
      !clade_a %in% ignore_clades,
      !clade_b %in% ignore_clades
    ) %>%
    group_by(clade_a, clade_b) %>%
    summarise(
      n                  = n(),
      mean_tmrca         = mean(tmrca,     na.rm = TRUE),
      median_tmrca       = median(tmrca,   na.rm = TRUE),
      sd_tmrca           = sd(tmrca,       na.rm = TRUE),
      mean_patristic_yrs = mean(patristic, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(n >= 5) %>%
    mutate(subtype = this) %>%
    rename(pre_clade = clade_a, post_clade = clade_b)
  
  # Overall summary
  stat_summary <- all_pairs_df %>%
    group_by(comparison) %>%
    summarise(
      n                  = n(),
      mean_tmrca         = mean(tmrca,     na.rm = TRUE),
      median_tmrca       = median(tmrca,   na.rm = TRUE),
      sd_tmrca           = sd(tmrca,       na.rm = TRUE),
      mean_patristic_yrs = mean(patristic, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(subtype = this)
  
  # Test: Pre-Post TMRCA vs Post-Post TMRCA
  # H0: TMRCA(Pre,Post) >= TMRCA(Post,Post)  (local continuity)
  # H1: TMRCA(Pre,Post) <  TMRCA(Post,Post)  (importation-driven discontinuity)
  t_postp <- df_post_post$tmrca[!is.na(df_post_post$tmrca)]
  t_cpp   <- df_pre_post$tmrca[!is.na(df_pre_post$tmrca)]
  
  wt_result <- list(statistic = NA_real_, p.value = NA_real_)
  if (length(t_postp) > 1 && length(t_cpp) > 1) {
    wt_result <- wilcox.test(t_cpp, t_postp, alternative = 'less', exact = FALSE)
    message(sprintf(
      '  [%s] Wilcoxon (Pre-Post < Post-Post): W=%.0f, p=%.2e  %s',
      this, wt_result$statistic, wt_result$p.value,
      ifelse(wt_result$p.value < 0.001, '***',
             ifelse(wt_result$p.value < 0.01, '**',
                    ifelse(wt_result$p.value < 0.05, '*', 'ns')))
    ))
  }
  
  # Effect size: Post-Post mean minus Pre-Post mean
  tmrca_gap <- mean(t_postp) - mean(t_cpp)
  message(sprintf(
    '  [%s] TMRCA Gap (Post-Post mean - Pre-Post mean): %.2f years',
    this, tmrca_gap
  ))

  list(
    subtype       = this,
    all_pairs_df  = all_pairs_df,
    clade_summary = clade_summary,
    stat_summary  = stat_summary,
    tmrca_gap     = tmrca_gap,
    wilcox_W      = wt_result$statistic,
    wilcox_p      = wt_result$p.value,
    n_pre         = length(pre_labels_all),
    n_post        = length(post_labels_all)
  )
}

# Run all subtypes
results <- map(these, analyze_subtype)
names(results) <- these
results <- Filter(Negate(is.null), results)

# Combine
all_pairs <- map_dfr(results, ~ .x$all_pairs_df)
all_clades <- map_dfr(results, ~ .x$clade_summary)
all_stats <- map_dfr(results, ~ .x$stat_summary)

# Subtype labels
label_subtype <- function(df) {
  df %>% mutate(
    subtype_label = case_when(
      subtype == 'H1N1' ~ 'A/H1N1pdm09',
      subtype == 'H3N2' ~ 'A/H3N2',
      subtype == 'B'    ~ 'B/Victoria'
    ) %>% factor(levels = c('A/H1N1pdm09', 'A/H3N2', 'B/Victoria'))
  )
}

# Plot data
plot_pairs <- all_pairs %>%
  label_subtype() %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c('Pre vs Pre', 'Pre vs Post', 'Post vs Post')
    )
  ) %>%
  filter(!is.na(tmrca))

# Fig 1: TMRCA density distributions
median_lines <- plot_pairs %>%
  group_by(subtype_label, comparison) %>%
  summarise(med_tmrca = median(tmrca, na.rm = TRUE), .groups = 'drop')

fig1 <- ggplot(plot_pairs, aes(x = tmrca, fill = comparison, color = comparison)) +
  annotate('rect', xmin = NPI_FROM, xmax = NPI_TO, ymin = -Inf, ymax = Inf, alpha = 0.09, fill = 'grey30') +
  annotate('text', x = mean(c(NPI_FROM, NPI_TO)), y = Inf, label = 'COVID-19\nloackdowns', vjust = 1.5, size = base.size / 5.8, family = base.family, color = 'grey40', lineheight = 0.8) +
  geom_density(alpha = 0.28, linewidth = 0.4, adjust = 1.1) +
  geom_vline(data = median_lines, aes(xintercept = med_tmrca, color = comparison), linetype = 'dashed', linewidth = 0.35) +
  geom_text(
    data = median_lines,
    aes(x = med_tmrca, y = Inf, label = round(med_tmrca, 1), color = comparison),
    vjust = -0.3, hjust = 0.5,
    size = base.size / 6, family = base.family,
    show.legend = FALSE
  ) +
  scale_fill_manual(values  = comp_colors, name = NULL) +
  scale_color_manual(values = comp_colors, name = NULL) +
  scale_x_continuous(
    name = 'Estimated TMRCA (Year)',
    limits = c(2009, 2025.5),
    breaks = seq(2009, 2025, 2)
  ) +
  scale_y_continuous(name = 'Density', expand = expansion(mult = c(0, 0.18))) +
  facet_wrap(~ subtype_label, ncol = 1, scales = 'free_y') +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position.inside = c(0.20, 0.92),
    legend.position    = 'inside',
    legend.text        = element_text(size = base.size * 0.6, color = base.col),
    legend.background  = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = base.size * 0.6, face = 'italic', color = base.col),
    panel.spacing.y    = unit(0.1, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    axis.text          = element_text(color = base.col, size = base.size * 0.6),
    axis.title         = element_text(color = base.col, face = 'bold', size  = base.size * 0.7),
    axis.ticks         = element_line(color = '#000000', linewidth = 0.3),
    axis.title.x       = element_text(margin = margin(t = 10))
  ) + guides(
    fill  = guide_legend(override.aes = list(color = NA, fill = comp_colors, alpha = 0.3, linewidth = 0.6)),
    color = "none"
  )

ggsave(fig1, filename = glue('{fig.dir}/{ofig}_tmrca_distributions.pdf'), width = 5, height = 8, units = 'in', bg = '#FFFFFF')

# Fig 2: Clade-level lollipop plot
clade_plot <- all_clades %>%
  label_subtype() %>%
  filter(!is.na(mean_tmrca)) %>%
  mutate(
    clade_pair = glue('{pre_clade} to {post_clade}'),
    ci_lo = mean_tmrca - sd_tmrca,
    ci_hi = mean_tmrca + sd_tmrca
  )

fig2 <- ggplot(clade_plot, aes(y = reorder(clade_pair, mean_tmrca), x = mean_tmrca, color = subtype_label)) +
  annotate('rect', xmin = NPI_FROM, xmax = NPI_TO, ymin = -Inf, ymax = Inf, alpha = 0.09, fill = 'grey30') +
  annotate('text', x = mean(c(NPI_FROM, NPI_TO)), y = Inf, label = 'NPI Period', hjust = 0.5, vjust = -0.5, size = base.size / 6, family = base.family, color = 'grey40') +
  geom_segment(
    aes(xend = NPI_FROM, yend = reorder(clade_pair, mean_tmrca)),
    linewidth = 0.6, alpha = 0.45, linetype = 'dotted'
  ) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.45, linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 0.5, alpha = 0.92) +
  geom_text(aes(label = round(mean_tmrca, 1)), hjust = 1.45, vjust = 0.35, size = base.size / 5.8, family = base.family, color = base.col, show.legend = FALSE) +
  scale_x_continuous(
    name   = 'Mean TMRCA \u00b1 SD (Year)',
    limits = c(2009, 2026),
    breaks = seq(2009, 2025, 2)
  ) +
  ggsci::scale_color_npg(name = 'Subtype') +
  labs(y = 'Clade Pair\n(Pre-Pandemic to Post-Pandemic)') +
  facet_wrap(~ subtype_label, scales = 'free_y', space = 'free_y', ncol = 1) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position    = 'none',
    strip.background   = element_blank(),
    strip.text         = element_text(size = base.size * 0.6, face = 'italic', color = base.col),
    panel.spacing.y    = unit(0.1, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    axis.text          = element_text(color = base.col, size = base.size * 0.6),
    axis.title         = element_text(color = base.col, face = 'bold', size  = base.size * 0.7),
    axis.ticks         = element_line(color = '#000000', linewidth = 0.3),
    axis.title.x       = element_text(margin = margin(t = 10))
  )

ggsave(fig2, filename = glue('{fig.dir}/{ofig}_clade_tmrca.pdf'), width = 5, height = 8, units = 'in', bg = '#FFFFFF')

# Fig 3: Mean TMRCA with 95% CI
summ_plot <- all_stats %>%
  label_subtype() %>%
  mutate(
    comparison = factor(comparison, levels = c('Pre vs Pre', 'Pre vs Post', 'Post vs Post')),
    se    = sd_tmrca / sqrt(n),
    ci_lo = mean_tmrca - 1.96 * se,
    ci_hi = mean_tmrca + 1.96 * se
  )

npi_bg_df <- summ_plot %>%
  distinct(subtype_label) %>%
  mutate(ymin = NPI_FROM, ymax = NPI_TO)

fig3 <- ggplot(summ_plot, aes(x = comparison, y = mean_tmrca, color = comparison)) +
  
  geom_rect(data = npi_bg_df, aes(ymin = ymin, ymax = ymax), xmin = -Inf, xmax = Inf, inherit.aes = FALSE, alpha = 0.09, fill = 'grey30') +
  
  annotate(
    'text', x = 'Post vs Post', y = mean(c(NPI_FROM, NPI_TO)) - 1, label = 'COVID-19 loackdowns',
    hjust = 0.5, vjust = -0.8, size = base.size / 6, family = base.family, color = 'grey40'
  ) +
  
  geom_hline(
    data = summ_plot %>% filter(comparison == 'Pre vs Pre'), aes(yintercept = mean_tmrca), inherit.aes = FALSE,
    color = col.pre, linetype = 'dotted', linewidth = 0.9, alpha = 0.6
  ) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.10, linewidth = 0.5) +
  geom_point(size = 0.8) +
  geom_text(aes(label = round(mean_tmrca, 1)), vjust = -1.1, size = base.size / 5, family = base.family, color = base.col, show.legend = FALSE) +
  
  scale_color_manual(values = comp_colors, name = NULL) +
  scale_x_discrete(limits = c('Pre vs Pre', 'Pre vs Post', 'Post vs Post')) +
  scale_y_continuous(
    name   = 'Mean TMRCA (Year) \u00b1 95% CI',
    breaks = seq(2009, 2025, 2),
    limits = c(2009, 2025.5)
  ) +
  
  labs(x = NULL) +
  facet_wrap(. ~ subtype_label, nrow = 1) +
  
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position    = 'none',
    strip.background   = element_blank(),
    strip.text         = element_text(size = base.size * 0.6, face = 'italic', color = base.col),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    axis.text.x        = element_text(angle = 0, color = base.col, size = base.size * 0.6),
    axis.text.y        = element_text(color = base.col, size = base.size * 0.6),
    axis.ticks         = element_line(color = '#000000', linewidth = 0.3),
    axis.title         = element_text(color = base.col, face = 'bold', size = base.size * 0.7)
  )

ggsave(fig3, filename = glue('{fig.dir}/{ofig}_tmrca_summary.pdf'), width = 10, height = 3.5, units = 'in', bg = '#FFFFFF')

# Output tables
table_stats <- all_stats %>%
  label_subtype() %>%
  transmute(
    Subtype                      = as.character(subtype_label),
    Comparison                   = comparison,
    `N Pairs`                    = n,
    `Mean TMRCA (Year)`          = round(mean_tmrca,         2),
    `Median TMRCA (Year)`        = round(median_tmrca,       2),
    `SD`                         = round(sd_tmrca,           2),
    `Mean Patristic Dist (Yrs)`  = round(mean_patristic_yrs, 2)
  ) %>%
  arrange(Subtype, Comparison)

table_clades <- all_clades %>%
  label_subtype() %>%
  transmute(
    Subtype                      = as.character(subtype_label),
    `Pre-Pandemic Clade`         = pre_clade,
    `Post-Pandemic Clade`        = post_clade,
    `N Pairs`                    = n,
    `Mean TMRCA (Year)`          = round(mean_tmrca,         2),
    `Median TMRCA (Year)`        = round(median_tmrca,       2),
    `SD TMRCA`                   = round(sd_tmrca,           2),
    `Mean Patristic Dist (Yrs)`  = round(mean_patristic_yrs, 2)
  ) %>%
  arrange(Subtype, `Mean TMRCA (Year)`)

# Gap effect size: Post-Post mean minus Pre-Post mean
table_gap <- map_dfr(results, function(r) {
  postp <- r$stat_summary %>% filter(comparison == 'Post vs Post') %>% pull(mean_tmrca)
  cpp   <- r$stat_summary %>% filter(comparison == 'Pre vs Post')  %>% pull(mean_tmrca)
  tibble(
    Subtype                      = r$subtype,
    `Mean TMRCA: Post vs Post`   = round(postp,       2),
    `Mean TMRCA: Pre vs Post`    = round(cpp,         2),
    `TMRCA Gap (Years)`          = round(postp - cpp, 2),
    `Wilcoxon W`                 = round(r$wilcox_W,  0),
    `Wilcoxon p (one-sided)`     = signif(r$wilcox_p, 3),
    `Significance`               = ifelse(r$wilcox_p < 0.001, '***',
                                          ifelse(r$wilcox_p < 0.01, '**',
                                                 ifelse(r$wilcox_p < 0.05, '*', 'ns'))),
    `Interpretation`             = ifelse(
      (postp - cpp) > 2,
      'Pre-Post TMRCA substantially older than Post-Post: supports importation-driven lineage replacement',
      'Similar TMRCA: consistent with continuous evolution'
    )
  )
})

openxlsx::write.xlsx(
  list(
    'TMRCA Summary'     = table_stats,
    'Clade-Level TMRCA' = table_clades,
    'TMRCA Gap Effect'  = table_gap
  ),
  file = glue('{tbl.dir}/Table.Lineage_Continuity_Test.xlsx')
)