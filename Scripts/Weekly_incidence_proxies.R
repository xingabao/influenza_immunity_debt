# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lubridate)))
suppressMessages(suppressWarnings(library(cowplot)))
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
dat.fac <- readxl::read_excel(glue('{dat.dir}/Macau/FLU-AQ-CLIMATE.xlsx'))

# Prepare plot data
UDFAB <- dat.fac %>% dplyr::select(KID, FLUAB) %>% mutate(KID = ymd(KID))
UDFA <- dat.fac %>% dplyr::select(KID, FLUA) %>% mutate(KID = ymd(KID))
UDFB <- dat.fac %>% dplyr::select(KID, FLUB) %>% mutate(KID = ymd(KID))

colnames(UDFAB) <- c('KID', 'FLU')
colnames(UDFA) <- c('KID', 'FLU')
colnames(UDFB) <- c('KID', 'FLU')

prepare <- function(UDF) {
  colnames(UDF) <- c('KID', 'FLU')
  
  threshold <- quantile(UDF$FLU, 0.5)
  threshold.peak <- quantile(UDF$FLU, 0.85)
  
  UDF$KID <- as.POSIXct(UDF$KID)
  UDF <- UDF[order(UDF$KID), ]
  
  cross_index <- which(
    (head(UDF$FLU, -1) < threshold & tail(UDF$FLU, -1) >= threshold) |
      (head(UDF$FLU, -1) >= threshold & tail(UDF$FLU, -1) < threshold)
  )
  
  points <- data.frame(
    KID_prev = UDF$KID[cross_index],
    FLU_prev = UDF$FLU[cross_index],
    KID_next = UDF$KID[cross_index + 1],
    FLU_next = UDF$FLU[cross_index + 1]
  )
  
  points$KID_selected <- if_else(
    points$FLU_prev < points$FLU_next, 
    as.POSIXct(points$KID_prev), 
    as.POSIXct(points$KID_next)
  )
  
  points$FLU_selected <- if_else(
    points$FLU_prev < points$FLU_next, 
    points$FLU_prev, 
    points$FLU_next
  )
  
  KIDs <- points$KID_selected
  index_pairs <- data.frame(
    first = 1:(length(KIDs) - 1),
    second = 2:length(KIDs)
  )
  
  peaks <- data.frame()
  cindex <- 1
  for (index in 1:nrow(index_pairs)) {
    first <- KIDs[index_pairs[index, 'first']]
    second <- KIDs[index_pairs[index, 'second']]
    duration <- as.numeric(as.Date(second) - as.Date(first))
    if (duration > 30) {
      sub <- UDF[(UDF$KID >= first & (UDF$KID <= second)), ]
      if (max(sub$FLU) >= threshold.peak) {
        Middle <- sub[sub$FLU == max(sub$FLU), ]$KID
        tmp <- data.frame(Group = cindex, Start = first, Middle = Middle, End = second, Duration = duration)
        cindex <- cindex + 1
        peaks <- rbind(peaks, tmp)
      }
    }
  }
  
  peaks$Start <- as.POSIXct(peaks$Start)
  peaks$End <- as.POSIXct(peaks$End)
  
  res <- list()
  res$peaks <- peaks
  res$points <- points
  res$UDF <- UDF
  res$threshold <- threshold
  
  return(res)
}

# Prepare data
RESAB <- prepare(UDF = UDFAB)
RESA <- prepare(UDF = UDFA)
RESB <- prepare(UDF = UDFB)

# Plot function
plot. <- function(RES, ylab) {
  ggplot(RES$UDF, aes(x = KID, y = FLU)) +
    geom_line(color = "#141473") +
    geom_hline(yintercept = RES$threshold, linetype = "dashed", color = "#E50914") +
    geom_point(data = RES$points, aes(x = KID_selected, y = FLU_selected), color = '#E50914', size = 1.5, alpha = 0.5) +
    geom_rect(data = RES$peaks, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, alpha = 0.2, fill = "#00A087") +
    geom_text(
      data = RES$peaks,
      aes(x = as.POSIXct((as.numeric(Start) + as.numeric(End)) / 2), y = 0.55), label = RES$peaks$Group,
      family = base.family, color = base.col, size = base.size / 3.88
    ) +
    scale_x_datetime(
      breaks = seq(as.POSIXct("2010-01-01"), as.POSIXct("2026-01-01"), by = "1 year"),
      labels = format(seq(as.POSIXct("2010-01-01"), as.POSIXct("2026-01-01"), by = "1 year"), "%Y"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, .6)) +
    labs(y = ylab, x = NULL) +
    coord_cartesian(xlim = c(as.POSIXct("2010-01-01"), as.POSIXct("2026-04-01"))) +
    theme_classic() +
    theme(
      axis.text.x = element_text(margin = margin(t = 5), color = base.col, size = base.size * 0.9, angle = 30, hjust = 1, family = base.family),
      axis.text.y = element_text(margin = margin(r = 5), color = base.col, size = base.size * 0.9, family = base.family),
      axis.title.x = element_text(margin = margin(t = 15), color = base.col, size = base.size * 0.9, family = base.family),
      axis.title.y = element_text(margin = margin(r = 15, l = 10), color = base.col, size = base.size * 0.9, family = base.family)
    )
}

ggA <- plot.(RES = RESA, ylab = 'Incidence of Flu A') + theme(axis.title.x = element_blank(), axis.text.x  = element_blank())
ggB <- plot.(RES = RESB, ylab = 'Incidence of Flu B') + theme(axis.title.x = element_blank(), axis.text.x  = element_blank())
ggAB <- plot.(RES = RESAB, ylab = 'Incidence of Flu A + B')

openxlsx::write.xlsx(x = RESA$peaks, file = glue('{tbl.dir}/Table.Influenza.Peaks.A.xlsx'))
openxlsx::write.xlsx(x = RESB$peaks, file = glue('{tbl.dir}/Table.Influenza.Peaks.B.xlsx'))
openxlsx::write.xlsx(x = RESAB$peaks, file = glue('{tbl.dir}/Table.Influenza.Peaks.AB.xlsx'))

# Combination
main_gg <- (ggA / ggB / ggAB) + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = list(c('', '', ''))) & theme(plot.tag = element_text(family = base.family, color = base.col, size = base.size))

color_bar_grob <- ggplotGrob(
  ggplot() +
    annotate("rect", xmin = as.POSIXct("2011-11-01"), xmax = as.POSIXct("2020-03-31"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#99999999") +
    annotate("rect", xmin = as.POSIXct("2020-04-01"), xmax = as.POSIXct("2021-03-31"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#EE000055") +
    annotate("rect", xmin = as.POSIXct("2021-04-01"), xmax = as.POSIXct("2023-04-30"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#3B499255") +
    annotate("rect", xmin = as.POSIXct("2023-05-01"), xmax = as.POSIXct("2025-10-30"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#99999999") +
    scale_x_datetime(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    coord_cartesian(xlim = c(as.POSIXct("2010-01-01"), as.POSIXct("2025-12-31"))) +
    theme_void(base_size = base.size * 0.9, base_family = base.family) +
    theme(plot.margin = margin(0, 0, 0, 0))
)

font.size <- base.size / 3.88
lineheight <- 0.9
lwd.size <- 0.4

text_desc_grob <- ggplotGrob(
  ggplot() +
    annotate("text", x = as.POSIXct("2016-08-15"), y = 100, size = font.size * 0.9, label = "Pre-pandemic period", lineheight = lineheight, family = base.family, vjust = 0.4) +
    geom_segment(aes(x = as.POSIXct("2020-10-01"), xend = as.POSIXct("2020-04-01"), y = 70, yend = 70), arrow = arrow(length = unit(0.1, 'cm')), linewidth = lwd.size) +
    geom_segment(aes(x = as.POSIXct("2022-12-01"), xend = as.POSIXct("2023-04-15"), y = 70, yend = 70), arrow = arrow(length = unit(0.1, 'cm')), linewidth = lwd.size) +
    annotate("text", x = as.POSIXct("2021-11-01"), y = 100, size = font.size * 0.9, label = "Pandemic period", lineheight = lineheight, family = base.family, vjust = 0.4) +
    annotate("text", x = as.POSIXct("2020-10-01"), y = 25, size = font.size * 0.7, label = "(Acute phase)", color = "#EE0000") +
    annotate("text", x = as.POSIXct("2022-04-15"), y = 25, size = font.size * 0.7, label = "(Transition phase)", color = "#3B4992") +
    annotate("text", x = as.POSIXct("2024-08-30"), y = 100, size = font.size * 0.9, label = "Post-pandemic period", lineheight = lineheight, family = base.family, vjust = 0.4) +
    scale_x_datetime(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    coord_cartesian(xlim = c(as.POSIXct("2010-01-01"), as.POSIXct("2025-12-31")), clip = 'off') +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
)

jst <- 0.008
gg <- ggdraw(main_gg) +  
  draw_grob(text_desc_grob, x = 0, y = 0.500 + jst, width = 1, height = 0.150) + 
  draw_grob(color_bar_grob, x = 0, y = 0.947 + jst, width = 1, height = 0.025) + 
  draw_grob(color_bar_grob, x = 0, y = 0.639 + jst, width = 1, height = 0.025) + 
  draw_grob(color_bar_grob, x = 0, y = 0.332 + jst, width = 1, height = 0.025)

# Save to file
width = 10; height = 8
ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, bg = '#FFFFFF')