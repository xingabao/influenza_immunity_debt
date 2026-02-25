# Load R packages
suppressMessages(suppressWarnings(library(mem)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lubridate)))
Sys.setlocale("LC_TIME", "C")

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
this <- 'HK'
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Define functions
memsurveillance. <- function(
    i.current, i.epidemic.thresholds = NA, i.intensity.thresholds = NA, 
    i.mean.length = 10, i.force.length = FALSE, i.output = ".", 
    i.graph.title = "", i.graph.subtitle = "", i.graph.file = FALSE, 
    i.graph.file.name = "", i.week.report = NA, i.equal = FALSE, 
    i.pos.epidemic = FALSE, i.no.epidemic = FALSE, i.no.intensity = FALSE, 
    i.epidemic.start = NA, i.range.x = c(40, 20), i.range.x.53 = FALSE, 
    i.range.y = NA, i.no.labels = FALSE, i.start.end.marks = TRUE, 
    i.mem.info = TRUE,
    legend.x = NULL, legend.y = NULL, legend.cex = 1, legend.y.intersp = 0.8, legend.ncol = 1,
    axis.text.cex = 1, axis.text.line = 0.7, axis.text.t.col = '#000000', axis.text.l.col = '#000000',
    axis.title.x = "Week", axis.title.y = "Weekly rate",
    axis.title.cex = 1, axis.title.col = '#000000', axis.title.line.x = 2.5, axis.title.line.y = 1.8
) {
  if (is.null(dim(i.current))) {
    stop("Incorrect number of dimensions, input must be a data.frame.")
  } else if (!(ncol(i.current) == 1)) {
    stop("Incorrect number of dimensions, only one season required.")
  }
  if (!is.numeric(i.epidemic.thresholds) || length(i.epidemic.thresholds) == 1) i.epidemic.thresholds <- rep(NA, 2)
  if (!is.numeric(i.intensity.thresholds) || length(i.intensity.thresholds) == 1) i.intensity.thresholds <- rep(NA, 3)
  if (!is.numeric(i.range.x) || length(i.range.x) != 2) i.range.x <- c(40, 20)
  if (i.range.x.53) { esquema.temporadas.1 <- 53 } else { esquema.temporadas.1 <- 52 }
  if (i.range.x[1] == i.range.x[2]) i.range.x[2] <- i.range.x[1] - 1
  if (i.range.x[1] < i.range.x[2]) {
    esquema.temporadas.2 <- max(1, i.range.x[1])
    esquema.temporadas.3 <- min(esquema.temporadas.1, i.range.x[2])
    esquema.temporadas.4 <- c(esquema.temporadas.2:esquema.temporadas.3)
  } else {
    esquema.temporadas.2 <- min(esquema.temporadas.1, i.range.x[1])
    esquema.temporadas.3 <- max(1, i.range.x[2])
    esquema.temporadas.4 <- c(esquema.temporadas.2:esquema.temporadas.1, 1:esquema.temporadas.3)
  }
  semanas <- length(esquema.temporadas.4)
  esquema.semanas <- data.frame(numero.semana = 1:semanas, nombre.semana = esquema.temporadas.4)
  current.season <- i.current
  names(current.season) <- "rates"
  current.season$nombre.semana <- rownames(i.current)
  rownames(current.season) <- NULL
  current.season <- merge(esquema.semanas, current.season, by = "nombre.semana", all.x = TRUE)
  current.season <- current.season[order(current.season$numero.semana), ]
  rownames(current.season) <- NULL
  if (!is.na(i.week.report) && any(i.week.report == as.numeric(esquema.semanas$nombre.semana))) {
    semana.report <- ((1:semanas)[i.week.report == as.numeric(esquema.semanas$nombre.semana)])[1]
    if (!is.na(semana.report) && semana.report < semanas) 
      current.season$rates[(semana.report + 1):semanas] <- NA
  } else {
    if (all(is.na(current.season$rates))) { semana.report <- semanas } else { semana.report <- max((1:semanas)[!is.na(current.season$rates)], na.rm = TRUE) }
  }
  umbral.pre <- as.numeric(i.epidemic.thresholds[1])
  if (i.equal) { umbral.pos <- as.numeric(i.epidemic.thresholds[1]) } else { umbral.pos <- as.numeric(i.epidemic.thresholds[2]) }
  duracion.media <- i.mean.length
  if (!is.na(i.epidemic.start)) { semana.inicio.forzado <- ((1:semanas)[i.epidemic.start == as.numeric(esquema.semanas$nombre.semana)])[1] } else { semana.inicio.forzado <- NA }
  if (any(current.season$rates > umbral.pre, na.rm = TRUE)) { semana.inicio.real <- min((1:semanas)[current.season$rates > umbral.pre], na.rm = TRUE) } else { semana.inicio.real <- NA }
  if (!is.na(semana.inicio.forzado)) { if (semana.inicio.forzado > semana.report) semana.inicio.forzado <- NA }
  if (!is.na(semana.inicio.forzado) && !is.na(semana.inicio.real)) { if (semana.inicio.forzado == semana.inicio.real) semana.inicio.forzado <- NA }
  if (!is.na(semana.inicio.forzado)) { semana.inicio <- semana.inicio.forzado } else { semana.inicio <- semana.inicio.real }
  week.peak <- which.max(current.season$rates)
  if (!is.na(semana.inicio)) {
    if (i.force.length) {
      semana.fin <- semana.inicio + i.mean.length
      if (semana.fin > semanas) semana.fin <- NA
    } else {
      punto.de.busqueda <- max(semana.inicio, semana.inicio.real, week.peak, na.rm = TRUE)
      semana.fin.1 <- (1:semanas)[current.season$rates < umbral.pos & punto.de.busqueda < (1:semanas)]
      if (any(semana.fin.1, na.rm = TRUE)) { semana.fin <- min(semana.fin.1, na.rm = TRUE) } else { semana.fin <- NA }
    }
  } else {
    semana.fin <- NA
  }
  if (i.no.epidemic) {
    semana.inicio <- NA
    semana.fin <- NA
  }
  limites.niveles <- as.vector(i.intensity.thresholds)
  limites.niveles[limites.niveles < 0] <- 0
  if (is.na(semana.inicio)) {
    umbrales.1 <- rep(umbral.pre, semana.report + 1)
    umbrales.2 <- rep(NA, semanas)
    intensidades.1 <- array(dim = c(semanas, 3))
    intensidades.2 <- array(dim = c(semanas, 3))
  } else {
    if (is.na(semana.fin)) {
      umbrales.1 <- rep(umbral.pre, semana.inicio - 1)
      umbrales.2 <- rep(NA, max(duracion.media, semana.report - semana.inicio + 1))
      if (!i.no.intensity) {
        intensidades.1 <- array(dim = c(semana.inicio - 1, 3))
        intensidades.2 <- matrix(rep(limites.niveles, max(duracion.media, semana.report - semana.inicio + 1)), ncol = 3, byrow = TRUE)
      } else {
        intensidades.1 <- array(dim = c(semana.inicio - 1, 3))
        intensidades.2 <- array(dim = c(max(duracion.media, semana.report - semana.inicio + 1), 3))
      }
    } else {
      umbrales.1 <- rep(umbral.pre, semana.inicio - 1)
      umbrales.2 <- rep(NA, semana.fin - semana.inicio)
      if (!i.no.intensity) {
        intensidades.1 <- array(dim = c(semana.inicio - 1, 3))
        intensidades.2 <- matrix(rep(limites.niveles, semana.fin - semana.inicio), ncol = 3, byrow = TRUE)
      } else {
        intensidades.1 <- array(dim = c(semana.inicio - 1, 3))
        intensidades.2 <- array(dim = c(semana.fin - semana.inicio, 3))
      }
    }
  }
  if (i.pos.epidemic) { umbrales.3 <- rep(umbral.pos, semanas) } else { umbrales.3 <- rep(NA, semanas) }
  umbrales <- c(umbrales.1, umbrales.2, umbrales.3)[1:semanas]
  intensidades.3 <- array(dim = c(semanas, 3))
  intensidades <- rbind(intensidades.1, intensidades.2, intensidades.3)[1:semanas, ]
  dgraf <- as.data.frame(cbind(current.season$rates, umbrales, intensidades))
  names(dgraf) <- c("Value", "Epidemic threshold", paste("Intensidad", 1:3))
  if (i.graph.file.name == "") { graph.name <- "mem surveillance graph" } else { graph.name <- i.graph.file.name }
  etiquetas <- c("Weekly values", "Epidemic", "Medium", "High", "Very high")
  tipos <- c(1, 2, 2, 2, 2)
  anchos <- c(3, 2, 2, 2, 2)
  colores <- c("#808080", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")
  if (is.numeric(i.range.y)) { range.y.bus <- i.range.y } else { range.y.bus <- c(0, mem:::maxFixNA(dgraf)) }
  if (FALSE) {
    otick <- mem:::optimal.tickmarks(range.y.bus[1], range.y.bus[2], 10)
  } else {
    otick <- list()
    otick[['by']] <- 0.05
    otick[['number']] <- 12
    otick[['range']] <- c(0.00, 0.55)
    otick[['tickmarks']] <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55)
  }
  range.y <- c(otick$range[1], otick$range[2] + otick$by/2)
  if (i.graph.file) {
    tiff(filename = paste(i.output, "/", graph.name, ".tiff", sep = ""), width = 8, height = 6, units = "in", pointsize = "12", compression = "lzw", bg = "white", res = 300, antialias = "none")
  }
  matplot(
    1:semanas, dgraf, type = "l", lty = tipos, lwd = anchos, 
    col = colores, xlab = "", ylab = "", axes = FALSE, ylim = range.y, 
    main = i.graph.title
  )
  points(1:semanas, dgraf[, 1], pch = 19, type = "p", col = "#000000", cex = 1)
  if (i.start.end.marks) {
    if (!is.na(semana.inicio)) 
      points(x = semana.inicio, y = current.season$rates[semana.inicio], pch = 1, bg = "#FFFFFF", col = "#FF0000", lwd = 7)
    if (!is.na(semana.fin) && i.pos.epidemic) {
      if (is.na(current.season$rates[semana.fin])) {
        points(x = semana.fin, y = 0, pch = 13, bg = "#FFFFFF", col = "#40FF40", lwd = 7)
      } else {
        points(x = semana.fin, y = current.season$rates[semana.fin], pch = 1, bg = "#FFFFFF", col = "#40FF40", lwd = 7)
      }
    }
  }
  axis(1, at = seq(1, semanas, 1), labels = FALSE, cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
  axis(1, at = seq(1, semanas, 2), tick = FALSE, labels = esquema.semanas$nombre.semana[seq(1, semanas, 2)], cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
  axis(1, at = seq(2, semanas, 2), tick = FALSE, labels = esquema.semanas$nombre.semana[seq(2, semanas, 2)], cex.axis = axis.text.cex, line = axis.text.line, col.axis = axis.text.t.col, col = axis.text.l.col)
  axis(2, at = otick$tickmarks, lwd = 1, cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
  mtext(1, text = axis.title.x, line = axis.title.line.x, cex = axis.title.cex, col = axis.title.col)
  mtext(2, text = axis.title.y, line = axis.title.line.y, cex = axis.title.cex, col = axis.title.col)
  mtext(3, text = i.graph.subtitle, cex = axis.title.cex, col = axis.title.col)
  if (!i.no.labels) {
    if (is.na(semana.inicio)) {
      text.x <- semana.report - semanas/25
      text.y <- umbral.pre
      text.l <- round(umbral.pre, 2)
      text.p <- 3
      text.s <- 1
      text.c <- colores[2]
    } else {
      if (is.na(semana.fin)) { lugar.intensidad <- min(semanas, semana.inicio + max(duracion.media, semana.report - semana.inicio + 1)) } else { lugar.intensidad <- semana.fin - 1 }
      text.x <- c(semana.inicio, rep(lugar.intensidad, 3), semanas) - semanas/25
      text.y <- c(umbral.pre, limites.niveles, umbral.pos)
      text.l <- round(c(umbral.pre, limites.niveles, umbral.pos), 2)
      text.p <- rep(3, 5)
      text.s <- rep(1, 5)
      text.c <- colores[c(2:5, 2)]
      quitar.columnas <- numeric()
      if (!i.pos.epidemic || is.na(semana.fin)) quitar.columnas <- c(quitar.columnas, 5)
      if (i.no.intensity) quitar.columnas <- c(quitar.columnas, 2:4)
      if (length(quitar.columnas) > 0) {
        text.x <- text.x[-quitar.columnas]
        text.y <- text.y[-quitar.columnas]
        text.l <- text.l[-quitar.columnas]
        text.p <- text.p[-quitar.columnas]
        text.s <- text.s[-quitar.columnas]
        text.c <- text.c[-quitar.columnas]
      }
    }
    text(text.x, text.y, text.l, pos = text.p, col = text.c, cex = text.s)
  }
  etiquetas.leyenda <- c("End", "Start", etiquetas)
  tipos.leyenda <- c(NA, NA, tipos)
  anchos.leyenda <- c(7, 7, anchos)
  colores.leyenda <- c("#40FF40", "#FF0000", colores)
  puntos.leyenda <- c(1, 1, rep(NA, 5))
  bg.leyenda <- c("#FFFFFF", "#FFFFFF", rep(NA, 5))
  quitar.columnas <- numeric()
  if (!i.start.end.marks || !is.na(semana.inicio.forzado)) quitar.columnas <- c(quitar.columnas, 1:2)
  if (!i.pos.epidemic || is.na(semana.fin) || is.na(i.epidemic.thresholds[2])) quitar.columnas <- c(quitar.columnas, 1)
  if (is.na(semana.inicio)) quitar.columnas <- c(quitar.columnas, 2)
  if (i.no.epidemic || is.na(i.epidemic.thresholds[1])) quitar.columnas <- c(quitar.columnas, 4)
  if (i.no.intensity) quitar.columnas <- c(quitar.columnas, 5:7)
  quitar.columnas <- c(quitar.columnas, (5:7)[is.na(i.intensity.thresholds)])
  if (length(quitar.columnas) > 0) {
    etiquetas.leyenda <- etiquetas.leyenda[-quitar.columnas]
    tipos.leyenda <- tipos.leyenda[-quitar.columnas]
    anchos.leyenda <- anchos.leyenda[-quitar.columnas]
    colores.leyenda <- colores.leyenda[-quitar.columnas]
    puntos.leyenda <- puntos.leyenda[-quitar.columnas]
    bg.leyenda <- bg.leyenda[-quitar.columnas]
  }
  if (is.numeric(legend.x) & is.numeric(legend.y)) {
    xa = legend.x; ya = legend.y
  } else {
    if (is.na(semana.inicio) || is.na(semana.fin)) {
      xa <- "topright"
      ya <- NULL
    } else {
      if (semana.fin < 0.8 * semanas) { xa <- "topright" } else { xa <- "topleft" }
      ya <- NULL
    }
  }
  legend(
    x = xa, y = ya, inset = c(0, -0.05), xjust = 0, legend = rev(etiquetas.leyenda), 
    bty = "n", lty = rev(tipos.leyenda), lwd = rev(anchos.leyenda), 
    col = rev(colores.leyenda), pch = rev(puntos.leyenda), 
    bg = rev(bg.leyenda), cex = legend.cex, x.intersp = 0.5, y.intersp = legend.y.intersp, 
    text.col = "#000000", ncol = legend.ncol
  )
  if (i.graph.file) dev.off()
}

plot.mem. <- function(
    x, legend.x = NULL, legend.y = NULL, legend.cex = 1, legend.y.intersp = 0.8, legend.ncol = 1,
    axis.text.cex = 1, axis.text.line = 0.7, axis.text.t.col = '#000000', axis.text.l.col = '#000000',
    axis.title.x = "Week", axis.title.y = "Weekly rate",
    axis.title.cex = 1, axis.title.col = '#000000', axis.title.line.x = 2.5, axis.title.line.y = 1.8,
    plot = TRUE
) {
  semanas <- dim(x$data)[1]
  anios <- dim(x$data)[2]
  names.seasons <- sub("^.*\\(([^\\)]*)\\)$", "\\1", names(x$data), perl = TRUE)
  datos.graf <- x$moving.epidemics
  colnames(datos.graf) <- names(x$data)
  lab.graf <- (1:semanas) + x$ci.start[2, 2] - x$ci.start[1, 2]
  lab.graf[lab.graf > 52] <- (lab.graf - 52)[lab.graf > 52]
  lab.graf[lab.graf < 1] <- (lab.graf + 52)[lab.graf < 1]
  tipos <- rep(1, anios)
  anchos <- rep(2, anios)
  pal <- colorRampPalette(RColorBrewer::brewer.pal(4, "Spectral"))
  colores <- pal(anios)
  if (FALSE) {
    otick <- mem:::optimal.tickmarks(0, mem:::maxFixNA(datos.graf), 10)
  } else {
    otick <- list()
    otick[['by']] <- 0.05
    otick[['number']] <- 12
    otick[['range']] <- c(0.00, 0.55)
    otick[['tickmarks']] <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55)
  }
  range.y <- c(otick$range[1], otick$range[2] + otick$by/2)
  
  if (plot) {
    matplot(
      x = 1:semanas, 
      datos.graf,
      type = "l", 
      sub = NULL, 
      lty = tipos, 
      lwd = anchos,
      col = colores,
      xlim = c(1, semanas),
      xlab = "", 
      ylab = "",
      axes = FALSE,
      font.main = 2, 
      font.sub = 1,
      ylim = range.y
    )
    axis(1, at = seq(1, semanas, 1), labels = FALSE, cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
    axis(1, at = seq(1, semanas, 2), tick = FALSE, labels = as.character(lab.graf)[seq(1, semanas, 2)], cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
    axis(1, at = seq(2, semanas, 2), tick = FALSE, labels = as.character(lab.graf)[seq(2, semanas, 2)], cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col, line = axis.text.line)
    mtext(1, text = axis.title.x, line = axis.title.line.x, cex = axis.title.cex, col = axis.title.col)
    axis(2, at = otick$tickmarks, lwd = 1, cex.axis = axis.text.cex, col.axis = axis.text.t.col, col = axis.text.l.col)
    mtext(2, text = axis.title.y, line = axis.title.line.y, cex = axis.title.cex, col = axis.title.col)
    i.temporada <- x$ci.start[1, 2]
    f.temporada <- x$ci.start[1, 2] + x$mean.length - 1
    lines(x = rep(i.temporada - 0.5, 2), y = range.y, col = "#FF0000", lty = 2, lwd = 2)
    lines(x = rep(f.temporada + 0.5, 2), y = range.y, col = "#40FF40", lty = 2, lwd = 2)
    lines(x = c(1, semanas), y = rep(x$epidemic.thresholds[1], 2), col = "#8c6bb1", lty = 2, lwd = 2)
    if (is.numeric(legend.x) & is.numeric(legend.y)) {
      xa = legend.x; ya = legend.y
    } else {
      ya <- otick$range[2]
      if ((i.temporada - 1) <= (semanas - f.temporada)) { xa <- f.temporada + 1 } else { xa <- 1 }
    }
    legend(
      x = xa, y = ya, inset = c(0, 0), xjust = 0, seg.len = 1, 
      legend = names.seasons, bty = "n", lty = tipos, lwd = anchos, 
      col = colores, cex = legend.cex, x.intersp = 0.5, y.intersp = legend.y.intersp, 
      text.col = "#000000", ncol = legend.ncol
    )
  }
  
  return(list(x = x, semanas = semanas))
}

# Prepare plot data
if (TRUE) {
  UDF.A <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUA) %>%
    filter(Year > 2009, Year < 2026)
  
  # Arrange data
  UDF.A$KID <- as.Date(UDF.A$KID)
  UDF.A$Year <- as.integer(UDF.A$Year)
  UDF.A$Week <- as.integer(UDF.A$Week)
  
  #
  UDF.mat.A <- reshape2::dcast(UDF.A, Week ~ Year, value.var = "FLUA")
  UDF.mat.A <- UDF.mat.A[order(UDF.mat.A$Week), ]
  rownames(UDF.mat.A) <- UDF.mat.A$Week
  UDF.mat.A$Week <- NULL
  UDF.mat.A <- as.data.frame(UDF.mat.A)
  
  # Create a model using the `memmodel` function from the `mem` package
  mem.model.A <- memmodel(
    i.data = UDF.mat.A,
    i.seasons = 15,
    i.type.threshold = 4,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the MEM model
  summary(mem.model.A)
  mem.model.A$centered.start
  mem.model.A$centered.length
  mem.model.A$mean.start
  mem.model.A$mean.length
  
  thresholds.A <- memintensity(mem.model.A)
  print(thresholds.A)
  
  # Pre-Pandemic
  # ---------------------------------
  UDF.A1 <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUA) %>%
    filter(Year > 2015 & Year < 2020)
  
  UDF.A1$KID <- as.Date(UDF.A1$KID)
  UDF.A1$Year <- as.integer(UDF.A1$Year)
  UDF.A1$Week <- as.integer(UDF.A1$Week)
  
  UDF.mat.A1 <- reshape2::dcast(UDF.A1, Week ~ Year, value.var = "FLUA")
  UDF.mat.A1 <- UDF.mat.A1[order(UDF.mat.A1$Week), ]
  rownames(UDF.mat.A1) <- UDF.mat.A1$Week
  UDF.mat.A1$Week <- NULL
  UDF.mat.A1 <- as.data.frame(UDF.mat.A1)
  
  mem.model.A1 <- memmodel(
    i.data = UDF.mat.A1,
    i.seasons = 10,
    i.type.threshold = 2,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the model
  summary(mem.model.A1)
  mem.model.A1$centered.start
  mem.model.A1$centered.length
  mem.model.A1$mean.start
  mem.model.A1$mean.length
  
  thresholds.A1 <- memintensity(mem.model.A1)
  print(thresholds.A1)
  
  # Post-Pandemic
  # ---------------------------------
  UDF.A2 <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUA) %>%
    filter(Year >= 2023)
  
  UDF.A2$KID <- as.Date(UDF.A2$KID)
  UDF.A2$Year <- as.integer(UDF.A2$Year)
  UDF.A2$Week <- as.integer(UDF.A2$Week)
  
  UDF.mat.A2 <- reshape2::dcast(UDF.A2, Week ~ Year, value.var = "FLUA")
  UDF.mat.A2 <- UDF.mat.A2[order(UDF.mat.A2$Week), ]
  rownames(UDF.mat.A2) <- UDF.mat.A2$Week
  UDF.mat.A2$Week <- NULL
  UDF.mat.A2 <- as.data.frame(UDF.mat.A2)
  
  mem.model.A2 <- memmodel(
    i.data = UDF.mat.A2,
    i.seasons = 10,
    i.type.threshold = 2,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the model
  summary(mem.model.A2)
  mem.model.A2$centered.start
  mem.model.A2$centered.length
  mem.model.A2$mean.start
  mem.model.A2$mean.length
  
  thresholds.A2 <- memintensity(mem.model.A2)
  print(thresholds.A2)
}

if (TRUE) {
  UDF.B <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUB) %>%
    filter(Year > 2009, Year < 2026)
  
  # Arrange data
  UDF.B$KID <- as.Date(UDF.B$KID)
  UDF.B$Year <- as.integer(UDF.B$Year)
  UDF.B$Week <- as.integer(UDF.B$Week)
  
  #
  UDF.mat.B <- reshape2::dcast(UDF.B, Week ~ Year, value.var = "FLUB")
  UDF.mat.B <- UDF.mat.B[order(UDF.mat.B$Week), ]
  rownames(UDF.mat.B) <- UDF.mat.B$Week
  UDF.mat.B$Week <- NULL
  UDF.mat.B <- as.data.frame(UDF.mat.B)
  
  # Create a model using the `memmodel` function from the `mem` package
  mem.model.B <- memmodel(
    i.data = UDF.mat.B,
    i.seasons = 15,
    i.type.threshold = 4,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the MEM model
  summary(mem.model.B)
  mem.model.B$centered.start
  mem.model.B$centered.length
  mem.model.B$mean.start
  mem.model.B$mean.length
  
  thresholds.B <- memintensity(mem.model.B)
  print(thresholds.B)
  
  # Pre-Pandemic
  # ---------------------------------
  UDF.B1 <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUB) %>%
    filter(Year > 2015 & Year < 2020)
  
  UDF.B1$KID <- as.Date(UDF.B1$KID)
  UDF.B1$Year <- as.integer(UDF.B1$Year)
  UDF.B1$Week <- as.integer(UDF.B1$Week)
  
  UDF.mat.B1 <- reshape2::dcast(UDF.B1, Week ~ Year, value.var = "FLUB")
  UDF.mat.B1 <- UDF.mat.B1[order(UDF.mat.B1$Week), ]
  rownames(UDF.mat.B1) <- UDF.mat.B1$Week
  UDF.mat.B1$Week <- NULL
  UDF.mat.B1 <- as.data.frame(UDF.mat.B1)
  
  mem.model.B1 <- memmodel(
    i.data = UDF.mat.B1,
    i.seasons = 10,
    i.type.threshold = 4,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the model
  summary(mem.model.B1)
  mem.model.B1$centered.start
  mem.model.B1$centered.length
  mem.model.B1$mean.start
  mem.model.B1$mean.length
  
  thresholds.B1 <- memintensity(mem.model.B1)
  print(thresholds.B1)
  
  # Post-Pandemic
  # ---------------------------------
  UDF.B2 <- readxl::read_excel(glue('{dat.dir}/{this}/FLU-CL-AQ.xlsx')) %>%
    mutate(Year = isoyear(date)) %>%
    mutate(Week = isoweek(date)) %>%
    mutate(KID = ymd(date)) %>%
    dplyr::select(Year, Week, KID, FLUB) %>%
    filter(Year >= 2023)
  
  UDF.B2$KID <- as.Date(UDF.B2$KID)
  UDF.B2$Year <- as.integer(UDF.B2$Year)
  UDF.B2$Week <- as.integer(UDF.B2$Week)
  
  UDF.mat.B2 <- reshape2::dcast(UDF.B2, Week ~ Year, value.var = "FLUB")
  UDF.mat.B2 <- UDF.mat.B2[order(UDF.mat.B2$Week), ]
  rownames(UDF.mat.B2) <- UDF.mat.B2$Week
  UDF.mat.B2$Week <- NULL
  UDF.mat.B2 <- as.data.frame(UDF.mat.B2)
  
  mem.model.B2 <- memmodel(
    i.data = UDF.mat.B2,
    i.seasons = 10,
    i.type.threshold = 4,
    i.mem.info = FALSE
  )
  
  # Determine epidemic thresholds based on the model
  summary(mem.model.B2)
  mem.model.B2$centered.start
  mem.model.B2$centered.length
  mem.model.B2$mean.start
  mem.model.B2$mean.length
  
  thresholds.B2 <- memintensity(mem.model.B2)
  print(thresholds.B2)
}

# Save to file
width <- 9; height <- 6
mfrow <- c(2, 3)
mar.a <- c(2, 3, 0, 1.5) + 0.1
mar.b <- c(2, 1.5, 0, 0) + 0.1
mar.a2 <- c(4, 3, 1, 1.5) + 0.1
mar.b2 <- c(4, 1.5, 1, 0) + 0.1
x.just <- 0.06

pdf(glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, family = "Times", bg = '#FFFFFF')
opar <- par(mfrow = mfrow, mar = mar.a, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
men.list <- plot.mem.(
  x = mem.model.A,
  legend.x = 28.72,
  legend.y = 0.59,
  legend.cex = 1,
  legend.y.intersp = 1.0, 
  legend.ncol = 2,
  axis.text.line = 0.8
)
usr <- par("usr")
x.pos <- usr[1] - x.just*(usr[2] - usr[1])
y.pos <- usr[4] - 0.00*(usr[4] - usr[3])

if (TRUE) {
  opar <- par(mar = mar.b, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
  men.list <- plot.mem.(
    x = mem.model.A1,
    legend.x = 30.72,
    legend.y = 0.59,
    legend.cex = 1,
    legend.y.intersp = 1.0, 
    legend.ncol = 2,
    axis.text.line = 0.8, 
    plot = FALSE
  )
  x = men.list[['x']]
  semanas = men.list[['semanas']]
  memsurveillance.(
    i.current = data.frame(rates = x$typ.curve[, 2], row.names = rownames(x$data)), 
    i.epidemic.thresholds = x$epidemic.thresholds, 
    i.intensity.thresholds = x$intensity.thresholds, 
    i.mean.length = x$mean.length, 
    i.force.length = TRUE, 
    i.graph.file = FALSE, 
    i.week.report = NA, 
    i.equal = FALSE, 
    i.pos.epidemic = TRUE,
    i.no.epidemic = FALSE, 
    i.no.intensity = FALSE,
    i.epidemic.start = x$ci.start[2, 2], 
    i.range.x = as.numeric(c(rownames(x$data)[1], rownames(x$data)[semanas])), 
    i.range.x.53 = FALSE, 
    i.range.y = NA, 
    i.no.labels = FALSE, 
    i.start.end.marks = TRUE,
    legend.x = 29.72,
    legend.y = 0.58,
    legend.cex = 1,
    legend.y.intersp = 1.0,
    axis.text.line = 0.8
  )
  usr <- par("usr")
  x.pos <- usr[1] - x.just*(usr[2] - usr[1])
  y.pos <- usr[4] - 0.00*(usr[4] - usr[3])
}

if (TRUE) {
  opar <- par(mar = mar.b, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
  men.list <- plot.mem.(
    x = mem.model.A2,
    legend.x = 30.72,
    legend.y = 0.59,
    legend.cex = 1,
    legend.y.intersp = 1.0, 
    legend.ncol = 2,
    axis.text.line = 0.8, 
    plot = FALSE
  )
  x = men.list[['x']]
  semanas = men.list[['semanas']]
  memsurveillance.(
    i.current = data.frame(rates = x$typ.curve[, 2], row.names = rownames(x$data)), 
    i.epidemic.thresholds = x$epidemic.thresholds, 
    i.intensity.thresholds = x$intensity.thresholds, 
    i.mean.length = x$mean.length, 
    i.force.length = TRUE, 
    i.graph.file = FALSE, 
    i.week.report = NA, 
    i.equal = FALSE, 
    i.pos.epidemic = TRUE,
    i.no.epidemic = FALSE, 
    i.no.intensity = FALSE,
    i.epidemic.start = x$ci.start[2, 2], 
    i.range.x = as.numeric(c(rownames(x$data)[1], rownames(x$data)[semanas])), 
    i.range.x.53 = FALSE, 
    i.range.y = NA, 
    i.no.labels = FALSE, 
    i.start.end.marks = TRUE,
    legend.x = 29.72,
    legend.y = 0.58,
    legend.cex = 1,
    legend.y.intersp = 1.0,
    axis.text.line = 0.8
  )
  usr <- par("usr")
  x.pos <- usr[1] - x.just*(usr[2] - usr[1])
  y.pos <- usr[4] - 0.00*(usr[4] - usr[3])
}

opar <- par(mar = mar.a2, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
men.list <- plot.mem.(
  x = mem.model.B,
  legend.x = 25.72,
  legend.y = 0.50,
  legend.cex = 1,
  legend.y.intersp = 1.0, 
  legend.ncol = 2,
  axis.text.line = 0.8
)
usr <- par("usr")
x.pos <- usr[1] - x.just*(usr[2] - usr[1])
y.pos <- usr[4] - 0.00*(usr[4] - usr[3])

if (TRUE) {
  opar <- par(mar = mar.b2, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
  men.list <- plot.mem.(
    x = mem.model.B1,
    legend.x = 25.72,
    legend.y = 0.50,
    legend.cex = 1,
    legend.y.intersp = 1.0, 
    legend.ncol = 2,
    axis.text.line = 0.8,
    plot = FALSE
  )
  x = men.list[['x']]
  semanas = men.list[['semanas']]
  memsurveillance.(
    i.current = data.frame(rates = x$typ.curve[, 2], row.names = rownames(x$data)), 
    i.epidemic.thresholds = x$epidemic.thresholds, 
    i.intensity.thresholds = x$intensity.thresholds, 
    i.mean.length = x$mean.length, 
    i.force.length = TRUE, 
    i.graph.file = FALSE, 
    i.week.report = NA, 
    i.equal = FALSE, 
    i.pos.epidemic = TRUE,
    i.no.epidemic = FALSE, 
    i.no.intensity = FALSE,
    i.epidemic.start = x$ci.start[2, 2], 
    i.range.x = as.numeric(c(rownames(x$data)[1], rownames(x$data)[semanas])), 
    i.range.x.53 = FALSE, 
    i.range.y = NA, 
    i.no.labels = FALSE, 
    i.start.end.marks = TRUE,
    legend.x = 29.72,
    legend.y = 0.50,
    legend.cex = 1,
    legend.y.intersp = 1.0,
    axis.text.line = 0.8
  )
  usr <- par("usr")
  x.pos <- usr[1] - x.just*(usr[2] - usr[1])
  y.pos <- usr[4] - 0.00*(usr[4] - usr[3])
}

if (TRUE) {
  opar <- par(mar = mar.b2, mgp = c(3, 0.5, 0), xpd = TRUE, family = 'serif')
  men.list <- plot.mem.(
    x = mem.model.B2,
    legend.x = 25.72,
    legend.y = 0.50,
    legend.cex = 1,
    legend.y.intersp = 1.0, 
    legend.ncol = 2,
    axis.text.line = 0.8,
    plot = FALSE
  )
  x = men.list[['x']]
  semanas = men.list[['semanas']]
  memsurveillance.(
    i.current = data.frame(rates = x$typ.curve[, 2], row.names = rownames(x$data)), 
    i.epidemic.thresholds = x$epidemic.thresholds, 
    i.intensity.thresholds = x$intensity.thresholds, 
    i.mean.length = x$mean.length, 
    i.force.length = TRUE, 
    i.graph.file = FALSE, 
    i.week.report = NA, 
    i.equal = FALSE, 
    i.pos.epidemic = TRUE,
    i.no.epidemic = FALSE, 
    i.no.intensity = FALSE,
    i.epidemic.start = x$ci.start[2, 2], 
    i.range.x = as.numeric(c(rownames(x$data)[1], rownames(x$data)[semanas])), 
    i.range.x.53 = FALSE, 
    i.range.y = NA, 
    i.no.labels = FALSE, 
    i.start.end.marks = TRUE,
    legend.x = 29.72,
    legend.y = 0.50,
    legend.cex = 1,
    legend.y.intersp = 1.0,
    axis.text.line = 0.8
  )
  usr <- par("usr")
  x.pos <- usr[1] - x.just*(usr[2] - usr[1])
  y.pos <- usr[4] - 0.00*(usr[4] - usr[3])
}

par(opar)
dev.off()