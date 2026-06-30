# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(splines)))
Sys.setlocale('LC_TIME', 'C')

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
n.years <- 15
test.sta <- FALSE

# Load data
dat.raw <- readxl::read_excel(glue('{dat.dir}/Macau/FLU-CL-AQ-day.xlsx'))
dat.raw$dif.temp <- dat.raw$max.temp - dat.raw$min.temp

# Prepare plot data
envs <- c('ave.temp', 'dif.temp', 'max.temp', 'min.temp', 'min.temp', 'humidity', 'sunshine')

for (env in envs) {
  
  if (flu == 'FLUA') { nTest = 'nFLUA' } else { nTest = 'nFLUB' }
  
  dat.fac <- dat.raw %>% dplyr::select(date, holiday, all_of(env), nSample, all_of(nTest)) %>%
    mutate(nNegative = nSample - get(nTest)) %>%
    mutate(DOW = weekdays(date))

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Model Fitting
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cb.temp <- crossbasis(dat.fac[[env]], lag = 21, argvar = list(fun = 'ns', df = 4), arglag = list(fun = 'ns', df = 3))
  
  # Build a Generalized Linear Model (GLM)
  model.binom <- glm(cbind(get(nTest), nNegative) ~ cb.temp + DOW + holiday + ns(date, df = 7 * n.years), family = quasibinomial(), data = dat.fac)
  
  # Prediction & Plotting
  cen.temp <- median(dat.fac[[env]], na.rm = TRUE)
  pred.binom <- crosspred(cb.temp, model.binom, cen = cen.temp, by = 0.5)
  
  # Overall Cumulative Association
  if (test.sta) layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
  
  # Overall Effect of Temperature on Positive Rate
  if (!test.sta) pdf(glue('{fig.dir}/{ofig}_{flu}_{env}_splines.pdf'), width = 5, height = 5, family = 'serif', bg = '#FFFFFF')
  par(
    mar = c(3, 3, 1, 1),
    cex.axis = 1.0,
    cex.lab = 1.2,
    font.lab = 2,
    mgp = c(1.5, .5, 0)
  )
  plot(
    pred.binom, 'overall', 
    xlab = ifelse(env == 'humidity', 'Humidity (%)', ifelse(env == 'sunshine', 'Hour (h)', 'Temperature (°C)')), 
    ylab = 'Odds Ratio (Positive Rate)', 
    col = 'red', lwd = 2, 
    ci.arg = list(col = gray(0.8), density = NA)
  )
  
  # Add reference lines to a plot
  abline(h = 1, lty = 2)
  abline(v = cen.temp, lty = 2, col = 'blue')
  if (env == 'humidity') {
    text(cen.temp, 0.8, labels = bquote(paste('Ref Humidity: ', .(sprintf('%.1f', cen.temp)), ' (%)')), col = 'blue')
  } else if (env == 'sunshine') {
    text(cen.temp, 0.8, labels = bquote(paste('Ref Hour: ', .(sprintf('%.1f', cen.temp)), ' (h)')), col = 'blue')
  }  else {
    text(cen.temp, 0.8, labels = bquote(paste('Ref Temp: ', .(sprintf('%.1f', cen.temp)), ' ', degree, 'C')), col = 'blue')
  }
  if (!test.sta) dev.off()
  
  # 3D Exposure-Lag-Response Surface
  if (!test.sta) pdf(glue('{fig.dir}/{ofig}_{flu}_{env}_3d.pdf'), width = 5, height = 5, family = 'serif', bg = '#FFFFFF')
  par(
    mar = c(1, 1, 1, 1),
    cex.axis = 1.0,
    cex.lab = 1.2,
    font.lab = 2,
    tcl = -0.3
  )
  plot(
    pred.binom, 
    xlab = ifelse(env == 'humidity', 'Humidity (%)', ifelse(env == 'sunshine', 'Hour (h)', 'Temperature (°C)')),
    zlab = '\nOR', 
    ylab = 'Lag (Days)', 
    theta = 210, 
    phi = 30, 
    ltheta = -120
  )
  if (!test.sta) dev.off()
  
  # Reset layout
  if (test.sta) layout(1)

  # Print model summary
  summary(model.binom)
}


