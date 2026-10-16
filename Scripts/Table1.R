# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))

# Define season function
get.season <- function(month) {
  if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else if (month %in% c(9, 10, 11)) {
    return("Autumn")
  } else {
    return("Winter")
  }
}

# Define the classification function
# Young, Children, Adolescents, Youth, Middle-age, Middle-age Elderly, Elderly
classify.age <- function(age) {
  if (age >= 0 && age < 5) {
    return("[0, 5)")
  } else if (age >= 5 && age < 10) {
    return("[5, 10)")
  } else if (age >= 10 && age < 20) {
    return("[10, 20)")
  } else if (age >= 20 && age < 35) {
    return("[20, 35)")
  } else if (age >= 35 && age < 50) {
    return("[35, 50)")
  } else if (age >= 50 && age < 65) {
    return("[50, 65)")
  } else if (age >= 65) {
    return("65+")
  } else {
    return(NA)
  }
}

# instead of the threshold-style "P < 0.0001".
# Note: values below .Machine$double.eps (2.220e-16) are reported as 2.220 × 10^-16.
format_p_exact <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < .Machine$double.eps) {
      return('2.220 × 10^-16')
    }
    if (x >= 0.001) {
      return(formatC(x, format = 'f', digits = 3))
    }
    exponent <- floor(log10(x))
    sprintf('%.3f × 10^%d', x / 10^exponent, exponent)
  }, character(1))
}

# Prepare plot data
sdate <- '2020-04-01'
edate <- '2023-04-30'

# Load data
# The primary clinical dataset comprising individual-level electronic health records 
# from Kiang Wu Hospital is not publicly available due to patient privacy regulations
# and ethical restrictions regarding the protection of personal health information;
# however, anonymized data supporting the findings of this study may be made available 
# to qualified researchers upon reasonable request to the corresponding authors,
# subject to approval by the institutional review board and the execution of a data sharing agreement.
UDF <- data.table::fread(glue('{dat.dir}/Macau/INFLU-U2.csv.gz'))

# Arrange data
IDF <- UDF %>%
  dplyr::filter(FLUA != -1 & FLUB != -1) %>%
  dplyr::filter(KID >= as.Date('2010-01-01')) %>%
  dplyr::select(Local, Gender, Age, Month, FLUA, FLUB, Emergency) %>% 
  dplyr::filter(Gender %in% c('F', 'M'))

IDF[, Season := sapply(Month, get.season)]
IDF[, Age := sapply(Age, classify.age)]

dat <- IDF %>% dplyr::select(-c(Month))

dat$Age <- factor(dat$Age, levels =  c("[0, 5)", "[5, 10)", "[10, 20)", "[20, 35)", "[35, 50)", "[50, 65)", "65+"))
dat$Season <- factor(dat$Season, levels = c("Spring", "Summer", "Autumn", "Winter"))
dat$Gender <- as.factor(dat$Gender)
dat$Area <- factor(ifelse(dat$Local == 1, "Local", "Tourist"), levels = c("Local", "Tourist"))
dat$Department  <- factor(ifelse(dat$Emergency == 1, "Emergency", "Outpatient"), levels = c("Emergency", "Outpatient"))
dat$Local <- NULL

cols <- c('FLUA', 'FLUB')
vals <- c('Age', 'Gender', 'Department', 'Area', 'Season')

DF <- data.frame(Variable = 'Variable', Samples = 'Samples', FLUA.P = 'FLUA.P', FLUA.C = 'FLUA.C', FLUB.P = 'FLUB.P', FLUB.C = 'FLUB.C')
index <- 2
for (val in vals) {
  DF[index, 'Variable'] = val; index = index + 1
  
  subdat.A <- data.frame(VAL = dat[[rlang::sym(val)]], FLU = dat[['FLUA']]) %>% as_tibble()
  sumdat.A <- subdat.A %>%
    group_by(VAL) %>%
    summarise(
      Samples = n(),
      FLU = sum(FLU),
      FLUP = sprintf('%.2f', FLU / Samples * 100),
      .groups = 'drop'
    ) %>%
    na.omit()
  total.samples.A <- sum(sumdat.A$Samples)
  chisq.result.A <- chisq.test(rbind(sumdat.A$FLU, sumdat.A$Samples - sumdat.A$FLU), correct = FALSE)
  
  subdat.B <- data.frame(VAL = dat[[rlang::sym(val)]], FLU = dat[['FLUB']]) %>% as_tibble()
  sumdat.B <- subdat.B %>%
    group_by(VAL) %>%
    summarise(
      Samples = n(),
      FLU = sum(FLU),
      FLUP = sprintf('%.2f', FLU / Samples * 100),
      .groups = 'drop'
    ) %>%
    na.omit()
  total.samples.B <- sum(sumdat.B$Samples)
  chisq.result.B <- chisq.test(rbind(sumdat.B$FLU, sumdat.B$Samples - sumdat.B$FLU), correct = FALSE)
  
  for (ind in 1:nrow(sumdat.A)) {
    DF[index, 'Variable'] = as.character(sumdat.A[ind, ]$VAL)
    DF[index, 'Samples'] = as.character(sumdat.A[ind, ]$Samples)
    DF[index, 'FLUA.P'] = glue::glue('{sumdat.A[ind, ]$FLU} ({sumdat.A[ind, ]$FLUP})')
    if (ind == 1) DF[index, 'FLUA.C'] = sprintf("χ² = %.2f", chisq.result.A$statistic)
    if (ind == 2) DF[index, 'FLUA.C'] = sprintf('P = %s', format_p_exact(chisq.result.A$p.value))
    DF[index, 'FLUB.P'] = glue::glue('{sumdat.B[ind, ]$FLU} ({sumdat.B[ind, ]$FLUP})')
    if (ind == 1) DF[index, 'FLUB.C'] = sprintf("χ² = %.2f", chisq.result.B$statistic)
    if (ind == 2) DF[index, 'FLUB.C'] = sprintf('P = %s', format_p_exact(chisq.result.B$p.value))
    index = index + 1
  }
}

DF <- DF %>% mutate(across(everything(), ~ tidyr::replace_na(.x, "")))

writexl::write_xlsx(DF, glue('{tbl.dir}/Table-Information.xlsx'), col_names = FALSE)