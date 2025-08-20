library(here)
library(data.table)
library(magrittr)
library(readxl)
library(lubridate)
library(seasonal)


# function for seasonal 12-month adjustment using X-13ARIMA-SEATS
adj_seas <- function(data, varname, exclude_var = NULL) {
  full_ts <- data[, ..varname] %>%
    ts(c(year(min(data[, date])), month(min(data[, date]))),
       c(year(max(data[, date])), month(max(data[, date]))),
       12)
  
  if (!is.null(exclude_var)) {
    # Create outlier specifications for excluded observations
    excluded_dates <- data[get(exclude_var) == TRUE, date]
    
    # Convert dates to X-13 format
    outlier_specs <- paste0("AO", year(excluded_dates), ".", 
                            month.abb[month(excluded_dates)])
    
    # Run seas with outlier adjustments for excluded periods
    m <- seas(full_ts, 
              regression.variables = outlier_specs,
              outlier = NULL)
    
  } else {
    m <- seas(full_ts)
  }
  
  # extract the final adjusted series
  adjusted_series <- series(m, "s11")
  return(as.numeric(adjusted_series))
}


#### CES vintage data ####

# include industries: total nonfarm, total private, rail, insurance, real estate, hospitals, religious orgs
c_ind <- c('000000','050000','434820','555240','555310','656220','808130')

# function to read the vintage grids downloaded from https://www.bls.gov/web/empsit/cesvinall.zip
read_ces <- function(ind, sa = FALSE) {
  dt <- fread(here('data-in', paste0('tri_', ind, '_', ifelse(sa, 'SA', 'NSA'), '.csv'))) %>%
    .[, date := as.Date(paste0(year, '-', sprintf('%02d', month), '-01'))] %>%
    .[, `:=`(year = NULL, month = NULL)] %>%
    melt(id.var = 'date', variable.factor = FALSE) %>%
    .[,
      .(date,
        estimate_date = as.Date(paste0(substr(variable,1,3),
                                       ifelse(as.numeric(substr(variable,5,6)) >= 39, '19', '20'), substr(variable,5,6),
                                       '01'), '%b%Y%d'),
        ind = ind,
        value)] %>%
    .[, monthdiff := interval(estimate_date, date) %/% months(1)] %>%
    .[monthdiff %in% 0:3]
  
  return(dt)
}

# CES total nonfarm minus industries with high rates of non-covered workers
dt_ces <- lapply(c_ind, read_ces, sa = TRUE) %>%
  rbindlist() %>%
  dcast(date + estimate_date + monthdiff ~ paste0('ces_', ind), value.var = 'value') %>%
  .[,
    .(date,
      estimate_date,
      monthdiff,
      ces_tot = ces_000000 - ces_434820 - ces_555240 - ces_555310 - ces_656220 - ces_808130,
      ces_prv = ces_050000 - ces_434820 - ces_555240 - ces_555310 - ces_656220 - ces_808130)] %>%
  dcast(date ~ monthdiff, value.var = c('ces_tot','ces_prv')) %>%
  .[,
    .(date,
      ces_tot_otm = ces_tot_0 - ces_tot_1,
      ces_prv_otm = ces_prv_0 - ces_prv_1)]


#### ADP data ####
# vintages for two ADP series downloaded from ALFRED
dt_adp1 <- fread(here('data-in','NPPTTL_2','vintages_starting_2011-02-02.csv')) %>%
  melt(id.var = 'observation_date', variable.factor = FALSE, variable.name = 'vintage', value.name = 'adp_estimate') %>%
  .[, vintage := gsub('NPPTTL_','',vintage, fixed = TRUE)] %>%
  .[, adp_estimate := adp_estimate * 1000]

dt_adp2 <- fread(here('data-in','ADPMNUSNERSA_2','vintages_starting_2022-08-31.csv')) %>%
  melt(id.var = 'observation_date', variable.factor = FALSE, variable.name = 'vintage', value.name = 'adp_estimate') %>%
  .[, vintage := gsub('ADPMNUSNERSA_','',vintage, fixed = TRUE)]

# combine and select earliest two vintages for each month and calculate otm value
dt_adp <- rbind(dt_adp1, dt_adp2) %>%
  .[order(vintage, -observation_date)] %>%
  .[!is.na(adp_estimate)] %>%
  .[, sort := 0:(.N-1), vintage] %>%
  .[sort %in% 0:1] %>%
  .[, release_date := max(observation_date), vintage] %>%
  dcast(release_date ~ paste0('adp_est', sort), value.var = 'adp_estimate', fun.aggregate = min) %>%
  .[adp_est0 > 0 & adp_est1 > 0] %>%
  .[,
    .(date = release_date,
      adp = adp_est0 / 1000,
      adp_otm = (adp_est0 - adp_est1) / 1000)]

#### QCEW data ####
# series selected with the BLS series report tool
dt_qcew <- read_excel(here('data-in','SeriesReport-20250813205929_1ddb77.xlsx'), skip = 2) %>%
  data.table() %>%
  melt(id.var = 'Series ID', variable.factor = FALSE) %>%
  dcast(variable ~ `Series ID`, value.var = 'value') %>%
  .[,
    .(date = as.Date(paste0(variable, '01'), '%b\n%Y%d'),
      qcew_tot = (ENUUS00010010 - ENUUS00010511 + ENUUS0001051133 - ENUUS000105482 - ENUUS000105524 - ENUUS000105531 - ENUUS000105622 - ENUUS000105813 - ENUUS000105814) / 1000,
      qcew_prv = (ENUUS00010510 - ENUUS00010511 + ENUUS0001051133 - ENUUS000105482 - ENUUS000105524 - ENUUS000105531 - ENUUS000105622 - ENUUS000105813 - ENUUS000105814) / 1000)] %>%
  .[!is.na(date)] %>%
  .[order(date)]

# seasonally adjust QCEW series, manually excluding 2020 (auto outlier detection didn't seem to work)
dt_qcew[, covid := year(date) == 2020]
dt_qcew[, qcew_tot_sa := adj_seas(dt_qcew, 'qcew_tot', 'covid')]
dt_qcew[, qcew_tot_otm_sa := qcew_tot_sa - shift(qcew_tot_sa)]

# Seasonally adjust the private series
dt_qcew[, qcew_prv_sa := adj_seas(dt_qcew, 'qcew_prv', 'covid')]
dt_qcew[, qcew_prv_otm_sa := qcew_prv_sa - shift(qcew_prv_sa)]

# combine the three files
dt_all <- dt_qcew %>%
  dt_ces[., on = 'date'] %>%
  dt_adp[., on = 'date']

# file with data for dots
fwrite(dt_all, here('data-out', 'ces-accuracy.csv'))


#### regressions ####
# regression of CES on QCEW
mod1 <- dt_all %>%
  .[!is.na(ces_tot_otm) & year(date) != 2020] %>%
  lm(qcew_tot_otm_sa ~ ces_tot_otm, data = .) %>%
  summary()

dt_mod1 <- data.table(beta0 = mod1$coefficients[1,1],
                      beta1 = mod1$coefficients[2,1],
                      r2 = mod1$r.squared,
                      n = mod1$df[1] + mod1$df[2] + 1)

fwrite(dt_mod1, here('data-out', 'ces-accuracy-mod1.csv'))

# regression of CES on QCEW where ADP is available
mod2a <- dt_all %>%
  .[!is.na(adp_otm) & year(date) != 2020] %>%
  lm(qcew_prv_otm_sa ~ ces_prv_otm, data = .) %>%
  summary()

dt_mod2a <- data.table(beta0 = mod2a$coefficients[1,1],
                       beta1 = mod2a$coefficients[2,1],
                       r2 = mod2a$r.squared,
                       n = mod2a$df[1] + mod1$df[2] + 1)

fwrite(dt_mod2a, here('data-out', 'ces-accuracy-mod2a.csv'))

# regression of ADP on QCEW
mod2b <- dt_all %>%
  .[!is.na(adp_otm) & year(date) != 2020] %>%
  lm(qcew_prv_otm_sa ~ adp_otm, data = .) %>%
  summary()

dt_mod2b <- data.table(beta0 = mod2b$coefficients[1,1],
                       beta1 = mod2b$coefficients[2,1],
                       r2 = mod2b$r.squared,
                       n = mod2b$df[1] + mod1$df[2] + 1)

fwrite(dt_mod2b, here('data-out', 'ces-accuracy-mod2b.csv'))



# line chart of absolute errors
dt_line <- dt_all %>%
  .[!is.na(ces_tot_otm),
    .(date,
      ces_error = ces_tot_otm - qcew_tot_otm_sa,
      ces_abs_error = abs(ces_tot_otm - qcew_tot_otm_sa),
      adp_error = adp_otm - qcew_tot_otm_sa,
      adp_abs_error = abs(adp_otm - qcew_tot_otm_sa))] %>%
  .[, ces_error_6m := frollmean(ces_error, 6)] %>%
  .[, ces_abs_error_6m := frollmean(ces_abs_error, 6)] %>%
  .[, adp_error_6m := frollmean(adp_error, 6)] %>%
  .[, adp_abs_error_6m := frollmean(adp_abs_error, 6)]

fwrite(dt_line, here('data-out', 'ces-accuracy-line.csv'))


# Bias by presidential term
dt_all %>%
  .[, pres := 'Other'] %>%
  .[year(date) >= 2001, `:=`(pres = 'Bush', party = 'R')] %>%
  .[year(date) >= 2009, `:=`(pres = 'Obama', party = 'D')] %>%
  .[year(date) >= 2017, `:=`(pres = 'Trump (1st)', party = 'R')] %>%
  .[year(date) >= 2021, `:=`(pres = 'Biden', party = 'D')] %>%
  .[year(date) >= 2025, `:=`(pres = 'Trump (2nd)', party = 'R')] %>%
  .[!is.na(ces_tot_otm),
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_qcew_otm_sa = sum(qcew_tot_otm_sa) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N),
    .(party,
      pres)] %>%
  fwrite(here('data-out', 'ces-accuracy-pres.csv'))


# get some avg errors for different time periods
dt_all %>%
  .[!is.na(ces_tot_otm) & year(date) != 2020,
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N)]

dt_all %>%
  .[!is.na(ces_tot_otm),
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N)]

dt_all %>%
  .[!is.na(adp_otm) & year(date) != 2020,
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N)]

dt_all %>%
  .[!is.na(adp_otm) & year(date) != 2020,
    .(avg_error = sum(adp_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(adp_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N)]

dt_all %>%
  .[!is.na(ces_tot_otm),
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N)]

dt_all %>%
  .[!is.na(ces_tot_otm),
    .(avg_error = sum(ces_tot_otm - qcew_tot_otm_sa) / .N,
      avg_abs_error = sum(abs(ces_tot_otm - qcew_tot_otm_sa)) / .N,
      avg_abs_qcew_otm_sa = sum(abs(qcew_tot_otm_sa)) / .N),
    .(year(date) %in% 2020:2022)]

