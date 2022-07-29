# NOTE: We will probably want to replicate much of what we have published, for completeness.  But unless we get different findings, I think that will be more of an exercise to check results rather than something we will focus on.

######################################################################
# Author: Magali Blanco
# Date: 7/22/2022
# Purpose: Run dementia-TRAP survival models
######################################################################

######################################################################
# SETUP
######################################################################
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, kableExtra, survival) #survminer,

set.seed(1)

source("functions.R")
output_data_path <- file.path("Data", "Output")
if(!file.exists(output_data_path)) {dir.create(output_data_path)}

image_path <- file.path("..", "manuscript", "images")
######################################################################
# LOAD DATA
######################################################################
sur <- readRDS(file.path("Data", "Output", "sur.rda")) 

######################################################################
# COMMON VARIABLES
######################################################################
pnc_units <- 1000
bc_units <- 10
no2_units <- 5
pm25_units <- 1

sensitivity_yr <- 2005

######################################################################
# UPDATE DATA
######################################################################

sur <- sur %>%
  group_by(pollutant) %>%
  mutate(
    pollutant_prediction = case_when(
      grepl("ufp", pollutant) ~ pollutant_prediction/pnc_units,
      pollutant == "bc" ~pollutant_prediction/bc_units,
      pollutant == "no2" ~pollutant_prediction/no2_units,
      pollutant == "pm25" ~pollutant_prediction/pm25_units
    )) %>%
  ungroup()

sur <- drop_na(sur, apoe) 
  
# --> WHY DO SOME COVERAGE VARIABLES HAVE NAS BUT STILL HAVE PREDICTIONS?? DROP UNNECESSARY COLUMNS?

# for multipollutant modeling
sur_w <- sur %>%
  select(-exp_avg0) %>%
  pivot_wider(names_from = c(model, pollutant), values_from = pollutant_prediction) 

main_ap_models <- c(
  # single pollutant models
  as.list(str_subset(names(sur_w), "^MM|^SP|^ST")),
  # 2 pollutant models
  # --> NEED?
  list(paste("MM", c("no2", "pm25"), sep = "_")),
  list(paste("SP", c("no2", "pm25"), sep = "_")),
  list(paste("ST", c("no2", "pm25"), sep = "_")),
  # 3 pollutant models
  list(paste("MM", c("ufp_10_42", "bc", "no2"), sep = "_")),
  # 4 pollutant models
  list(paste("MM", c("ufp_10_42", "bc", "no2", "pm25"), sep = "_"))
)
# TEST
# main_ap_models <- c("MM_bc", "MM_no2")
######################################################################
# FUNCTIONS
######################################################################

# --> also pull out person-years & ties?

# fn runs coxph() models and returns HR results
# dt = group_split(sur_w, exposure_duration)[[1]]
# event_indicator = "dementia_now"
# pollutant_predictors = c("MM_no2", "MM_bc")
# other_predictors = c("cal_2yr", "male", "race_white", "degree", "income_cat")
run_cox <- function(dt, event_indicator, pollutant_predictors, other_predictors = c("cal_2yr", "male", "race_white", "degree", "income_cat")) {
  model_description <- paste(first(dt$exposure_duration), "yr", paste(pollutant_predictors, collapse="+"))
  exposure_model <- paste(unique(str_extract(pollutant_predictors, "MM|SP|ST")), collapse = " + ")
  message(model_description)
  
  #create a survival object
  surv_object <- Surv(time = dt$age_start_exposure,
                      time2 = dt$age_end_exposure,
                      event = dt[[event_indicator]])
  
  predictors <- c(pollutant_predictors, other_predictors)
  cox_model <- coxph(as.formula(paste("surv_object ~" ,paste(predictors, collapse = "+"), "+ strata(apoe)")), 
                     data=dt,
                     cluster = study_id, #robust = T,
                     weights = model_wt) %>%
    summary()
  
  # --> this uses robust SE to calculate 95% conf int, right? 
  result <- cox_model$conf.int %>% 
    as.data.frame() %>%
    rownames_to_column(var = "covariate") %>%
    mutate(pollutant_predictors = paste(substr(pollutant_predictors, 4, str_length(pollutant_predictors)), collapse = " + "), #first(dt$pollutant),
           exposure_duration = first(dt$exposure_duration),
           exposure_model = exposure_model,
           outcome = event_indicator) %>%
    rename_all(~make.names(.)) %>%
    select(-exp..coef.) %>%
    rename(hr = exp.coef.,
           lower95 = lower..95,
           upper95 = upper..95)
  
  return(result)
  }

run_cox_many_times <- function(dt, pollutant_predictors, ...) {
  hr_results <- data.frame()
  
  for(i in c("dementia_now", "ad_now")) {
    message(paste("outcome:", i))
    
    for(m in pollutant_predictors) {
      #m <- unlist(m)
      temp <- lapply(group_split(dt, exposure_duration), run_cox, event_indicator = i, pollutant_predictors=m, ...) %>% 
        bind_rows()
      
      hr_results <- rbind(hr_results, temp)
    }}
  return(hr_results)
  }

######################################################################
# MODELS
######################################################################
# first set of models
hrs_main <- run_cox_many_times(dt = sur_w, pollutant_predictors = main_ap_models) %>%
  mutate(description = "Main Models")

# 2005+ (person-years start in 2006 since exposure is always 1 yr earlier)
hrs_more_recent <- sur_w %>%
  filter(exposure_year >= sensitivity_yr) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., pollutant_predictors = main_ap_models) %>%
  mutate(description = paste0(sensitivity_yr, "+"))
 
# no IPW
hrs_no_ipw <- sur_w %>%
  mutate(model_wt = 1) %>%
  run_cox_many_times(., pollutant_predictors = main_ap_models) %>%
  mutate(description = "No IPW")



######################################################################
# ADDITIONAL EXPLORATORY ANALYSES
######################################################################






######################################################################
# SAVE HAZARD RATIOS
######################################################################
hrs <- rbind(hrs_main, hrs_more_recent, hrs_no_ipw)
saveRDS(hrs, file.path(output_data_path, "hazard_ratios.rda"))
