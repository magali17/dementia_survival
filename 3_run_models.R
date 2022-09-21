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

# citation("survival")
#knitr::write_bib(c(.packages(), "survival"), file.path("~", "Desktop", "packages.ris"))


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
bc_units <- 100
no2_units <- 5
pm25_units <- 1

save(pnc_units, bc_units, no2_units, pm25_units, file = file.path(output_data_path, "trap_model_units.rda"))

main_exposure_duration <- 10
first_exposure_year <- load(file.path(output_data_path, "first_exposure_year.rda")) #2005

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
  list(c("MM_ufp_10_42", "SP_pm25")),
  list(c("MM_ufp_10_70", "SP_pm25")),
  list(c("MM_ufp_20_1k", "SP_pm25")),
  list(c("MM_ufp_36_1k", "SP_pm25")),
  list(c("MM_bc", "SP_pm25")),
  list(c("MM_no2", "SP_pm25")),
  
  # 3 pollutant models
  #list(paste("MM", c("ufp_10_42", "bc", "no2"), sep = "_")),
  # 4 pollutant models
  list(c(paste("MM", c("ufp_10_42", "bc", "no2"), sep = "_"), "SP_pm25"))
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
run_cox <- function(dt, 
                    start_time = "age_start_exposure", end_time = "age_end_exposure", event_indicator, 
                    pollutant_predictors, 
                    other_predictors = c("cal_2yr", "male", "race_white", "degree", "income_cat")) {
  
  model_description <- paste(first(dt$exposure_duration), "yr", paste(pollutant_predictors, collapse="+"))
  message(model_description)
  exposure_model <- paste(unique(str_extract(pollutant_predictors, "MM|SP|ST")), collapse = " + ")
  
  
  #create a survival object
  surv_object <- Surv(time = dt[[start_time]],
                      time2 = dt[[end_time]],
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

run_cox_many_times <- function(dt, pollutant_predictors = main_ap_models, ...) {
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
hrs_main <- sur_w %>%
  filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(.) %>%
  mutate(description = paste0("Main Models (", first_exposure_year, "+)"))

# all years
hrs_full_cohort <- sur_w %>%
  filter(exposure_duration == main_exposure_duration) %>%
  run_cox_many_times(.) %>%
  mutate(description = "Full Cohort (1994+)")

# fewer years
hrs_2010 <- sur_w %>%
  filter(exposure_year >= first_exposure_year+5,
         exposure_duration == main_exposure_duration
         ) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(.) %>%
  mutate(description = paste0(first_exposure_year+5, "+"))
 
# no IPW
hrs_no_ipw <- sur_w %>%
  filter(exposure_year >= first_exposure_year,
         exposure_duration == main_exposure_duration
         ) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  mutate(model_wt = 1) %>%
  run_cox_many_times(.) %>%
  mutate(description = "No IPW")

# adjust for intake year
hrs_intake_age <- sur_w %>%
  filter(exposure_year >= first_exposure_year,
         exposure_duration == main_exposure_duration
         ) %>% 
  mutate(cal_2yr = droplevels(cal_2yr),
         intakeage = cut(intakeage, seq(min(intakeage), max(intakeage), 5), right = F)
         ) %>%  
  run_cox_many_times(., other_predictors = c("intakeage", "cal_2yr", "male", "race_white", "degree", "income_cat")) %>%
  mutate(description = "Intake Age adjustment")

# calendar time-axis, adjust 2-year age group
# --> did this correctly? see Rachel's work
hrs_calendar_axis <- sur_w %>%
  filter(exposure_year >= first_exposure_year,
         exposure_duration == main_exposure_duration
         ) %>%  
  mutate(cal_2yr = droplevels(cal_2yr),
         exposure_year0 = exposure_year-355/356, #start is 1 day into the previous year. may be fine to just use 1 for simplicity
         age_start_exposure = cut(age_start_exposure, seq(min(age_start_exposure), max(age_start_exposure), 2), right = F)
         ) %>% 
  run_cox_many_times(., start_time = "exposure_year0", end_time = "exposure_year", 
                     other_predictors = c("age_start_exposure", "cal_2yr", "male", "race_white", "degree", "income_cat")
                     ) %>%
  mutate(description = "Time Axis, 2yr age adj")


#--> also stratify by sex. or make this the main model?
# --> also adjust for marital/living status or social engagement

######################################################################
# INTERACTION MODELS
######################################################################

## --> 




######################################################################
# SAVE HAZARD RATIOS
######################################################################
hrs <- rbind(hrs_main, hrs_full_cohort, hrs_2010, hrs_no_ipw, hrs_intake_age, hrs_calendar_axis)
saveRDS(hrs, file.path(output_data_path, "hazard_ratios.rda"))
