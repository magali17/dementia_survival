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
#sur <- readRDS(file.path("Data", "Output", "sur.rda")) 
sur0 <- readRDS(file.path("Data", "Output", "sur_full_cohort.rda")) %>%
  filter(!is.na(apoe))

######################################################################
# COMMON VARIABLES
######################################################################
pnc_units <- 1000
bc_units <- 100
no2_units <- 5
pm25_units <- 1

save(pnc_units, bc_units, no2_units, pm25_units, file = file.path(output_data_path, "trap_model_units.rda"))

main_exposure_duration <- 10
#this is only being used for sensitivity analyses now
load(file.path(output_data_path, "1_prep_data", "first_exposure_year.rda")) 

freeze_date <- readRDS(file.path(output_data_path, "freeze_date.rda"))
######################################################################
# UPDATE DATA
######################################################################

sur <- sur0 %>%
  group_by(pollutant) %>%
  mutate(
    pollutant_prediction = case_when(
      grepl("ufp|ns_|pnc", pollutant) ~ pollutant_prediction/pnc_units,
      pollutant == "bc" | grepl("_bc", pollutant) ~ pollutant_prediction/bc_units,
      #grepl("bc", pollutant) ~ pollutant_prediction/bc_units,
      pollutant == "no2" ~pollutant_prediction/no2_units,
      pollutant == "pm25" ~pollutant_prediction/pm25_units, #TRUE ~ pollutant
      TRUE ~ pollutant_prediction
    )) %>%
  ungroup()

#sur <- drop_na(sur, apoe) 
  
# --> WHY DO SOME COVERAGE VARIABLES HAVE NAS BUT STILL HAVE PREDICTIONS?? DROP UNNECESSARY COLUMNS?

# for multipollutant modeling
sur_w <- sur %>%
  select(-exp_avg0) %>%
  pivot_wider(names_from = c(model, pollutant), values_from = pollutant_prediction) 
save(sur_w, file = file.path("Data", "Output", "sur_w.rda"))


#outcomes <- str_subset(names(sur_w), "_now$" )
outcomes <- c("dementia_now", "ad_now", "non_ad_now")

# 2 pollutant models
basic_ap_models <- c(
  list(c("MM_ufp_10_42", "SP_pm25")),
  list(c("MM_bc", "SP_pm25")),
  list(c("MM_no2", "SP_pm25")))


all_ap_models <- c(
  # single pollutant models
  as.list(c("MM_ufp_10_42", "MM_bc", "MM_no2")),
  # compare to Rachel's findings
  list(c("ST_pm25")),
  # other SP 2.5 on its own
  list(c("SP_pm25")),
  
  # 2 pollutant models
  basic_ap_models,
  list(c("MM_ufp_10_70", "SP_pm25")),
  list(c("MM_ufp_20_1k", "SP_pm25")),
  
  #bins
  list(c("MM_ns_11.5", "SP_pm25")),
  list(c("MM_ns_10_100", "SP_pm25")),
  list(c("MM_ns_86.6", "SP_pm25")),
  list(c("MM_ns_64.9", "SP_pm25")),
  list(c("MM_ns_48.7", "SP_pm25")),
  list(c("MM_ns_154.0", "SP_pm25")),
  list(c("MM_ns_36.5", "SP_pm25")),
  list(c("MM_ns_27.4", "SP_pm25")),
  list(c("MM_ns_20.5", "SP_pm25")),
  list(c("MM_ns_15.4", "SP_pm25")),
  list(c("MM_ns_115.5", "SP_pm25")),
  
  #ptrak diff
  list(c("MM_pnc_20_36", "SP_pm25")),
  list(c("MM_ufp_36_1k", "SP_pm25")),
  
  #pt size
  list(c("MM_pmdisc_sz", "SP_pm25")),
  
  #onroad
  list(c("MM_pnc_onrd", "SP_pm25")),
  

  list(c("SP_no2", "SP_pm25")),
  list(c("ST_no2", "SP_pm25")),
  
  # 4 pollutant models
  list(c(paste("MM", c("ufp_10_42", "bc", "no2"), sep = "_"), "SP_pm25")))

# TEST
# all_ap_models <- c("MM_bc", "SP_pm25")
######################################################################
# FUNCTIONS
######################################################################

# --> also pull out person-years & ties?

# fn runs coxph() models and returns HR results
# dt = group_split(sur_w, exposure_duration)[[1]]
# event_indicator = "dementia_now"
# pollutant_predictors = c("MM_no2", "SP_pm25")
# other_predictors = c("cal_2yr", "degree", "nses_z_cx", "male", "race_white")
# stratification = c("apoe")
# start_time = "age_start_exposure"
# end_time = "age_end_exposure"

run_cox <- function(dt, 
                    start_time = "age_start_exposure", end_time = "age_end_exposure", event_indicator, 
                    pollutant_predictors, 
                    other_predictors = c("cal_2yr", "degree", "nses_z_cx", "male", "race_white"), #c("cal_2yr", "degree", "nses_z_cx"),
                    stratification = c("apoe") #c("apoe", "male", "race_white")
                    ) {
  
  model_description <- paste(first(dt$exposure_duration), "yr", paste(pollutant_predictors, collapse="+"))
  message(model_description)
  exposure_model <- paste(unique(str_extract(pollutant_predictors, "MM|SP|ST")), collapse = " + ")
  
  
  #create a survival object
  surv_object <- Surv(time = dt[[start_time]],
                      time2 = dt[[end_time]],
                      event = dt[[event_indicator]])
  strata_predictors <- paste0("strata(", stratification, ")", collapse = "+")
  
  predictors <- c(pollutant_predictors, other_predictors)
  cox_model <- coxph(as.formula(paste("surv_object ~" ,paste(predictors, collapse = "+"), 
                                      "+", strata_predictors
                                      #"+ strata(apoe) + strata(male) + strata(race_white)"
                                      )), 
                     data=dt,
                     cluster = study_id, #robust = T,
                     weights = model_wt) 
  print(cox_model)
  
  cox_model <- cox_model %>%
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

# run a cox model for each outcome, pollutant predictor set, and exposure duration combination
run_cox_many_times <- function(dt, pollutant_predictors = basic_ap_models, outcomes.="dementia_now", exposure_duration. = main_exposure_duration, ...) {
  hr_results <- data.frame()
  
  dt <- filter(dt, exposure_duration %in% exposure_duration.)
  
  for(i in outcomes.) {
    message(paste("outcome:", i))
    
    for(m in pollutant_predictors) {
      #m <- unlist(m)
      temp <- lapply(group_split(dt, exposure_duration), run_cox, event_indicator = i, pollutant_predictors=m, ...) %>% 
        bind_rows()
      
      hr_results <- rbind(hr_results, temp)
    }}
  return(hr_results)
  }

save(run_cox, run_cox_many_times, file = file.path("Data","Output", "cox_model_fns.rda"))

######################################################################
# MODELS
######################################################################
hrs_main <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., pollutant_predictors = all_ap_models)  %>%
  mutate(description = ifelse(grepl("pnc_20_36|ufp_36_1k|ns_", pollutant_predictors), "PNC by Size",
                              ifelse(grepl("pmdisc_sz", pollutant_predictors), "Particle Size",
                                     ifelse(grepl("pnc_onrd", pollutant_predictors), "PNC, Onroad Model", "Main Analysis"
                                            ))))  

# other outcomes
hrs_other_outcomes <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., outcomes.= c("ad_now", "non_ad_now")) %>%
  mutate(description = paste0("Dementia Subtype Outcomes"))

# other exposure periods
hrs_other_exposure_periods <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., exposure_duration. = c(1,5)) %>%
  mutate(description = paste0("Shorter Exposure Period"))

# # all years
# hrs_extended_cohort <- sur_w %>%
#   run_cox_many_times(.) %>%
#   mutate(description = "Extended Cohort (1994+)")

# most recent years
hrs_2010 <- sur_w %>%
  filter(exposure_year < first_exposure_year+10) %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%
  run_cox_many_times(.) %>%
  mutate(description = paste0("Restricted Cohort (", first_exposure_year+10, "+)"))

# drop last 2 years w/ administrative censoring - incidence is artificially high b/c person w/o visits for 
# that particular year (~50%?) are not included in the denominator
start_of_bad_yrs <- 2019 

hrs_drop_last_2yr <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%
  run_cox_many_times(.) %>%
  mutate(description = paste0("Restricted Cohort (Pre ",start_of_bad_yrs,")"))

 
# no IPW
hrs_no_ipw <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  mutate(model_wt = 1) %>%
  run_cox_many_times(.) %>%
  mutate(description = "No IPW")

# adjust for baseline age
hrs_baseline_age <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>% 
  mutate(cal_2yr = droplevels(cal_2yr)) %>%  
  group_by(study_id) %>%
  mutate(baseline_age = min(age_start_exposure)) %>% 
  ungroup() %>%
  mutate(baseline_age = cut(baseline_age, seq(min(baseline_age), max(baseline_age), 5), right = F)) %>%
  run_cox_many_times(., other_predictors = c("intakeage", "cal_2yr", "degree", "nses_z_cx", "male", "race_white")) %>% 
  mutate(description = "Baseline Age Adjustment")


# --> note, this produces very unstable age range estimates, probably b/c of very small sample sizes

# calendar time-axis, adjust 2-year age group
hrs_calendar_axis <- sur_w %>%
  #filter(exposure_year >= first_exposure_year) %>%  
  mutate(cal_2yr = droplevels(cal_2yr),
         exposure_year0 = exposure_year-355/356, #start is 1 day into the previous year. may be fine to just use 1 for simplicity
         age_start_exposure = cut(age_start_exposure, seq(min(age_start_exposure), max(age_start_exposure), 2), right = F)) %>% 
  run_cox_many_times(., start_time = "exposure_year0", end_time = "exposure_year", 
                     other_predictors = c("age_start_exposure", "cal_2yr", "degree", "nses_z_cx", "male", "race_white")) %>%  
  mutate(description = "Calendar Time Axis & 2yr Age Adjustment")

######################################################################
# INTERACTION MODELS
######################################################################

# apoe, male, race_white, degree (non/GED/HS/other vs BS/MS,PhD), NSES (nses_z_cx)
var <- c("male", "race_white", "apoe", "high_degree", "low_ndi")

effect_modification_ap_models <- c(
  list(c("MM_ufp_10_42_var", "MM_ufp_10_42_not_var", "SP_pm25")),
  list(c("MM_bc_var", "MM_bc_not_var", "SP_pm25")),
  list(c("MM_no2_var", "MM_no2_not_var", "SP_pm25")))

hrs_effect_modification <- lapply(var, function(i){
   sur_w %>%
    #filter(exposure_year >= first_exposure_year) %>% 
    mutate(cal_2yr = droplevels(cal_2yr),
           
           high_degree = ifelse(degree %in% c(0:2,6), 0, ifelse(degree %in% c(3:5), 1, NA)), #) %>% select(degree, high_degree) %>% group_by(high_degree) %>% summarize(paste(unique(degree)))
           #high_income = ifelse(income_cat %in% c(1:2), 0, ifelse(income_cat %in% c(3:4), 1, NA)),
           # NDI z-score can be interpreted as the tract-level disadvantage relative to other tracts that year
           low_ndi = ifelse(nses_z_cx < 0, 1, 0),
           
           var = as.numeric(as.character(get(i))),
           
           MM_ufp_10_42_var = var*MM_ufp_10_42,
           MM_ufp_10_42_not_var = (1-var)*MM_ufp_10_42,
           MM_bc_var = var*MM_bc,
           MM_bc_not_var = (1-var)*MM_bc,
           MM_no2_var = var*MM_no2,
           MM_no2_not_var = (1-var)*MM_no2,
    ) %>%  
    run_cox_many_times(., pollutant_predictors = effect_modification_ap_models) %>%
    mutate(description = paste(i, "Interaction"),
           #covariate = gsub("_var", paste0("_", i), covariate)
           )
  }) %>%
  bind_rows()



# same as above but to get the p-value

# --> need to recode to just pull out p-val, not HR & CI

hrs_em_pval_ap_models <- c(
  list(c("MM_ufp_10_42*var", "SP_pm25")),
  list(c("MM_bc*var", "SP_pm25")),
  list(c("MM_no2*var", "SP_pm25")))

hrs_effect_modification_pval <- lapply(var, function(i){
  sur_w %>%
    #filter(exposure_year >= first_exposure_year) %>%
    mutate(cal_2yr = droplevels(cal_2yr),
           high_degree = ifelse(degree %in% c(0:2,6), 0, ifelse(degree %in% c(3:5), 1, NA)), #) %>% select(degree, high_degree) %>% group_by(high_degree) %>% summarize(paste(unique(degree)))
           #high_income = ifelse(income_cat %in% c(1:2), 0, ifelse(income_cat %in% c(3:4), 1, NA)),
           low_ndi = ifelse(nses_z_cx < 0, 1, 0),
           
           var = as.numeric(as.character(get(i)))
           ) %>%
    run_cox_many_times(., pollutant_predictors = hrs_em_pval_ap_models) %>%
    mutate(description = paste(i, "Interaction P-Value"),
           #covariate = gsub("var", i, covariate)
           )
}) %>%
  bind_rows()
 
######################################################################
# UFP SIZE -  EXPLORATORY
######################################################################
  
# # adjust for ufp size # - don't see anything when we linearly adjsut
# hrs_ufp_size <- sur_w %>%
#   mutate(cal_2yr = droplevels(cal_2yr),
#          MM_pmdisc_sz = cut(MM_pmdisc_sz, seq(min(MM_pmdisc_sz, na.rm = T), max(MM_pmdisc_sz, na.rm = T), length.out=4))
#          
#          ) %>%   
#   group_by(study_id) %>%
#   ungroup() %>%
#   run_cox_many_times(., other_predictors = c("intakeage", "cal_2yr", "degree", "nses_z_cx", "male", "race_white", "MM_pmdisc_sz")) %>% 
#   mutate(description = "Adjusted for Pt Size")

######################################################################
# SAVE HAZARD RATIOS
######################################################################
hrs <- rbind(hrs_main, hrs_other_outcomes, hrs_other_exposure_periods, #hrs_extended_cohort, 
             hrs_no_ipw, hrs_baseline_age, hrs_calendar_axis, hrs_2010,
             hrs_drop_last_2yr,
             hrs_effect_modification, hrs_effect_modification_pval)
saveRDS(hrs, file.path(output_data_path, "hazard_ratios.rda"))

message("done with 3_run_models.R")
