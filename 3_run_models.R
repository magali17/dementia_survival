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
sur0 <- readRDS(file.path("Data", "Output", "sur_full_cohort.rda")) %>%
  filter(!is.na(apoe))

# pm25_coverage <- 0.95
sur0_pm25 <- readRDS(file.path("Data", "Output", "sur_full_cohort_no_coverage_restriction.rda")) %>%
  filter(!is.na(apoe),
         
         # -->this causes NAs in MM below.
         
         # --> relabel/organize original survival dataset to long format???
         
         #exp_wks_coverage10_yr_ST >= pm25_coverage,
         #exp_avg0=="exp_avg10_yr_ST"
         )
# sur0_pm25 %>% filter(study_id == first(study_id)) %>% View()

######################################################################
# COMMON VARIABLES
######################################################################
pnc_units <- 2000 #1000
bc_units <- 100
no2_units <- 2 #5
pm25_units <- 1

save(pnc_units, bc_units, no2_units, pm25_units, file = file.path(output_data_path, "trap_model_units.rda"))

main_exposure_duration <- 10
#this is only being used for sensitivity analyses now
load(file.path(output_data_path, "1_prep_data", "first_exposure_year.rda")) 

freeze_date <- readRDS(file.path(output_data_path, "1_prep_data", "freeze_date.rda"))

# want to drop last 2 years w/ administrative censoring - incidence is artificially high b/c person w/o visits for 
start_of_bad_yrs <- readRDS(file.path("Data", "Output", "1_prep_data", "start_of_bad_yr.rda"))
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
  select(-exp_avg0,
         -ends_with(c("yr_MM", "yr_SP", "yr_ST"))  #this causes issues when making wide otherwise
         ) %>%
  pivot_wider(names_from = c(model, pollutant), values_from = pollutant_prediction) 

# 7/21/23: upated this b/c was having issues loading in Rmd
#save(sur_w, file = file.path("Data", "Output", "sur_w.rda"))
saveRDS(sur_w, file = file.path("Data", "Output", "sur_w.rda"))



# sur_w %>%
#   filter(exposure_duration==10,
#          exposure_year < start_of_bad_yrs
#          ) %>%  
#   ggplot(aes(
#     x=ST_pm25, y=MM_bc,
#     #x=SP_pm25, y= ST_pm25, 
#              col=exposure_year, group= exposure_year)) +
#   geom_point(alpha=0.1) +
#   geom_smooth() #+ geom_abline(slope = 1, intercept = 0)


# # correlation between 10 yr MM BC & ST PM2.5 exposure. range ~ 0.35-0.5, so not super high. there's no temporal trend
# sur_w %>%
#   filter(exposure_duration==10,
#          exposure_year < start_of_bad_yrs
#   ) %>% 
#   group_by(exposure_year) %>%
#   summarize(
#     r = cor(ST_pm25, MM_bc)
#   ) %>%
#   ggplot(aes(x=exposure_year, y=r)) + 
#   geom_point() + 
#   geom_smooth()


######################################################################
# SAME FOR PM2.5
sur_pm25 <- sur0_pm25 %>%
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

# for multipollutant modeling

sur_w_pm25 <- sur_pm25 %>%
  select(-exp_avg0,
         -ends_with(c("yr_MM", "yr_SP", "yr_ST"))  #this causes issues when making wide otherwise
         ) %>%
  pivot_wider(names_from = c(model, pollutant), values_from = pollutant_prediction) 

# sur_w_pm25 %>% filter(exposure_duration==10) %>% View()

saveRDS(sur_w_pm25, file = file.path("Data", "Output", "sur_w_pm25.rda"))


######################################################################
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

# PM2.5 models; compare to Rachel's findings
st_pm2.5_models <- c(
  list(c("ST_pm25")), # could drop this since already doing
  # 2 pollutant models
  list(c("MM_ufp_10_42", "ST_pm25")),
  list(c("MM_bc", "ST_pm25")),
  list(c("MM_no2", "ST_pm25")),
  list(c("MM_ufp_10_42", "MM_bc", "MM_no2", "ST_pm25"))
  )



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
                    other_predictors = c("cal_2yr", "degree", "nses_z_cx", "male", "race_white"), 
                    stratification = c("apoe") 
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
                                      "+", strata_predictors)), 
                     data=dt,
                     cluster = study_id, 
                     weights = model_wt) 
  print(cox_model)
  
  cox_model <- cox_model %>%
    summary()
  
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

#save(run_cox, run_cox_many_times, file = file.path("Data","Output", "cox_model_fns.rda"))
save(run_cox, run_cox_many_times, file = file.path("Data","Output", "cox_model_fns.RData"))

######################################################################
# MODELS
######################################################################
hrs_main <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%   
  run_cox_many_times(., pollutant_predictors = all_ap_models)  %>%
  mutate(description = ifelse(grepl("pnc_20_36|ufp_36_1k|ns_", pollutant_predictors), "PNC by Size",
                              ifelse(grepl("pmdisc_sz", pollutant_predictors), "Particle Size",
                                     ifelse(grepl("pnc_onrd", pollutant_predictors), "PNC, Onroad Model", "Main Analysis"
                                            ))))  

hrs_dont_drop_last_2yr <- sur_w %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%
  run_cox_many_times(.) %>%
  mutate(#description = paste0("Restricted Cohort (Pre ",start_of_bad_yrs,")")
    description = "Non-Restricted Cohort (1994-2020)"
         )

# other outcomes
hrs_other_outcomes <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., outcomes.= c("ad_now", "non_ad_now")) %>%
  mutate(description = paste0("Dementia Subtype Outcomes"))

# other exposure periods
hrs_other_exposure_periods <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., exposure_duration. = c(1,5)) %>%
  mutate(description = paste0("Shorter Exposure Period"))

# most recent years
hrs_2010 <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs,
         exposure_year < first_exposure_year+10) %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%
  run_cox_many_times(.) %>%
  mutate(description = paste0("Restricted Cohort (", first_exposure_year+10, "+)")
         #description = paste0("Restricted Cohort (", first_exposure_year+10, "-", start_of_bad_yrs-1, ")")
         )

 
# no IPW
hrs_no_ipw <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  mutate(model_wt = 1) %>%
  run_cox_many_times(.) %>%
  mutate(description = "No IPW")

# adjust for baseline age
hrs_baseline_age <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
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
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr),
         exposure_year0 = exposure_year-355/356, #start is 1 day into the previous year. may be fine to just use 1 for simplicity
         age_start_exposure = cut(age_start_exposure, seq(min(age_start_exposure), max(age_start_exposure), 2), right = F)) %>% 
  run_cox_many_times(., start_time = "exposure_year0", end_time = "exposure_year", 
                     other_predictors = c("age_start_exposure", #"cal_2yr", 
                                          "degree", "nses_z_cx", "male", "race_white")) %>%  
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
    filter(exposure_year < start_of_bad_yrs)  %>%
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
    filter(exposure_year < start_of_bad_yrs)  %>%
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
# SENSITIVITY: USE ST PM2.5 (INSTEAD OF SP) PM2.5 IN MODELS
######################################################################
hrs_st_pm25 <- sur_w %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>% 
  run_cox_many_times(., pollutant_predictors = st_pm2.5_models) %>%
  mutate(
    description = "ST PM2.5"
  )

saveRDS(hrs_st_pm25, file.path(output_data_path, "hazard_ratios_alternative_st_pm2.5.rda"))

######################################################################
# SAVE HAZARD RATIOS
######################################################################
hrs <- rbind(hrs_main, hrs_other_outcomes, hrs_other_exposure_periods, #hrs_extended_cohort, 
             hrs_no_ipw, hrs_baseline_age, hrs_calendar_axis, hrs_2010,
             hrs_dont_drop_last_2yr,
             hrs_effect_modification, hrs_effect_modification_pval)

saveRDS(hrs, file.path(output_data_path, "hazard_ratios.rda"))


######################################################################
# ALTERNATIVE ANALYSIS: ST PM2.5 MODELS
######################################################################
# thes HRs use the full ST PM2.5 dataset, without special considerations for e.g. MM coverage variables

# --> note: need to update code above to include an ST PM2.5 coverage threshold if we publish this

pm2.5_hrs <- sur_w_pm25 %>%
  filter(exposure_year < start_of_bad_yrs)  %>%
  mutate(cal_2yr = droplevels(cal_2yr)) %>%
  run_cox_many_times(., pollutant_predictors = st_pm2.5_models)  %>%
  mutate(#description = ifelse(grepl("pnc_20_36|ufp_36_1k|ns_", pollutant_predictors), "PNC by Size",
  #                             ifelse(grepl("pmdisc_sz", pollutant_predictors), "Particle Size",
  #                                    ifelse(grepl("pnc_onrd", pollutant_predictors), "PNC, Onroad Model", "Main Analysis"
  #                                    )))
  description = "ST PM2.5"
)

saveRDS(pm2.5_hrs, file.path(output_data_path, "hazard_ratios_st_pm2.5.rda"))

######################################################################
message("done with 3_run_models.R")


