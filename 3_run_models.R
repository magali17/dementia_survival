# NOTE: We will probably want to replicate much of what we have published, for completeness.  But unless we get different findings, I think that will be more of an exercise to check results rather than something we will focus on.

######################################################################
# Author: Magali Blanco
# Date: 7/22/2022
# Purpose: Run dementia-TRAP survival model
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

source("0_functions.R")
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

predictors <- c("pollutant_prediction", "cal_2yr", "male", "race_white", "degree", "income_cat")

pnc_units <- 1000
bc_units <- 10
no2_units <- 5
pm25_units <- 1

sur <- sur %>%
  group_by(pollutant) %>%
  mutate(
    pollutant_prediction2 = case_when(
      grepl("ufp", pollutant) ~ pollutant_prediction/pnc_units,
      pollutant == "bc" ~pollutant_prediction/bc_units,
      pollutant == "no2" ~pollutant_prediction/no2_units,
      pollutant == "pm25" ~pollutant_prediction/pm25_units
    )) %>%
  ungroup()

sur <- drop_na(sur, apoe, predictors)
   
######################################################################
# FUNCTIONS
######################################################################

# dt = group_split(sur, pollutant, exposure_duration, model)[[2]]

run_cox <- function(dt, event_indicator = "dementia_now") {
  model_description <- paste(first(dt$pollutant), first(dt$exposure_duration), first(dt$model))
  message(model_description)
  
  # print("data rows")
  # print(nrow(dt))
  
  #create a survival object
  surv_object <- Surv(time = dt$age_start_exposure,
                      time2 = dt$age_end_exposure,
                      event = dt[[event_indicator]])
  
  # print("surv_object length")
  # print(length(surv_object))
  
  #print("fitting model")
  cox_model <- coxph(as.formula(paste("surv_object ~" ,paste(predictors, collapse = "+"), "+ strata(apoe)")), 
                     data=dt,
                     cluster = study_id, #robust = T,
                     weights = model_wt) %>%
    summary()
  
  #print("summarizing results")
  # --> this uses robust SE to calculate 95% conf int, right? 
  result <- cox_model$conf.int %>% 
    as.data.frame() %>%
    rownames_to_column(var = "covariate") %>%
    mutate(pollutant = first(dt$pollutant),
           exposure_duration = first(dt$exposure_duration),
           model = first(dt$model),
           outcome = event_indicator
           ) %>%
    rename_all(~make.names(.)) %>%
    select(-exp..coef.) %>%
    rename(hr = exp.coef.,
           lower95 = lower..95,
           upper95 = upper..95,
           )
  
  return(result)
  }

######################################################################
# SINGLE POLLUTANT MODELS
######################################################################
# x = group_split(sur, pollutant, exposure_duration, model)[[1]]
one_pollutant <- data.frame()
for(i in c("dementia_now", "ad_now")) {
  message(i)
  
  temp <- lapply(group_split(sur, pollutant, exposure_duration, model), run_cox, event_indicator = i) %>%
     bind_rows()
   one_pollutant <- rbind(one_pollutant, temp)
}
# dementia_results <- lapply(group_split(sur, pollutant, exposure_duration, model), run_cox, event_indicator = "dementia_now") %>%
#   bind_rows()
# ad_results <- lapply(group_split(sur, pollutant, exposure_duration, model), run_cox, event_indicator = "ad_now") %>%
#   bind_rows()


######################################################################
# MULTI-POLLUTANT MODELS
######################################################################





