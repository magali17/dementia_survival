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

pnc_units <- 1000
bc_units <- 10
no2_units <- 5
pm25_units <- 1

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

sur <- drop_na(sur, apoe, predictors)
  
# for multipollutant modeling
sur_w <- sur %>%
  #filter(pollutant %in% c("pm25", "no2")) %>%
  select(-exp_avg0) %>%
  pivot_wider(names_from = c(model, pollutant), values_from = pollutant_prediction) 

######################################################################
# FUNCTIONS
######################################################################

# dt = group_split(sur, pollutant, exposure_duration, model)[[2]]

run_cox <- function(dt, event_indicator, predictors) {
  model_description <- paste(first(dt$pollutant), first(dt$exposure_duration), first(dt$model))
  message(model_description)
  
  #create a survival object
  surv_object <- Surv(time = dt$age_start_exposure,
                      time2 = dt$age_end_exposure,
                      event = dt[[event_indicator]])
  
  cox_model <- coxph(as.formula(paste("surv_object ~" ,paste(predictors, collapse = "+"), "+ strata(apoe)")), 
                     data=dt,
                     cluster = study_id, #robust = T,
                     weights = model_wt) %>%
    summary()
  
  # --> this uses robust SE to calculate 95% conf int, right? 
  result <- cox_model$conf.int %>% 
    as.data.frame() %>%
    rownames_to_column(var = "covariate") %>%
    mutate(pollutant = first(dt$pollutant),
           exposure_duration = first(dt$exposure_duration),
           model = first(dt$model),
           outcome = event_indicator) %>%
    rename_all(~make.names(.)) %>%
    select(-exp..coef.) %>%
    rename(hr = exp.coef.,
           lower95 = lower..95,
           upper95 = upper..95)
  
  return(result)
  }

######################################################################
# SINGLE POLLUTANT MODELS
######################################################################
#predictors other than strata(apoe)
one_pollutant_predictors <- c("pollutant_prediction", "cal_2yr", "male", "race_white", "degree", "income_cat")

one_pollutant_models <- data.frame()
for(i in c("dementia_now", "ad_now")) {
  message(i)
  # x = group_split(sur, pollutant, exposure_duration, model)[[1]]
  temp <- lapply(group_split(sur, pollutant, exposure_duration, model), run_cox, event_indicator = i, predictors=one_pollutant_predictors) %>% 
    bind_rows()
  one_pollutant_models <- rbind(one_pollutant_models, temp)
}

######################################################################
# PM2.5+NO2 MODELS FOR MM, SP, and ST
######################################################################
model_types <- c("MM", "SP", "ST")
pollutants <- c("no2", "pm25")
two_pollutant_models <- data.frame()

for(m in model_types) {
  # m = "MM"
  message(paste("model: ", m))
  model_predictors <- c(paste(m, pollutants, sep = "_"), 
                        "cal_2yr", "male", "race_white", "degree", "income_cat")
  
  for(i in c("dementia_now", "ad_now")) {
    # i = "dementia_now"
    message(paste("outcome: ", i))
    # x = group_split(sur_w, exposure_duration)[[1]]
    temp <- lapply(group_split(sur_w, exposure_duration), function(x) {
      x$model <- m
      x$pollutant <- paste(pollutants, collapse = "+")
      temp <- run_cox(dt=x, event_indicator = i, predictors=model_predictors)}) %>% 
      bind_rows()
      
      two_pollutant_models <- rbind(two_pollutant_models, temp)
    }
  }
 

######################################################################
# PNC+BC+NO2 FOR MM
######################################################################
m <- "MM"
pollutants <- c("ufp_10_42", "bc", "no2")
three_pollutant_mm_models <- data.frame()

model_predictors <- c(paste(m, pollutants, sep = "_"), 
                      "cal_2yr", "male", "race_white", "degree", "income_cat")

for(i in c("dementia_now", "ad_now")) {
  message(paste("outcome: ", i))
  temp <- lapply(group_split(sur_w, exposure_duration), function(x) {
    x$model <- m
    x$pollutant <- paste(pollutants, collapse = "+")
    temp <- run_cox(dt=x, event_indicator = i, predictors=model_predictors)}) %>% 
    bind_rows()
  
  three_pollutant_mm_models <- rbind(three_pollutant_mm_models, temp)
}


######################################################################
# PNC+BC+NO2+PM2.5 FOR MM
######################################################################
pollutants <- c("ufp_10_42", "bc", "no2", "pm25")
four_pollutant_mm_models <- data.frame()

model_predictors <- c(paste(m, pollutants, sep = "_"), 
                      "cal_2yr", "male", "race_white", "degree", "income_cat")

for(i in c("dementia_now", "ad_now")) {
  message(paste("outcome: ", i))
  temp <- lapply(group_split(sur_w, exposure_duration), function(x) {
    x$model <- m
    x$pollutant <- paste(pollutants, collapse = "+")
    temp <- run_cox(dt=x, event_indicator = i, predictors=model_predictors)}) %>% 
    bind_rows()
  
  four_pollutant_mm_models <- rbind(four_pollutant_mm_models, temp)
}


######################################################################
# SAVE MODELS
######################################################################
hazard_ratios <- rbind(one_pollutant_models, two_pollutant_models, three_pollutant_mm_models, four_pollutant_mm_models)
saveRDS(hazard_ratios, file.path(output_data_path, "hazard_ratios.rda"))
