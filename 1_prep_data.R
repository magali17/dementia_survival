######################################################################
# Author: Magali Blanco
# Date: 6/4/2022
# Purpose: Clean and prepare data for a dementia surival analsyis
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
pacman::p_load(tidyverse, lubridate)    

set.seed(1)

source("0_functions.R")

######################################################################
# LOAD DATA
######################################################################
folder_date <- "20220610"
dt_path <- file.path("Data", "Raw", folder_date, "issue_001.rda")

if(file.exists(dt_path)) {
  sur0 <- readRDS(dt_path)
  } else {
     sur0 <- haven::read_sas(gsub(".rda", ".sas7bdat", dt_path) , NULL)
     saveRDS(sur0, dt_path)
   }

#remove SAS labels
#labelled::var_label(sur0) <- NULL

######################################################################
# CLEAN VARIABLES
######################################################################
sur <- sur0 %>%
  #don't use non-final variables
  select(-c(nindx)) %>%
  # --> ADD ST 2019
  filter(pollutant %in% c("ufp_20_1k","bc", "no2", "pm25"),
         #drop ppl w/ only baseline visits
         last_visit > intakedt #801 ppl dropped
         ) %>%
  mutate(
    #2 yr cal categories except for years 2016-2018
    cal_2yr = factor(ifelse(exposure_year==2018, 2016, floor(exposure_year/2)*2)),
    # Unless someone was suspected of X and was evaluated for it, X will be blank. ‘0’ is a confirmation 'no' while blank is a presumed 'no'
    corrected_anydementia = ifelse(is.na(corrected_anydementia), 0, corrected_anydementia),
    corrected_dsmivdx = ifelse(is.na(corrected_dsmivdx), 0, corrected_dsmivdx),
    final_nindx = ifelse(is.na(final_nindx), 0, final_nindx),
    #Time-varying covariates. Could also use age_at_exposure, although it will be less accurate
    ##starting age is age on Jan 1st of a given year if participant was enrolled, or on intake date, whichever came later
    age_start_exposure = ifelse(year(intakedt) == exposure_year, 
                                as.numeric(intakedt-birthdt)/365,
                                as.numeric(as.Date(paste0(exposure_year, "-01-01")) - birthdt)/365), 
    ##ending age is age at end of year or last visit date, whichever came first 
    age_end_exposure = ifelse(format(last_visit, "%Y") == exposure_year,
                              as.numeric(last_visit - birthdt)/365,
                              as.numeric(as.Date(paste0(exposure_year, "-12-31")) - birthdt)/365),
    ### need to change age_last_visit to age_act for sensitivity analyses
    #dementia_now = ifelse(age_end_exposure < age_last_visit, 0, anydementia),
    ad_nincds = ifelse(final_nindx %in% c(1,2), 1, 0),
    #ad_now = ifelse(age_end_exposure < age_last_visit, 0, ad_nincds),
    ## filter out invalid casi scores 
    casi_irt = ifelse(casi_valid ==1, casi_irt, NA),
    income_cat = ifelse(tr_med_inc_hshld < 35000, 1,  
                        ifelse(tr_med_inc_hshld >= 35000 & tr_med_inc_hshld < 50000, 2,
                               ifelse(tr_med_inc_hshld >= 50000 & tr_med_inc_hshld < 75000, 3, 4))),  
    race_white = ifelse(race == 1, 1, 0),
    # --> should 9 be a 6 instead? but would change Rachel's models
    # make unknown degree=9 "none" 
    degree = ifelse(degree == 9, 0, degree),
    # combine GED and HS 
    degree = ifelse(degree %in% c(1:2), 1, degree)
    ) %>%
  #time-varying covariates
  group_by(study_id) %>%
  # --> need to change this if we use act onset age b/c it is 1 yr before the last visit
  mutate(
    dementia_now = ifelse(exposure_year == max(exposure_year), corrected_anydementia, 0),
    ad_now = ifelse(exposure_year == max(exposure_year), ad_nincds, 0),
  ) %>%
  ungroup() %>%
  add_factor_refs()

# --> why is income coded as 1-6,9 AND A-F (6)? miscoded?

# --> need to change this if use Onset_Age_ACT (1 yr less for cases)? or will only the age change
# --> if use act onset age, drop last row for cases

# individual level dataset
fixed_covars <- names(sur)[!grepl("exp_avg|exp_wks|coverage|imp_qual|cal_|exposure|pollutant|dementia_now|ad_now|onsetage|onsetdate", names(sur))] #DIAGCRIT
sur_person <- select(sur, all_of(fixed_covars)) %>% distinct()
 
######################################################################
# MISSING PERSON LEVEL COVARIATES
######################################################################
missing <- sur_person %>%
  summarize_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "covariate", values_to = "count") %>%
  mutate(prop_missing = count/nrow(sur_person)) %>%
  filter(prop_missing >0) %>%
  arrange(-prop_missing) 

saveRDS(missing, file.path("Data", "Output", "missing_table.rda"))

# drop columns with too much missingness, defined as more than apoe (12%)
apoe_missing <- missing$prop_missing[missing$covariate=="apoe"] #0.1220489
too_much_missingness <- missing %>%
  filter(prop_missing>apoe_missing) %>%
  pull(covariate)

sur <- select(sur, -too_much_missingness)
sur_person <- select(sur_person, -too_much_missingness)

# update missing covariates that we want to fill
#missing_vars <- missing$covariate [missing$covariate != "apoe"]
missing_vars <- setdiff(setdiff(missing$covariate, too_much_missingness), "apoe")
#summary(missing$prop_missing[missing$covariate %in% missing_vars]) # mean (range): 0.0097453 (0.0002108, 0.0944351)

impute_value <- function(dt) {
  # remove factor status
  dt <- as.numeric(as.character(dt))
  is_integer <- all(round(dt) == dt, na.rm = T)
  
  # categorical variables are coded as integers. round to make sure these continue to be integers (~mode)
  if(is_integer) {mean <- round(mean(dt, na.rm=T))}  
  if(!is_integer) {mean <- mean(dt, na.rm=T)}  #casi_irt, grip_strength
  
  dt <- ifelse(is.na(dt), mean, dt)
  return(dt)
  }

# manually impute occupation since values are not ordinal. only 5 people
occupation_mode <- as.data.frame(table(sur_person$occupation)) %>% 
  filter(Freq==max(Freq)) %>% pull(Freq)
sur_person <- mutate(sur_person, occupation = ifelse(is.na(occupation), occupation_mode, occupation))

missing_vars_imputed <- apply(sur_person[missing_vars], 2, impute_value) %>%
  as.data.frame()
 
sur_person <- cbind(select(sur_person, -missing_vars), missing_vars_imputed )
sur <- left_join(select(sur, -missing_vars), cbind(study_id=sur_person$study_id, missing_vars_imputed) )

# variables w/o missingness (for IPW later)
still_missing <- sur_person %>%
  summarize_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "covariate", values_to = "count") %>%
  mutate(prop_missing = count/nrow(sur_person)) %>%
  #filter(prop_missing >0) %>%
  arrange(-prop_missing) 

######################################################################
# IPW FOR APOE MISSINGNESS
######################################################################
apoe <- mutate(sur_person, apoe_available = !is.na(apoe)) %>%
  select(-apoe) %>%
  drop_na() %>%
  mutate_at(c("birth_cohort", "birth_cohort_2yr", "final_nindx", "cohort", "casi_valid", "marital", "livingsb", "iadl_flag", "iadl_sum", "IADL_Bills", "IADL_House", "IADL_Meals", "IADL_phone", "IADL_shopping", "cholmed", "Chair_Able", "Gait_able", "grip_able", "gripsmall", "halfmile", "heavyhouse", "lift10lbs", "stairs", "degree", "race", "htnmed"), as.factor)

# use all columns except study_id and apoe_available as predictors
lasso_results <- lasso_fn(dt = apoe, x_names = names(apoe)[2:(ncol(apoe)-1)], y_name = "apoe_available", family. = "binomial", lambda. = 0.01) #lambda 0.005

apoe_missing_cov <- lasso_results$results$cov

# model.matrix creates dummy variables like those selected from lasso
apoe <- model.matrix(apoe_available~., data=apoe) %>% as.data.frame() %>%
     cbind(apoe_available=apoe$apoe_available)

#don't include subcategories w/ few counts (produces very unstable weights; e.g., probability of apoe being available is ~0)
drop_vars <- apoe %>%
  group_by(apoe_available) %>%
  select(apoe_missing_cov) %>%
  summarize_all(~sum(.)) %>%  
  select(-apoe_available) %>%
  apply(.,2, function(x) any(x<=1)) %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  filter(.==TRUE) %>%
  pull(rowname)

apoe_missing_cov <- setdiff(apoe_missing_cov, drop_vars)

apoe_model <- glm(as.formula(paste("apoe_available ~", paste(apoe_missing_cov, collapse = "+"))), family = "binomial", data=apoe)
# denominator: probability of APOE being present
apoe_available_prob <- predict(apoe_model, type = "response") 
# summary(1/apoe_available_prob) #check that probability isn't Tiny (like when used "race")

#numerator: stabilizer 
apoe_male_m <- glm(apoe_available ~ male1, family = "binomial", data = apoe) 
apoe_male_prob <- predict(apoe_male_m, type = "response") 
# summary(apoe_male_prob/apoe_available_prob)

apoe_wts <- select(apoe, study_id) %>%
  mutate(model_wt = apoe_male_prob/apoe_available_prob)
# summary(apoe_wts$model_wt)

# add to datasets 
sur <- left_join(sur, apoe_wts)
sur_person <- left_join(sur_person, apoe_wts)

######################################################################
# RECODE FACTORS FOR MODELING
######################################################################
sur <- add_factor_refs(sur)

# don't need this for sur_person b/c this is just for descriptives 

######################################################################
# SAVE DATA
######################################################################
saveRDS(sur, file.path("Data", "Output", "sur.rda"))
saveRDS(sur_person, file.path("Data", "Output", "sur_person.rda"))
