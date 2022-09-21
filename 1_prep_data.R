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
pacman::p_load(tidyverse, lubridate, kableExtra)    

set.seed(1)

source("functions.R")
output_data_path <- file.path("Data", "Output", "1_prep_data")
if(!file.exists(output_data_path)) {dir.create(output_data_path)}

image_path <- file.path("..", "manuscript", "images")
first_exposure_year <- 2005
save(first_exposure_year, file = file.path(output_data_path, "first_exposure_year.rda"))

print_diagnostics <- FALSE

# exclusion counts
exclusion_table <- tibble(description = as.character(),  
                        persons = as.numeric(),
                        persons_dropped = as.numeric(),
                        persons_pct_dropped = as.numeric(),
                        person_years = as.numeric(),
                        person_years_dropped = as.numeric(),
                        person_years_pct_dropped = as.numeric()
                        )

# fn adds the number of unique persons and person-years for any given dataset to the exclusion tibble
# dt=sur0
# description. = "test"
count_remaining_sample <- function(dt, description.) {
  temp <- distinct(dt, study_id, exposure_year)
  
  exclusion_table <- add_row(exclusion_table,
                       description = description.,
                       persons = length(unique(temp$study_id)),
                       person_years = nrow(temp))
  
  exclusion_table <- mutate(exclusion_table,
                            persons_dropped = lag(persons) - persons,
                            persons_pct_dropped = round(persons_dropped/lag(persons), 2)*100,
                            person_years_dropped = lag(person_years) - person_years,
                            person_years_pct_dropped = round(person_years_dropped/lag(person_years), 2)*100
                            )
  
  return(exclusion_table)
}
######################################################################
# LOAD DATA
######################################################################
folder_date <- "20220623" #"20220610"
dt_path <- file.path("Data", "Raw", folder_date, "issue_001.rda")

if(file.exists(dt_path)) {
  sur0 <- readRDS(dt_path)
  } else {
     sur0 <- haven::read_sas(gsub(".rda", ".sas7bdat", dt_path) , NULL)
     saveRDS(sur0, dt_path)
   }

exclusion_table <- count_remaining_sample(sur0, description. = "Full cohort/raw data")

######################################################################
# MM vs ST COVERAGE
######################################################################

# --> WHY DOES MM sometimes have more coverage than ST?

if(print_diagnostics == TRUE) {
  ggplot(data=sur0, aes(x=exp_wks_coverage01_yr_MM, y = exp_wks_coverage01_yr_ST)) + geom_point(alpha=0.1) + geom_abline(slope = 1, intercept = 0, color="yellow") + geom_vline(xintercept = 0.95, color="red")
  ggsave(file.path(image_path, "Other", "mm_vs_st_01yr_coverage.png"), width=6, height=6)
  prop.table(table(sur0$exp_wks_coverage01_yr_MM > sur0$exp_wks_coverage01_yr_ST))
  
  ggplot(data=sur0, aes(x=exp_wks_coverage10_yr_MM, y = exp_wks_coverage10_yr_ST)) + geom_point(alpha=0.1) + geom_abline(slope = 1, intercept = 0, color="yellow") + geom_vline(xintercept = 0.95, color="red")
  ggsave(file.path(image_path, "Other", "mm_vs_st_10yr_coverage.png"), width=6, height=6)
  # prop of time where MM > ST. # 35%
  prop.table(table(sur0$exp_wks_coverage10_yr_MM > sur0$exp_wks_coverage10_yr_ST))
}

######################################################################
# CLEAN VARIABLES
######################################################################
sur <- sur0 %>%
  #don't use non-final variables
  select(-c(nindx)) %>%
  filter(grepl("ufp|bc|no2|pm25", pollutant),
         #drop ppl w/ only baseline visits
         last_visit > intakedt #801 ppl dropped# last intakedt is thus 2017-08-01 (vs 2018-09-27)
         ) %>%
  mutate(
    #2 yr cal categories except for years 2016-2018
    cal_2yr = factor(floor(exposure_year/2)*2),
    # Unless someone was suspected of X and was evaluated for it, X will be blank. ‘0’ is a confirmation 'no' while blank is a presumed 'no'
    corrected_anydementia = ifelse(is.na(corrected_anydementia), 0, corrected_anydementia),
    corrected_dsmivdx = ifelse(is.na(corrected_dsmivdx), 0, corrected_dsmivdx),
    vad_dementia = ifelse(corrected_dsmivdx==2, 1, 0),
    mixed_dementia = ifelse(corrected_dsmivdx==5, 1, 0),
    other_dementia = ifelse(corrected_dsmivdx %in% c(3,4,6), 1, 0),
    non_ad_dementia = ifelse(corrected_dsmivdx %in% c(2:6), 1, 0),
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
    ### may need to change age_last_visit to age_act for sensitivity analyses
    #dementia_now = ifelse(age_end_exposure < age_last_visit, 0, anydementia),
    ad_nincds = ifelse(final_nindx %in% c(1,2), 1, 0),
    #ad_now = ifelse(age_end_exposure < age_last_visit, 0, ad_nincds),
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
  # would need to change this if we use act onset age b/c it is 1 yr before the last visit
  mutate(
    dementia_now = ifelse(exposure_year == max(exposure_year), corrected_anydementia, 0),
    ad_now = ifelse(exposure_year == max(exposure_year), ad_nincds, 0),
    vad_now = ifelse(exposure_year == max(exposure_year), vad_dementia, 0),
    mixed_now = ifelse(exposure_year == max(exposure_year), mixed_dementia, 0),
    other_now = ifelse(exposure_year == max(exposure_year), other_dementia, 0),
    non_ad_now = ifelse(exposure_year == max(exposure_year), non_ad_dementia, 0),
  ) %>%
  ungroup() %>%
  add_factor_refs()
# update remaining data counts
exclusion_table <- count_remaining_sample(sur, description. = "people with 1+ follow-up visits")

# save semi-raw dataset with group categories used in analysis later (e.g., race, ad diagnosis)
saveRDS(sur, file.path("Data", "Output", "sur0_redefined_categories.rda"))

# --> need to change this if use Onset_Age_ACT (1 yr less for cases)? or will only the age change
# --> if use act onset age, drop last row for cases

# individual level dataset
fixed_covars <- names(sur)[!grepl("exp_avg|exp_wks|coverage|imp_qual|cal_|exposure|pollutant|_now$|onsetage|onsetdate", names(sur))] #DIAGCRIT
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

S(missing, file.path("Data", "Output", "missing_table.rda"))

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

# people with missing
if(print_diagnostics==TRUE) {
  ggplot(apoe, aes(x=last_visit, fill=apoe_available)) + geom_histogram()  
}

# --> is the goal with IPW to have a representative sample in terms of things that could be confounders? or just in "general" - to say that our study cohort is "representative" of the "general" population? This may determine what predictors we choose to use in IPW

# see what categories have few counts, stratified by apoe availability
if(print_diagnostics == TRUE) {
  library(table1)

  # table comparing counts/distribution of covariates in people w/ and w/o apoe values
  apoe %>%
    select(Onset_Age_ACT:last_col()) %>%
    mutate(apoe_available = as.factor(apoe_available)) %>%
    table1(~. | apoe_available, data=., overall = F, extra.col = list('P-value' = pvalue))
}

# use all columns except study_id and apoe_available as predictors
# fn uses model.matrix and then selects overall predictors
lasso_results <- lasso_fn(dt = apoe, x_names = names(apoe)[2:(ncol(apoe)-1)], y_name = "apoe_available", family. = "binomial")#, lambda. = 0.01)
                           

apoe_missing_cov <- lasso_results$results$cov

# # model.matrix creates dummy variables like those selected from lasso
# apoe <- model.matrix(apoe_available~., data=apoe) %>% as.data.frame() %>%
#      cbind(apoe_available=apoe$apoe_available)

# only include predictors where >= X individuals have that characteristic. Lower values indicate that this is rare in the cohort anyways, and leads to model instability while calculating weights
# e.g., 1-2 people out of ~5k have characteristic X and one is dropped b/c of no APOE status.
count_threshold <- 50 #if >~1% of the population has trait Y, use in IPW  #3

apoe_matrix <- apoe %>%
  model.matrix(apoe_available~., data=.) %>% as.data.frame() %>%
  cbind(apoe_available=apoe$apoe_available) %>%
  # don't need or already have (e.g., birth cohort, last_predx_visit is similar to last_visit)
  select(-c("(Intercept)", "study_id", "birthdt", "intakedt", "last_predx_visit")) #"last_visit",

keep_vars <- apoe_matrix %>%
  group_by(apoe_available) %>%
  summarize_all(~sum(.)) %>%
  summarize_all(~all(.>=count_threshold)) %>%
  select_if(~.==TRUE) %>% 
  names() %>%
  c(., "apoe_available")

# make sure variables have enough counts and were selected by Lasso
apoe_missing_cov <- keep_vars[keep_vars %in% apoe_missing_cov]
saveRDS(apoe_missing_cov, file.path(output_data_path, "lasso_predictors.rda"))

apoe_model <- glm(as.formula(paste("apoe_available ~", paste(apoe_missing_cov, collapse = "+"))), family = "binomial", data=apoe_matrix)
# denominator: probability of APOE being present
apoe_available_prob <- predict(apoe_model, type = "response") 
# summary(1/apoe_available_prob) #check that probability isn't Tiny (like when used "race")

#numerator: stabilizer 
apoe_male_m <- glm(apoe_available ~ male1, family = "binomial", data = apoe_matrix) 
apoe_male_prob <- predict(apoe_male_m, type = "response") 
# summary(apoe_male_prob/apoe_available_prob)

apoe_wts <- select(apoe, study_id) %>%
  mutate(model_wt = apoe_male_prob/apoe_available_prob)
# summary(apoe_wts$model_wt)
# table(apoe_wts$model_wt>3) #ppl w/ high weights

# add to datasets 
sur <- left_join(sur, apoe_wts) %>% ungroup()
sur_person <- left_join(sur_person, apoe_wts) %>% ungroup() %>%
  # drop model weights that won't actually be used in the models b/c apoe is missing
  mutate(model_wt = ifelse(is.na(apoe), NA, model_wt))

######################################################################
# RECODE FACTORS FOR MODELING
######################################################################
sur <- add_factor_refs(sur)
# don't need this for sur_person b/c this is just for descriptives 

######################################################################
# COVERAGE VARIABLE
######################################################################
coverage_threshold <- 0.95

# drop person-years if below the coverage_threshold
# since ST predictions will only be used as sensitivity analyses to see if our findings are similar, we want to make sure the person-years used are the same as MM
# thus, we will drop ST predictions if MM predictions are also dropped due to low covarege in the MM region
sur <- sur %>%
  # modify MM, ST, SP predictions based on MM coverage since these are sensitivity analyses
  mutate_at(vars(starts_with("exp_avg10_yr_")), ~ifelse(exp_wks_coverage10_yr_MM < coverage_threshold, NA, .)) %>%
  mutate_at(vars(starts_with("exp_avg05_yr_")), ~ifelse(exp_wks_coverage05_yr_MM < coverage_threshold, NA, .)) %>%
  mutate_at(vars(starts_with("exp_avg01_yr_")), ~ifelse(exp_wks_coverage01_yr_MM < coverage_threshold, NA, .)) %>%
  # make all coverage values NA if there are NA predictions, using MM since this is the same for ST and SP
  mutate_at(vars(starts_with("exp_wks_coverage10_yr_")), ~ifelse(is.na(exp_avg10_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("exp_wks_coverage05_yr_")), ~ifelse(is.na(exp_avg05_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("exp_wks_coverage01_yr_")), ~ifelse(is.na(exp_avg01_yr_MM), NA, .)) %>%
  
  mutate_at(vars(starts_with("exact_coverage10_yr_")), ~ifelse(is.na(exp_avg10_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("exact_coverage05_yr_")), ~ifelse(is.na(exp_avg05_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("exact_coverage01_yr_")), ~ifelse(is.na(exp_avg01_yr_MM), NA, .)) %>%
  
  mutate_at(vars(starts_with("imputed_coverage10_yr_")), ~ifelse(is.na(exp_avg10_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("imputed_coverage05_yr_")), ~ifelse(is.na(exp_avg05_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("imputed_coverage01_yr_")), ~ifelse(is.na(exp_avg01_yr_MM), NA, .)) %>%
  
  mutate_at(vars(starts_with("imp_qual10_yr_")), ~ifelse(is.na(exp_avg10_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("imp_qual05_yr_")), ~ifelse(is.na(exp_avg05_yr_MM), NA, .)) %>%
  mutate_at(vars(starts_with("imp_qual01_yr_")), ~ifelse(is.na(exp_avg01_yr_MM), NA, .)) 
 
# update remaining data counts. this needs to be done below after reformatting, otherwise the same "counts" are reported
#exclusion_table <- count_remaining_sample(sur, description. = "drop person-years if below the coverage threshold")

 
  
# # drop prediction if lower than a threshold
# sur$exp_avg10_yr_MM[sur$exp_wks_coverage10_yr_MM < coverage_threshold] <- NA
# sur$exp_avg10_yr_ST[sur$exp_wks_coverage10_yr_ST < coverage_threshold] <- NA
# sur$exp_avg10_yr_SP[sur$exp_wks_coverage10_yr_SP < coverage_threshold] <- NA
# sur$exp_avg05_yr_MM[sur$exp_wks_coverage05_yr_MM < coverage_threshold] <- NA
# sur$exp_avg05_yr_ST[sur$exp_wks_coverage05_yr_ST < coverage_threshold] <- NA
# sur$exp_avg05_yr_SP[sur$exp_wks_coverage05_yr_SP < coverage_threshold] <- NA
# sur$exp_avg01_yr_MM[sur$exp_wks_coverage01_yr_MM < coverage_threshold] <- NA
# sur$exp_avg01_yr_ST[sur$exp_wks_coverage01_yr_ST < coverage_threshold] <- NA
# sur$exp_avg01_yr_SP[sur$exp_wks_coverage01_yr_SP < coverage_threshold] <- NA
# #recalculate indicator variables for 10 yr MM
# sur$exp_wks_coverage10_yr_MM[is.na(sur$exp_avg10_yr_MM)] <- NA
# sur$exact_coverage10_yr_MM[is.na(sur$exp_avg10_yr_MM)] <- NA
# sur$imputed_coverage10_yr_MM[is.na(sur$exp_avg10_yr_MM)] <- NA
# sur$imp_qual10_yr_MM[is.na(sur$exp_avg10_yr_MM)] <- NA

######################################################################
# SUBSET DATA
######################################################################
sur2 <- sur %>% 
  select(study_id, birth_cohort, intakeage, last_visit_age, onsetage, Onset_Age_ACT, corrected_anydementia, corrected_dsmivdx, dementia_now, ad_now, age_at_exposure, age_start_exposure, age_end_exposure, exposure_year, 
         pollutant, starts_with("exp_avg10_yr_"), starts_with("exp_avg05_yr_"), starts_with("exp_avg01_yr_"),
         exp_wks_coverage10_yr_MM, exact_coverage10_yr_MM, imputed_coverage10_yr_MM, imp_qual10_yr_MM,
    male, hispanic, race, race_white, apoe, bmi4, degree, income_cat, cal_2yr, model_wt
    ) %>%
  pivot_longer(starts_with("exp_avg"), names_to = "exp_avg0", values_to = "pollutant_prediction") %>%
  mutate(
    model = substr(exp_avg0, nchar(exp_avg0)-1, nchar(exp_avg0)),
    exposure_duration = ifelse(grepl("10_yr", exp_avg0), 10, ifelse(grepl("05_yr", exp_avg0), 05, ifelse(grepl("01_yr", exp_avg0), 1, NA)))
  ) %>%
  select(study_id, age_start_exposure, age_start_exposure, age_end_exposure, exposure_year, pollutant, exp_avg0, exposure_duration, model, pollutant_prediction, contains("coverage|imp_"), everything())

# full dataset that still has rows w/ missing/unused exposure values
initial_py <- sur2 %>%
  #filter(grepl("10_yr", exp_avg0)) %>%
  group_by(pollutant, exp_avg0) %>%
   
  summarize(
    initial_persons = length(unique(study_id[!all(is.na(pollutant_prediction))])),
    initial_person_years = length(pollutant_prediction[!all(is.na(pollutant_prediction))])
  ) %>% 
 filter(initial_persons>0)

# drop person-years without predictions after reformatting 
sur2 <- sur2 %>% drop_na(pollutant_prediction) 
#remaining data
exclusion_table <- sur2 %>%
  filter(model=="MM", exposure_duration==10) %>%
  count_remaining_sample(description. = "person-years above the coverage threshold, for 10 yr MM; after formatting data")

saveRDS(sur2, file.path("Data", "Output", "sur_full_cohort.rda"))

######################################################################
# REDUCED ANALYSIS COHORT 
######################################################################
# people with APOE
sur3 <- filter(sur2, !is.na(apoe))
# remaining data
exclusion_table <- sur3 %>%
  filter(model=="MM", exposure_duration==10) %>%
  count_remaining_sample(description. = "cohort with APOE values")

# 2005+ cohort 
sur3 <- filter(sur3, exposure_year >= first_exposure_year)
# remaining data
exclusion_table <- sur3 %>%
  filter(model=="MM", exposure_duration==10) %>%
  count_remaining_sample(description. = "2005+ Cohort")

## counts for 2010+
exclusion_table <- filter(sur3, exposure_year >= 2010) %>%
  filter(model=="MM", exposure_duration==10) %>%
  count_remaining_sample(description. = "2010+ Cohort")

######################################################################
# FINAL INDIVIDUAL LEVEL FOLLOW UP YEARS (WITH EXPOSURE)
######################################################################

#complete_fu_yrs <- sur2 %>%
complete_fu_yrs <- sur3 %>%
  filter(pollutant == first(pollutant),
         model == first(model),
         exposure_duration == first(exposure_duration)
  ) %>% 
  group_by(study_id) %>%
  mutate(follow_up_years = n()) %>% 
  distinct(study_id, follow_up_years)  %>%
  ungroup()

# left join to complete_fu_yrs to only keep people remaining
sur_person <- left_join(complete_fu_yrs, sur_person)

# summary(apoe_wts$model_wt)

######################################################################
# SAVE DATA
######################################################################
# all data
saveRDS(sur3, file.path("Data", "Output", "sur.rda"))
#saveRDS(sur, file.path("Data", "Output", "sur_all_data.rda"))
saveRDS(sur_person, file.path("Data", "Output", "sur_person.rda"))
saveRDS(exclusion_table, file.path("Data", "Output", "exclusion_table.rda"))

