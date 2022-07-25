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

source("0_functions.R")
# output_data_path <- file.path("Data", "Output", "1_prep_data")
# if(!file.exists(output_data_path)) {dir.create(output_data_path)}

image_path <- file.path("..", "manuscript", "images")
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


######################################################################
# MM vs ST COVERAGE
######################################################################

# --> WHY DOES MM sometiems have more coverage than ST?

ggplot(data=sur0, aes(x=exp_wks_coverage01_yr_MM, y = exp_wks_coverage01_yr_ST)) + geom_point(alpha=0.1) + geom_abline(slope = 1, intercept = 0, color="yellow") + geom_vline(xintercept = 0.95, color="red")
ggsave(file.path(image_path, "Other", "mm_vs_st_01yr_coverage.png"), width=6, height=6)
prop.table(table(sur0$exp_wks_coverage01_yr_MM > sur0$exp_wks_coverage01_yr_ST))

ggplot(data=sur0, aes(x=exp_wks_coverage10_yr_MM, y = exp_wks_coverage10_yr_ST)) + geom_point(alpha=0.1) + geom_abline(slope = 1, intercept = 0, color="yellow") + geom_vline(xintercept = 0.95, color="red")
ggsave(file.path(image_path, "Other", "mm_vs_st_10yr_coverage.png"), width=6, height=6)
# prop of time where MM > ST. # 35%
prop.table(table(sur0$exp_wks_coverage10_yr_MM > sur0$exp_wks_coverage10_yr_ST))



######################################################################
# CLEAN VARIABLES
######################################################################
sur <- sur0 %>%
  #don't use non-final variables
  select(-c(nindx)) %>%
  filter(#pollutant %in% c("ufp_20_1k","bc", "no2", "pm25"),
          grepl("ufp|bc|no2|pm25", pollutant),
          
         #drop ppl w/ only baseline visits
         last_visit > intakedt #801 ppl dropped# last intakedt is thus 2017-08-01 (vs 2018-09-27)
         ) %>%
  mutate(
    #2 yr cal categories except for years 2016-2018
    # --> update if get years > 2018
    #cal_2yr = factor(ifelse(exposure_year==2018, 2016, floor(exposure_year/2)*2)),
    cal_2yr = factor(floor(exposure_year/2)*2),
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
    ## filter out invalid casi scores - not using this
    #casi_irt = ifelse(casi_valid ==1, casi_irt, NA),
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
sur <- left_join(sur, apoe_wts) %>%
  ungroup()
sur_person <- left_join(sur_person, apoe_wts) %>%
  ungroup() %>%
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

# drop predictions & indicator variables if below the coverage_threshold
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
 

#test %>% select(exp_avg10_yr_MM, exp_wks_coverage10_yr_MM) %>% View() 
  
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

# drop person-years without predictions
sur2 <- sur2 %>% drop_na(pollutant_prediction) 
# --> note, to look at indicator variables, filter to model == "MM". This value repeates for all other models

# final dataset
final_py <- sur2 %>%
  drop_na(pollutant_prediction) %>%
  # take person-years
  group_by(pollutant, exp_avg0) %>%
  summarize(
    final_persons = length(unique(study_id)),
    final_person_years = n()) 

# table of initial & final person years 
left_join(initial_py, final_py) %>%
  mutate(
    prop_persons_dropped = (initial_persons - final_persons)/initial_persons,
    prop_person_years_dropped = (initial_person_years - final_person_years)/initial_person_years
  ) %>%
  kable(caption = "Persons & person-years before and after dropping person-years with insufficient coverage", digits = 3) %>% 
  kable_styling() %>%
  save_kable(file.path(image_path, "Other", "sample_person_yrs.pdf"))

######################################################################
# FINAL INDIVIDUAL LEVEL FOLLOW UP YEARS (WITH EXPOSURE)
######################################################################

complete_fu_yrs <- sur2 %>%
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

######################################################################
# SAVE DATA
######################################################################
# all data
saveRDS(sur2, file.path("Data", "Output", "sur.rda"))
saveRDS(sur, file.path("Data", "Output", "sur_all_data.rda"))
saveRDS(sur_person, file.path("Data", "Output", "sur_person.rda"))
