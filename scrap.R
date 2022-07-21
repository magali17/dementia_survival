##########################################################################################
# NOAH TABLE 1
##########################################################################################









##########################################################################################
# OLD TABLE 1
##########################################################################################

#function returns Table 1 summary statistics for a given dataset, with a column name describing group 

library(qwraps2)
# define the markup language we are working in.
options(qwraps2_markup = "markdown")

t1.fn <- function(data, 
                  #description variables
                  column.name = "") {  
  
  #pacman::p_load(kableExtra)
  
  round.var <- 0
  
  
  t1 <- data %>%
    dplyr::summarize(
      N = n(),
      N_perc = round(N/nrow(.)*100, round.var),
      
      pys_total = sum(fu_yrs), #qwraps2::n_perc(fu_yrs, digits = 0, show_denom = "never", na_rm = TRUE), #tot_py[1], 
      
      follow_up_yrs_mean = round(mean(fu_yrs), round.var), 
      follow_up_yrs_sd = round(sd(fu_yrs), round.var),
      
      birth_cohort1_n = sum(birth_cohort==1895, na.rm=T),
      birth_cohort1_pct = round(birth_cohort1_n/N*100, round.var),
      birth_cohort2_n = sum(birth_cohort==1910, na.rm=T),
      birth_cohort2_pct = round(birth_cohort2_n/N*100, round.var),
      birth_cohort3_n = sum(birth_cohort==1915, na.rm=T),
      birth_cohort3_pct = round(birth_cohort3_n/N*100, round.var),
      birth_cohort4_n = sum(birth_cohort==1920, na.rm=T),
      birth_cohort4_pct = round(birth_cohort4_n/N*100, round.var),
      birth_cohort5_n = sum(birth_cohort==1925, na.rm=T),
      birth_cohort5_pct = round(birth_cohort5_n/N*100, round.var),
      birth_cohort6_n = sum(birth_cohort==1930, na.rm=T),
      birth_cohort6_pct = round(birth_cohort6_n/N*100, round.var),
      birth_cohort7_n = sum(birth_cohort==1935, na.rm=T),
      birth_cohort7_pct = round(birth_cohort7_n/N*100, round.var),
      
      exposure = mean_sd(x = no2_10yr, digits = 0, show_n = "never", denote_sd = "paren", na_rm = T),
      
      dementia_cases = n_perc(anydementia==1, digits = 0),
      ad_cases = n_perc(ad_nincds==1, digits = 0),
      
      age_entry_median = round(median(age_intake), round.var),
      age_entry_iqr = round(IQR(age_intake), round.var),
      female_n = sum(male == "0"), #sum(!male),
      female_pct = round(female_n/N*100, round.var),
      race_known_n = sum(!is.na(race_white)),
      race_white_n = sum(race_white =="1", na.rm = T),
      race_white_pct = round(race_white_n/race_known_n*100, round.var),
      edu_known_n = sum(!is.na(degree)),
      degree_known_n = sum(!is.na(degree)),
      degree0_n = sum(degree == "0", na.rm = T),
      degree0_pct = round(degree0_n/degree_known_n*100, round.var),
      #degree==1 is for GED or HS
      degree1_n = sum(degree == "1", na.rm = T),
      degree1_pct = round(degree1_n/degree_known_n*100, round.var),
      # degree2_n = sum(degree ==2, na.rm = T),
      # degree2_pct = round(degree2_n/degree_known_n*100, 1),
      degree3_n = sum(degree =="3", na.rm = T),
      degree3_pct = round(degree3_n/degree_known_n*100, round.var),
      degree4_n = sum(degree =="4", na.rm = T),
      degree4_pct = round(degree4_n/degree_known_n*100, round.var),
      degree5_n = sum(degree =="5", na.rm = T),
      degree5_pct = round(degree5_n/degree_known_n*100, round.var),
      degree6_n = sum(degree =="6", na.rm = T),
      degree6_pct = round(degree6_n/degree_known_n*100, round.var),
      income_cat_known_n = sum(!is.na(income_cat)),
      income1_n = sum(income_cat=="1", na.rm=T),
      income1_pct = round(income1_n/income_cat_known_n*100, round.var),
      income2_n = sum(income_cat=="2", na.rm=T),
      income2_pct = round(income2_n/income_cat_known_n*100, round.var),
      income3_n = sum(income_cat=="3", na.rm=T),
      income3_pct = round(income3_n/income_cat_known_n*100, round.var),
      income4_n = sum(income_cat=="4", na.rm=T),
      income4_pct = round(income4_n/income_cat_known_n*100, round.var),
      apoe_known_n = sum(!is.na(apoe)),
      apoe_carrier_n = sum(apoe, na.rm=T),
      apoe_carrier_pct = round(apoe_carrier_n/apoe_known_n*100, round.var),
      smoke_known_n = sum(!is.na(smoke)),
      smoke_never_n = sum(smoke=="0", na.rm = T),
      smoke_never_pct = round(smoke_never_n/smoke_known_n*100, round.var),
      smoke_former_n = sum(smoke=="1", na.rm = T),
      smoke_former_pct = round(smoke_former_n/smoke_known_n*100, round.var),
      smoke_current_n = sum(smoke=="2", na.rm = T),
      smoke_current_pct = round(smoke_current_n/smoke_known_n*100, round.var),
      exercise_regular_known_n = sum(!is.na(exercise_reg)),
      exercise_regular_n = sum(exercise_reg=="1", na.rm = T),
      exercise_regular_pct = round(exercise_regular_n/exercise_regular_known_n*100, round.var),
      bmi_known_n = sum(!is.na(bmi)),
      bmi_under_n = sum(bmi=="0", na.rm = T),
      bmi_under_pct = round(bmi_under_n/bmi_known_n*100, round.var),
      bmi_normal_n = sum(bmi=="1", na.rm = T),
      bmi_normal_pct = round(bmi_normal_n/bmi_known_n*100, round.var),
      bmi_over_n = sum(bmi=="2", na.rm = T),
      bmi_over_pct = round(bmi_over_n/bmi_known_n*100, round.var),
      bmi_obese_n = sum(bmi=="3", na.rm = T),
      bmi_obese_pct = round(bmi_obese_n/bmi_known_n*100, round.var),
      hypertension_known_n = sum(!is.na(Hypertension)),
      hypertension_n = sum(Hypertension == "1", na.rm=T),
      hypertension_pct = round(hypertension_n/hypertension_known_n*100, round.var),
      diabetes_known_n = sum(!is.na(Diabetes)),
      diabetes_n = sum(Diabetes == "1", na.rm=T),
      diabetes_pct = round(diabetes_n/diabetes_known_n*100, round.var),
      cv_dis_known_n = sum(!is.na(CV_DIS)),
      cv_dis_n = sum(CV_DIS == "1", na.rm=T),
      cv_dis_pct = round(cv_dis_n/cv_dis_known_n*100, round.var),
      heart_dis_known_n = sum(!is.na(Heart_Dis)),
      heart_dis_n = sum(Heart_Dis == "1", na.rm=T),
      heart_dis_pct = round(heart_dis_n/heart_dis_known_n*100, round.var),
      
      act_cohort_known = sum(!is.na(cohort)),
      act_cohort1_n = sum(cohort == "1", na.rm=T),
      act_cohort1_pct = round(act_cohort1_n/act_cohort_known*100, round.var),
      act_cohort2_n = sum(cohort == "2", na.rm=T),
      act_cohort2_pct = round(act_cohort2_n/act_cohort_known*100, round.var),
      act_cohort3_n = sum(cohort == "3", na.rm=T),
      act_cohort3_pct = round(act_cohort3_n/act_cohort_known*100, round.var),
      
      casi_irt_mean = round(mean(casi_irt, na.rm=T), 2), 
      casi_irt_sd = round(sd(casi_irt, na.rm=T), 2),
      
      model_wt_mean_sd = qwraps2::mean_sd(model_wt, digits = 2, na_rm = T, denote_sd = "paren")
    ) 
  
  t1 <- t1 %>%
    mutate(
      "Participants (n, %)" = paste0(N, " (", N_perc, "%)" ),
      "Person-years (n, %)" = pys_total,
      
      "Entry age, years (median, IQR)" = paste0(age_entry_median, " (", age_entry_iqr, ")"),
      "Follow-up years (mean, SD)" = paste0(follow_up_yrs_mean, " (", follow_up_yrs_sd, ")"),
      
      #exposure
      "Exposure at Entry, ppb (mean, SD)" = exposure,
      
      # cases
      "Dementia Cases (n, %)" = dementia_cases,
      "Alzheimer's Disease Cases (n, %)" = ad_cases,
      
      # Demographics 
      "Birth Cohort (n, %)" = "",
      "1909 or earlier" =  paste0(birth_cohort1_n, " (", birth_cohort1_pct, "%)"),
      "1910-1914" =  paste0(birth_cohort2_n, " (", birth_cohort2_pct, "%)"),
      "1915-1919" =  paste0(birth_cohort3_n, " (", birth_cohort3_pct, "%)"),
      "1920-1924" =  paste0(birth_cohort4_n, " (", birth_cohort4_pct, "%)"), 
      "1925-1929" =  paste0(birth_cohort5_n, " (", birth_cohort5_pct, "%)"), 
      "1930-1934" =  paste0(birth_cohort6_n, " (", birth_cohort6_pct, "%)"), 
      "1935 or later" =  paste0(birth_cohort7_n, " (", birth_cohort7_pct, "%)"), 
      
      "Female (n, %)" = paste0(female_n, " (", female_pct, "%)"),
      "White Race (n, %)" = paste0(race_white_n, " (", race_white_pct, "%)"),
      "Education (n, %)" = "",
      "Less than High School" = paste0(degree0_n, " (", degree0_pct, "%)"),
      "High School or GED" = paste0(degree1_n, " (", degree1_pct, "%)"),
      "Bachelor's Degree" = paste0(degree3_n, " (", degree3_pct, "%)"),
      "Master's Degree" = paste0(degree4_n, " (", degree4_pct, "%)"),
      "Doctorate Degree" = paste0(degree5_n, " (", degree5_pct, "%)"),
      "Other" = paste0(degree6_n, " (", degree6_pct, "%)"),
      "Census Tract Income (n, %)" = "",
      "< $35k" = paste0(income1_n, " (", income1_pct, "%)"),
      "$35-50k" = paste0(income2_n, " (", income2_pct, "%)"),
      "$50-75k" = paste0(income3_n, " (", income3_pct, "%)"),
      "> $75k" = paste0(income4_n, " (", income4_pct, "%)"),
      
      # Health Indicators
      "Smoker (n, %)" = "",
      "Never" = paste0(smoke_never_n, " (", smoke_never_pct, "%)"),
      "Former" = paste0(smoke_former_n, " (", smoke_former_pct, "%)"),
      "Current" = paste0(smoke_current_n, " (", smoke_current_pct, "%)"),
      "Regular exercise (n, %)" = paste0(exercise_regular_n, " (", exercise_regular_pct, "%)"),
      "BMI (n, %)" = "",
      "Underweight" = paste0(bmi_under_n, " (", bmi_under_pct, "%)"),
      "Normal" = paste0(bmi_normal_n, " (", bmi_normal_pct, "%)"),
      "Overweight" = paste0(bmi_over_n, " (", bmi_over_pct, "%)"),
      "Obese" = paste0(bmi_obese_n, " (", bmi_obese_pct, "%)"),
      "Vascular Health (n, %)" = "",
      "Hypertension" = paste0(hypertension_n, " (", hypertension_pct, "%)"),
      "Diabetes" = paste0(diabetes_n, " (", diabetes_pct, "%)"),
      "Cardiovascular Disease" = paste0(cv_dis_n, " (", cv_dis_pct, "%)"),
      "Heart Disease" = paste0(heart_dis_n, " (", heart_dis_pct, "%)"),
      
      # Cognitive fn indicator
      "APOE carrier (n, %)" = paste0(apoe_carrier_n, " (", apoe_carrier_pct, "%)"),
      "Baseline CASI (mean, SD)" = paste0(casi_irt_mean, " (", casi_irt_sd, ")"),
      
      # Study Indicator
      "ACT Cohort (n, %)" = "",
      "Original" =  paste0(act_cohort1_n, " (", act_cohort1_pct, "%)"),
      "Expansion" =  paste0(act_cohort2_n, " (", act_cohort2_pct, "%)"),
      "Replacement" =  paste0(act_cohort3_n, " (", act_cohort3_pct, "%)"),
      
      #IPW weights
      "IPW Weight, (mean, SD)" = model_wt_mean_sd
      
    ) %>%
    #get rid of repeat columns
    select(-c(N:model_wt_mean_sd #, N_perc, follow_up_yrs_mean:model_wt_mean_sd
    ))  
  
  # transpose
  t1 <- t(t1) %>% 
    as.data.frame() %>% 
    # separate N & %
    #separate(col = V1, into = c("N", "pct"), sep = " ", fill = "right") %>%
    #rename_all(~paste0(column.name, "_", .)) %>%
    
    rownames_to_column(var = "Variable")  
  
  # replace NAs w/ ""
  #t1[3][is.na(t1[3])] <- ""
  
  # #rename "V1" column  
  # if(column.name != "") {
  #   names(t1) <- column.name
  #   }
  
  return(t1) 
}
