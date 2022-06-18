######################################################################
# Author: Magali Blanco
# Date: 6/4/2022
# Purpose: Script for common functions
######################################################################

######################################################################
# Lasso
######################################################################
pacman::p_load(glmnet)

lasso_fn <- function(dt, x_names, y_name, family. = "gaussian", lambda. = "") {
    
  x <- model.matrix(as.formula(paste(y_name, "~", paste(x_names, collapse = " + "))), dt)[, -1]
  y <- dt[[y_name]]   
  
  #select lambda through CV if not supplied
  if(lambda. == ""){
    cv.out <- cv.glmnet(x = x, y = y, alpha=1, family= family., standardize=T)
    lambda. <- cv.out$lambda.min
    }
  
  # run Lasso
  lasso.m <- glmnet(x = x, y = y, alpha = 1, family= family., standardize = T)
  
  #save coefficient estimates
  lasso_coef <- predict(lasso.m, type= "coefficients", s= lambda.)[1:(ncol(x)+1),] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(cov = rowname,
           coef = ".") %>%
    #keep coefficients that are not 0 or intercept values
    filter(coef != 0,
           cov != "(Intercept)")
  
  return(list(results = lasso_coef,
              lambda = lambda.))
  
}

######################################################################
# CODE FACTORS
######################################################################
# fn factorizes and creates reference categories 
#for model stability, the reference categories are those with the largest counts

add_factor_refs <- function(dt) {
  dt <- mutate(dt,
               apoe = relevel(factor(apoe), ref = "0"),
               male = relevel(factor(male), ref = "0"),
               degree = relevel(factor(degree), ref = "1"),
               income_cat = relevel(factor(income_cat), ref = "3"),
               race_white = relevel(factor(race_white), ref = "1"),
               cal_2yr = relevel(factor(cal_2yr), ref = "1996"), #was 2004 before
               )
  return(dt)
  }


######################################################################
#  
######################################################################




