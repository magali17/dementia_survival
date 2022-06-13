###################################
# Author: Magali Blanco
# Date: 6/4/2022
# Purpose: Scritp for common functions
###################################

###################################
# SETUP
###################################
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

# pacman::p_load(tidyverse,
#                lubridate
# )    

set.seed(1)

###################################
# Lasso
###################################

#load library for lasso
pacman::p_load(glmnet)


# dt = apoe
# x_names = names(apoe)[2:(ncol(apoe)-1)]
# y_name = "apoe_available"
# family. = "binomial"
# lambda. = ""


lasso_fn <- function(dt, x_names, y_name, family. = "gaussian", lambda. = "") {
  
  # # don't include variables that don't vary
  # no_variation <- apply(apoe, 2, function(x) length(unique(x))) %>%
  #   data.frame() %>%
  #   filter(.==1) %>%
  #   rownames()
  # x_names <- setdiff(x_names, no_variation)
  
    
  x <- model.matrix(as.formula(paste(y_name, "~", paste(x_names, collapse = " + "))), dt)[, -1]
  y <- dt[[y_name]]   
  
  # nrow(x)
  # length(y)
  
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
