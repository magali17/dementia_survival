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



# --> DELETE FN?
######################################################################
# compare indicator variables across pollutants
######################################################################
# dt = sur
# var= "exp_wks_coverage10_yr_"
compare_indicator_var <- function(var, dt = sur) {
  
  dt %>%
    select(pollutant, exposure_year, starts_with(var)) %>%
    pivot_longer(cols = starts_with(var), names_to = "var", values_to = "value") %>%
    drop_na(value) %>%  
    group_by(var, pollutant) %>%
    summarize(
      n = n(),
      min = min(value),
      mean = mean(value),
      median = median(value),
      max = max(value)
    )
  
}

######################################################################
# summarize a variable
######################################################################

summary_table <- function(dt, var) {

  dt %>%
    rename(var = var) %>%
    summarize(
      N = length(var),
      Min = min(var),
      Mean = mean(var),
      SD = sd(var),
      Median = median(var),
      IQR = IQR(var),
      Max = max(var)
    )
  
}

###############################################################################################################
# ALTERNATIVE BOXPLOTS
###############################################################################################################
# function returns a boxplot with different whisker definitions to avoid plotting extreme/outlier points
alt_boxplot <- function(df, var, min_q=0.025, max_q=0.975){
  df <- df %>%
    rename(var = var) %>%
    summarize(
      N = n(),
      Min = min(var),
      Qmin = quantile(var, min_q),
      Q25 = quantile(var, 0.25),
      Q50 = quantile(var, 0.50),
      Q75 = quantile(var, 0.75),
      Qmax = quantile(var, max_q),
      Max = max(var)
    )
  
  names(df)[names(df)==var] <- var
  
  return(df) 
  
}

######################################################################
# LABEL POLLUTANTS
######################################################################
# fn relabels pollutants, adds units

label_pollutants <- function(dt) {
  dt <- dt %>%
    mutate(pollutant = case_when(
      pollutant == "ufp_10_42" ~ "PNC (pt/cm3), NanoScan",
      pollutant == "ufp_20_1k" ~ "PNC (pt/cm3), P-TRAK",
      pollutant == "ufp_36_1k" ~ "PNC (pt/cm3), Screened P-TRAK",
      pollutant == "ufp_10_70" ~ "PNC (pt/cm3), DiSCmini",
      pollutant == "bc" ~ "BC (ng/m3)",
      pollutant == "no2" ~ "NO2 (ppb)",
      pollutant == "pm25" ~ "PM2.5 (ug/m3)",
      pollutant == "co2" ~ "CO2 (ppm)"
      ),
      pollutant = factor(pollutant, levels = c("PNC (pt/cm3), NanoScan",
                                               "PNC (pt/cm3), P-TRAK",
                                               "PNC (pt/cm3), Screened P-TRAK",
                                               "PNC (pt/cm3), DiSCmini",
                                               "BC (ng/m3)",
                                               "NO2 (ppb)",
                                               "PM2.5 (ug/m3)",
                                               "CO2 (ppm)"))
      )
  return(dt)
}


#######################################################################################################################
# facet_wrap_equal() and facet_grid_equal() functions act like facet_wrap() and facet_grid() in ggplot but it sets the axes ranges (min/max) of each facet to same scale so that the 1-1 line is always down the middle :D !!
#######################################################################################################################
# code source: https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/ 

FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


#same as above but for facet_grid()
FacetEqualGrid <- ggproto(
  "FacetEqualGrid", FacetGrid,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetGrid, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_grid_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_grid(...)
  
  ggproto(NULL, FacetEqualGrid,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


######################################################################
# COX MODELS
######################################################################






