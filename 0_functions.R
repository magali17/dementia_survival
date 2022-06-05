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
#  
###################################
