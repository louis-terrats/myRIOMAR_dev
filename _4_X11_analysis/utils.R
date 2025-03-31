list_of_packages <- c("plyr", "tidyverse", "ggpubr", "viridis", "doParallel", "zoo", "ggnewscale", "ggpubr", "scales")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

