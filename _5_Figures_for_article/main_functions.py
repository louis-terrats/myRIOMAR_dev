#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:29:09 2025

@author: terrats
"""

import rpy2.robjects as robjects

from utils import (load_csv_files_in_the_package_folder)
from _5_Figures_for_article.utils import save_files_for_Figure_1

def Figure_1(where_are_saved_satellite_data, where_to_save_the_figure) : 
        
    save_files_for_Figure_1(where_are_saved_satellite_data, 
                            where_to_save_the_figure,
                            date_of_the_map = "2011/02/02",
                            coordinates_of_the_map = {"lat_min" : 41,"lat_max" : 51.5,"lon_min" : -7,"lon_max" : 9.5})
    
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_1']
    
    # Call the R function
    r_function(where_to_save_the_figure = robjects.StrVector([where_to_save_the_figure]))
    
    
# def Figure_2(where_are_saved_regional_maps, where_to_save_the_figure) : 
