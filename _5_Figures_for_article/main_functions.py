#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:29:09 2025

@author: terrats
"""

import rpy2.robjects as robjects
import os

from utils import (load_csv_files_in_the_package_folder, path_to_fill_to_where_to_save_satellite_files)
from _5_Figures_for_article.utils import save_files_for_Figure_1, load_the_regional_maps_and_save_them_for_plotting

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
    
    
def Figure_2(where_are_saved_regional_maps, where_to_save_the_figure) : 

    dates_for_each_zone = {'GULF_OF_LION' : '2014-01-05',
                           'BAY_OF_BISCAY' : '2009-04-22',
                           'SOUTHERN_BRITTANY' : '2009-04-22',
                           'BAY_OF_SEINE' : '2018-02-25'}
    
    load_the_regional_maps_and_save_them_for_plotting(where_are_saved_regional_maps,
                                                      where_to_save_the_figure,
                                                      dates_for_each_zone)
    
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_2']
    
    # Call the R function
    r_function(where_to_save_the_figure = robjects.StrVector([where_to_save_the_figure]))
    
    
    
