#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:29:09 2025

@author: terrats
"""

import rpy2.robjects as robjects
import os, pickle, re

from utils import (load_csv_files_in_the_package_folder, path_to_fill_to_where_to_save_satellite_files, align_bathymetry_to_resolution)
from _5_Figures_for_article.utils import save_files_for_Figure_1, load_the_regional_maps_and_save_them_for_plotting
from _3_plume_detection.utils import (define_parameters, reduce_resolution,  
                                      preprocess_annual_dataset_and_compute_land_mask, create_polygon_mask,
                                      Create_the_plume_mask)

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
    

def Figure_4(where_are_saved_regional_maps, where_to_save_the_figure) : 
    
    Zone = 'BAY_OF_SEINE'
    Date = '2005-01-04'
    
    parameters = define_parameters(Zone)

    path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'RESULTS', Zone, 'SEXTANT', 'SPM', 'merged', 
                                                     'Standard', 'MAPS', 'DAILY', Date[:4], f'{Date}.pkl')

    # Open and load the file (binary file assumed to contain data)
    with open(path_to_the_satellite_file_to_use, 'rb') as f:
        ds = pickle.load(f)['Basin_map']['map_data']  

    # Reduce the resolution of the dataset to the specified latitude and longitude resolutions
    ds_reduced = (reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution']) 
                  if parameters['lat_new_resolution'] is not None
                  else ds) 

    bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced, f'{where_are_saved_regional_maps}/RESULTS/{Zone}/Bathy_data.pkl')

    (_, land_mask) = preprocess_annual_dataset_and_compute_land_mask( (path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + "/RESULTS/" + Zone)
                                                                       .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]', 
                                                                                'SEXTANT/SPM/merged/Standard/MAPS/MULTIYEAR')
                                                                       .replace('[YEAR]/[MONTH]/[DAY]', 'Averaged_over_multi-years.pkl')
                                                                       ), parameters)

    inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)
    
    the_plume = Create_the_plume_mask(ds_reduced, 
                                       bathymetry_data_aligned_to_reduced_map,
                                       land_mask,
                                       parameters,
                                       'Seine')
    the_plume.do_R_plot( where_to_save_the_plot = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURES', 'FIGURE_4'), 
                         name_of_the_plot = 'A' )
    
    the_plume.determine_SPM_threshold()
    the_plume.do_R_plot()

    the_plume.do_a_raw_plume_detection()
    the_plume.do_R_plot()

    the_plume.include_cloudy_regions_to_plume_area()

    the_plume.remove_the_areas_with_sediment_resuspension(maximal_bathymetry = parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
                                                          minimal_distance_from_estuary = parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][plume_name])

    the_plume.remove_shallow_waters()
    the_plume.do_R_plot()

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])
    the_plume.do_R_plot()                
    
    the_plume.dilate_the_main_plume_area_to_merge_close_plume_areas()

    the_plume.remove_small_shapes_that_do_not_meet_a_minimum_size_criterion()

    the_plume.set_pixels_to_False_if_outside_of_the_searching_area(inside_polygon_mask)

    the_plume.identify_the_main_plume_shape_based_on_the_plume_core_location()

    the_plume.remove_shallow_waters()

    if not np.isin(plume_name, ['Seine']) :
        the_plume.remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase()

    the_plume.do_R_plot()
