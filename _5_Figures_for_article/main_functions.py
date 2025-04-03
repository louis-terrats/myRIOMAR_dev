#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:29:09 2025

@author: terrats
"""

import rpy2.robjects as robjects
import os, pickle, re
import numpy as np
from functools import reduce


from utils import (load_csv_files_in_the_package_folder, path_to_fill_to_where_to_save_satellite_files, align_bathymetry_to_resolution, define_parameters)
from _5_Figures_for_article.utils import save_files_for_Figure_1, load_the_regional_maps_and_save_them_for_plotting
from _3_plume_detection.utils import (reduce_resolution,  
                                      preprocess_annual_dataset_and_compute_land_mask, create_polygon_mask,
                                      Create_the_plume_mask, Pipeline_to_delineate_the_plume)

def Figure_1(where_are_saved_satellite_data, where_to_save_the_figure) : 
        
    save_files_for_Figure_1(where_are_saved_satellite_data, 
                            where_to_save_the_figure,
                            date_of_the_map = "2011/02/02",
                            coordinates_of_the_map = {"lat_min" : 42,"lat_max" : 51.5,"lon_min" : -6,"lon_max" : 8})
    
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_1']
    
    # Call the R function
    r_function(where_to_save_the_figure = robjects.StrVector([where_to_save_the_figure]))
    
    
def Figure_2(where_are_saved_regional_maps, where_to_save_the_figure) : 

    dates_for_each_zone = {'GULF_OF_LION' : '2014-01-04',
                           'BAY_OF_BISCAY' : '2009-04-22',
                           'SOUTHERN_BRITTANY' : '2022-01-21',
                           'BAY_OF_SEINE' : '2018-02-25'}
    
    load_the_regional_maps_and_save_them_for_plotting(where_are_saved_regional_maps,
                                                      where_to_save_the_figure,
                                                      dates_for_each_zone)
    
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_2']
    
    # Call the R function
    r_function(where_to_save_the_figure = robjects.StrVector([where_to_save_the_figure]))
    

def Figure_5(where_are_saved_regional_maps, where_to_save_the_figure) : 

    dates_for_each_zone = {'GULF_OF_LION' : '2014-01-04',
                           'BAY_OF_BISCAY' : '2009-04-22',
                           'SOUTHERN_BRITTANY' : '2022-01-21',
                           'BAY_OF_SEINE' : '2018-02-25'}
    
    where_to_save_the_figure_5 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_5")
    os.makedirs(where_to_save_the_figure_5, exist_ok=True)
    
    for Zone, Date in dates_for_each_zone.items() : 

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
                
        all_mask_area = []
        all_river_mouth_to_remove = []
        thresholds = {key: None for key in parameters['starting_points']}
        # Loop through each plume starting point to process plume detection
        for plume_name, starting_point in parameters['starting_points'].items() : 
                        
            the_plume = Pipeline_to_delineate_the_plume(ds_reduced, 
                                               bathymetry_data_aligned_to_reduced_map,
                                               land_mask,
                                               parameters,
                                               plume_name,
                                               inside_polygon_mask,
                                               use_dynamic_threshold = True)
            
            thresholds[plume_name] = the_plume.SPM_threshold
            all_mask_area.append(the_plume.plume_mask)
            if "close_river_mouth_mask" in vars(the_plume) : 
                all_river_mouth_to_remove.append(the_plume.close_river_mouth_mask)
            
        
        # Combine all detected plume areas using logical OR
        final_mask_area = reduce(np.logical_or, all_mask_area)
        final_close_river_mouth_area = reduce(np.logical_or, all_river_mouth_to_remove)
            
        coordinates_of_the_map = define_parameters(Zone)    
            
        final_mask_area = (final_mask_area
                               .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0], 
                                              coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]), 
                                    lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                              coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])) )
                   
        ds_reduced = (ds_reduced
                   .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0], 
                                  coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]), 
                        lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                  coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])) ) 

        SPM_map = ds_reduced.to_dataframe().reset_index()
        SPM_map['plume'] = final_mask_area.values.flatten()
        SPM_map.to_csv(where_to_save_the_figure_5 + f"/DATA/{Zone}.csv")
            
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_5']
    
    # Call the R function
    r_function(where_to_save_the_figure = robjects.StrVector([where_to_save_the_figure_5]))


def Figure_4(where_are_saved_regional_maps, where_to_save_the_figure) : 
    
    Zone = 'GULF_OF_LION'
    plume_name = 'Grand Rhone'
    Date = '2007_01_02'
        
    parameters = define_parameters(Zone)

    path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'RESULTS', Zone, 'SEXTANT', 'SPM', 'merged', 
                                                     'Standard', 'MAPS', 'WEEKLY', Date[:4], f'Averaged_over_{Date[5:]}.pkl')

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
    
    where_to_save_the_figure_4 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_4")
    
    
    the_plume = Create_the_plume_mask(ds_reduced, 
                                       bathymetry_data_aligned_to_reduced_map,
                                       land_mask,
                                       parameters,
                                       plume_name)
    the_plume.do_R_plot( where_to_save_the_plot = where_to_save_the_figure_4, 
                         name_of_the_plot = 'A' )
    
    the_plume.determine_SPM_threshold()
    the_plume.do_R_plot( where_to_save_the_plot = where_to_save_the_figure_4, 
                         name_of_the_plot = 'B' )

    the_plume.do_a_raw_plume_detection()
    the_plume.do_R_plot( where_to_save_the_plot = where_to_save_the_figure_4, 
                         name_of_the_plot = 'C' )
    
    the_plume.include_cloudy_regions_to_plume_area()
    
    the_plume.remove_the_areas_with_sediment_resuspension(maximal_bathymetry = parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
                                                          minimal_distance_from_estuary = parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][plume_name])

    the_plume.remove_shallow_waters()
    the_plume.do_R_plot( where_to_save_the_plot = where_to_save_the_figure_4, 
                         name_of_the_plot = 'D' )

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])
    
    the_plume.dilate_the_main_plume_area_to_merge_close_plume_areas()

    the_plume.remove_small_shapes_that_do_not_meet_a_minimum_size_criterion()

    the_plume.set_pixels_to_False_if_outside_of_the_searching_area(inside_polygon_mask)

    the_plume.identify_the_main_plume_shape_based_on_the_plume_core_location()

    the_plume.remove_shallow_waters()

    if not np.isin(plume_name, ['Seine']) :
        the_plume.remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase()

    the_plume.do_R_plot( where_to_save_the_plot = where_to_save_the_figure_4, 
                         name_of_the_plot = 'E' )
