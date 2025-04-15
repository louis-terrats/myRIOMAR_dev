#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:41:48 2025

@author: terrats
"""
#%% Import libraries

# working_directory = "/data/home/terrats/" # '/data/home/terrats/' # on rs-pro machine  
path_to_RIOMAR_package = "/home/terrats/Desktop/RIOMAR/PACKAGE/myRIOMAR_dev"
# path_to_RIOMAR_package = "/data/home/terrats/PACKAGE/myRIOMAR_dev"
# path_to_RIOMAR_package = "/media/shared_storage/Documents/CODING/PROJECTS/myRIOMAR_dev"

import sys
sys.path.append(path_to_RIOMAR_package)

import myRIOMAR_dev as myRIOMAR

import matplotlib

# Set matplotlib backend to prevent plots from displaying
# matplotlib.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM memory !)
matplotlib.use('agg') # Prevent showing plot in the Plot panel (this saves RAM memory)


#%% Set core parameters
core_arguments = {  'Data_sources': ['ODATIS'], # ODATIS / SEXTANT / EUMETSAT
                    'Sensor_names':["modis"], 
                    'Satellite_variables':['NRRS412', 'NRRS443', 'NRRS488', 'NRRS531', 'NRRS547', 'NRRS551', 'NRRS667', 'NRRS748', 'SPM-G'], 
                    'Atmospheric_corrections':["polymer"], # Standard (SEXTANT / EUMETSAT) /  ...
                    'Temporal_resolution': ["DAILY"], # A enlever ! 
                    'start_day' : '2002/07/04',
                    'end_day' : '2021/01/31'}

core_arguments = {  'Data_sources': ['SEXTANT'], # ODATIS / SEXTANT / EUMETSAT
                    'Sensor_names':["merged"], 
                    'Satellite_variables':['SPM'], 
                    'Atmospheric_corrections':["Standard"], # Standard (SEXTANT / EUMETSAT) /  ...
                    'Temporal_resolution': ["WEEKLY"], # A enlever ! 
                    'start_day' : '1998/01/01',
                    'end_day' : '2024/12/31'}

#%% Downlaod data

myRIOMAR._0_data_downloading.Download_satellite_data(core_arguments,
                                                     nb_of_cores_to_use = 3,
                                                     overwrite_existing_satellite_files = False,
                                                     # where_to_save_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA"
                                                     where_to_save_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/")

#%% Plot some maps to check their content ?
myRIOMAR._0_data_downloading.Plot_and_Save_the_map(core_arguments,
                                                   nb_of_cores_to_use = 6,
                                                   where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
                                                   start_day_of_maps_to_plot = '2018/07/31',
                                                   end_day_of_maps_to_plot = '2018/08/01')


#%% Add my study zone

# myRIOMAR._Add_my_zone(name = 'THAY', min_lat = 9, max_lat = X)

#%% Match-up with SOMLIT and REPHY stations

myRIOMAR._1_data_validation.Match_up_with_insitu_measurements(core_arguments, 
                                                              # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                                              zones = ['FRANCE'],
                                                              redo_the_MU_database = False, 
                                                              nb_of_cores_to_use = 6,
                                                              where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
                                                              where_to_save_Match_Up_data = "/home/terrats/Desktop/RIOMAR/MANH")
                                                              # where_to_save_Match_Up_data = "/home/terrats/Desktop/RIOMAR/TEST")


#%% Make regional maps

myRIOMAR._2_regional_maps.create_regional_maps(core_arguments,
                                               Zones = ['SOUTHERN_BRITTANY'],
                                               overwrite_existing_regional_maps = True,
                                               save_map_plots_of_which_time_frequency = {'DAILY' : False, 'WEEKLY' : False, 'MONTHLY' : True, 'ANNUAL' : True},
                                               nb_of_cores_to_use = 4,
                                               where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
                                               where_to_save_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/NEW_TEST")

#%% Do the QC of these regional maps

# myRIOMAR._2_regional_maps.QC_of_regional_maps(core_arguments,
#                                                Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY'],
#                                                nb_of_cores_to_use = 4,
#                                                where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS")

#%% Detect the plumes

myRIOMAR._3_plume_detection.apply_plume_mask(core_arguments,
                                             # Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                             Zones = ['SOUTHERN_BRITTANY'],
                                             detect_plumes_on_which_temporal_resolution_data = 'WEEKLY',
                                             nb_of_cores_to_use = 4,
                                             use_dynamic_threshold = True,
                                             where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS",
                                             where_to_save_plume_results = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD")


#%% Plot time-series of plume surface

myRIOMAR._3_plume_detection.make_and_plot_time_series_of_plume_areas(core_arguments,
                                                     Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                                     # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                                     nb_of_cores_to_use = 4,
                                                     on_which_temporal_resolution_the_plumes_have_been_detected = 'WEEKLY',
                                                     where_are_saved_plume_results = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/FIXED_THRESHOLD/",
                                                     where_to_save_plume_time_series = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/FIXED_THRESHOLD/")


#%% Do X11 analysis on plume time-series

myRIOMAR._4_X11_analysis.Apply_X11_method_on_time_series(core_arguments,
                                            Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                            # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                            nb_of_cores_to_use = 4,
                                            on_which_temporal_resolution_the_plumes_have_been_detected = 'WEEKLY',
                                            where_are_saved_plume_time_series = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD",
                                            where_to_save_X11_results = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD",
                                            include_river_flow = True)


#%% Figures_for_the_article

#%% Figure 1

myRIOMAR._5_Figures_for_article.Figure_1(where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
                                         where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

myRIOMAR._5_Figures_for_article.Figure_2(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
                                         where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

myRIOMAR._5_Figures_for_article.Figure_4(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
                                         where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

myRIOMAR._5_Figures_for_article.Figure_5(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
                                         where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

myRIOMAR._5_Figures_for_article.Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD/",
                                           where_are_saved_plume_results_with_fixed_threshold = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/FIXED_THRESHOLD/",
                                           where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

myRIOMAR._5_Figures_for_article.Figure_8_9_10(
                                           where_are_saved_X11_results = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD/",
                                           where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")
