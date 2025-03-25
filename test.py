#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:41:48 2025

@author: terrats
"""
#%% Import libraries

# working_directory = "/data/home/terrats/" # '/data/home/terrats/' # on rs-pro machine  
path_to_RIOMAR_package = "/home/terrats/Desktop/RIOMAR/PACKAGE/myRIOMAR_dev"
# path_to_RIOMAR_package = "/media/shared_storage/Documents/CODING/PROJECTS/myRIOMAR_dev"

import sys
sys.path.append(path_to_RIOMAR_package)

import myRIOMAR_dev as myRIOMAR

import matplotlib

# Set matplotlib backend to prevent plots from displaying
# mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM memory !)
matplotlib.use('agg') # Prevent showing plot in the Plot panel (this saves RAM memory)


#%% Set core parameters
core_arguments = {  'Data_sources': ['ODATIS'], # ODATIS / SEXTANT / EUMETSAT
                    'Sensor_names':["modis"], 
                    'Satellite_variables':['NRRS412', 'NRRS443', 'NRRS488', 'NRRS531', 'NRRS547', 'NRRS551', 'NRRS667', 'NRRS748', 'SPM-G'], 
                    'Atmospheric_corrections':["polymer", "nirswir"], # Standard (SEXTANT / EUMETSAT) /  ...
                    'Temporal_resolution': ["DAILY"], # A enlever ! 
                    'start_day' : '2002/07/04',
                    'end_day' : '2021/01/31'}

# core_arguments = {  'Data_sources': ['SEXTANT'], # ODATIS / SEXTANT / EUMETSAT
#                     'Sensor_names':["merged"], 
#                     'Satellite_variables':['CHLA', 'SPM'], 
#                     'Atmospheric_corrections':["Standard"], # Standard (SEXTANT / EUMETSAT) /  ...
#                     'Temporal_resolution': ["DAILY"], # A enlever ! 
#                     'start_day' : '1998/01/01',
#                     'end_day' : '2024/12/31'}

#%% Downlaod data

myRIOMAR._0_data_downloading.Download_satellite_data(core_arguments,
                                                     nb_of_cores_to_use = 3,
                                                     overwrite_existing_satellite_files = False,
                                                     # where_to_save_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA"
                                                     where_to_save_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/")

#%% Downlaod data
myRIOMAR._0_data_downloading.Plot_and_Save_the_map(core_arguments,
                                                   nb_of_cores_to_use = 6,
                                                   where_are_saved_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA",
                                                   start_day_of_maps_to_plot = '2018/07/31',
                                                   end_day_of_maps_to_plot = '2018/08/01')

#%%

myRIOMAR._Add_my_zone(name = 'THAY', min_lat = 9, max_lat = X)

myRIOMAR._1_data_validation.Match_up_with_insitu_measurements(core_arguments, 
                                                              redo_the_MU_database = False, 
                                                              nb_of_cores_to_use = 6,
                                                              where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
                                                              where_to_save_Match_Up_data = "/home/terrats/Desktop/RIOMAR/TEST/MATCH_UP_DATA/")


#%%

myRIOMAR._2_regional_maps.create_regional_maps(core_arguments,
                                               Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY'],
                                               overwrite_existing_regional_maps = True,
                                               save_map_plots_of_which_time_frequency = {'DAILY' : False, 'WEEKLY' : False, 'MONTHLY' : True, 'ANNUAL' : True},
                                               nb_of_cores_to_use = 4,
                                               where_are_saved_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA",
                                               where_to_save_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS")


#%%

myRIOMAR._2_regional_maps.QC_of_regional_maps(core_arguments,
                                               Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY'],
                                               nb_of_cores_to_use = 4,
                                               where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS")

#%%

myRIOMAR._3_plume_detection.apply_plume_mask(core_arguments,
                                             Zones = ['GULF_OF_LION'],
                                             detect_plumes_on_which_temporal_resolution_data = 'WEEKLY',
                                             nb_of_cores_to_use = 6,
                                             where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS")

# myRIOMAR.plume_detection.plot_time_series_of_plume_areas(work_dir = work_dir,
#                                                      Zones = ['BAY_OF_SEINE'],
#                                                      Data_sources = ['SEXTANT'],
#                                                      Satellite_sensors = ['modis', 'merged'],
#                                                      Atmospheric_corrections = ['Standard'],
#                                                      Years = range(2018,2022),
#                                                      Time_resolutions = ['DAILY'])