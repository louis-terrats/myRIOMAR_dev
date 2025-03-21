#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:41:48 2025

@author: terrats
"""

# =============================================================================
# ### Import libraries and set paths
# =============================================================================

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

# =============================================================================
# ### Apply the functions
# =============================================================================    

core_arguments = {  'Data_sources': ['EUMETSAT'],
                    'Sensor_names':["olcia"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':["polymer"],
                    'Temporal_resolution': ["DAILY"],
                    'start_day' : '2018/07/31',
                    'end_day' : '2018/08/01'}


myRIOMAR._0_data_downloading.Download_satellite_data(core_arguments,
                                                     nb_of_cores_to_use = 2,
                                                     overwrite_existing_satellite_files = False,
                                                     where_to_save_satellite_data = "/media/shared_storage/Documents/CODING/TEST/SATELLITE_DATA")

myRIOMAR._0_data_downloading.Plot_and_Save_the_map(core_arguments,
                                                   nb_of_cores_to_use = 6,
                                                   where_are_saved_satellite_data = "/media/shared_storage/Documents/CODING/TEST/SATELLITE_DATA",
                                                   start_day_of_maps_to_plot = '2018/07/31',
                                                   end_day_of_maps_to_plot = '2018/08/01')

myRIOMAR._1_data_validation.Match_up_with_insitu_measurements(core_arguments, 
                                                              redo_the_MU_database = False, 
                                                              nb_of_cores_to_use = 6,
                                                              where_to_save_Match_Up_data = "/home/terrats/Desktop/RIOMAR/TEST/MATCH_UP_DATA/")


myRIOMAR._2_regional_maps.create_regional_maps(core_arguments,
                                               Zones = ['BAY_OF_BISCAY', 'BAY_OF_SEINE', 'GULF_OF_LION'],
                                               overwrite_existing_regional_maps = True,
                                               plot_the_daily_regional_maps = False)

# myRIOMAR._2_regional_maps.QC_of_regional_maps(arguments.update({'nb_of_cores_to_use' : 4,
#                                                                  "Zones" : ['BAY_OF_BISCAY']}) or arguments)


myRIOMAR._3_plume_detection.apply_plume_mask(arguments.update({'nb_of_cores_to_use' : 4,
                                                               'Temporal_resolution' : ['WEEKLY']}) or arguments)

# myRIOMAR.plume_detection.plot_time_series_of_plume_areas(work_dir = work_dir,
#                                                      Zones = ['BAY_OF_SEINE'],
#                                                      Data_sources = ['SEXTANT'],
#                                                      Satellite_sensors = ['modis', 'merged'],
#                                                      Atmospheric_corrections = ['Standard'],
#                                                      Years = range(2018,2022),
#                                                      Time_resolutions = ['DAILY'])