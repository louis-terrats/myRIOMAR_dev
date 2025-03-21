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

import sys
sys.path.append(path_to_RIOMAR_package)

import myRIOMAR_dev as myRIOMAR

import multiprocessing, matplotlib

# Set matplotlib backend to prevent plots from displaying
# mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM memory !)
matplotlib.use('agg') # Prevent showing plot in the Plot panel (this saves RAM memory)

# =============================================================================
# ### Apply the functions
# =============================================================================    

core_arguments = {  'Data_sources': ['SEXTANT'],
                    'Sensor_names':["modis", "merged"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':["polymer"],
                    'Temporal_resolution': ["DAILY"],
                    'start_day' : '2018/07/31',
                    'end_day' : '2018/08/31'}


myRIOMAR._0_data_downloading.Download_satellite_data(core_arguments,
                                                     nb_of_cores_to_use = 2,
                                                     overwrite_existing_satellite_files = False,
                                                     where_to_save_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA/")

myRIOMAR._0_data_downloading.Plot_and_Save_the_map(core_arguments,
                                                   nb_of_cores_to_use = 6,
                                                   where_are_saved_satellite_data = "/home/terrats/Desktop/RIOMAR/TEST/SATELLITE_DATA/",
                                                   start_day_of_maps_to_plot = '2018/07/31',
                                                   end_day_of_maps_to_plot = '2018/08/05')

myRIOMAR._1_data_validation.Match_up_with_insitu_measurements(core_arguments, 
                                                              redo_the_MU_database = False, 
                                                              nb_of_cores_to_use = 6,
                                                              where_to_save_Match_Up_data = "/home/terrats/Desktop/RIOMAR/TEST/MATCH_UP_DATA/")


myRIOMAR._2_regional_maps.create_regional_maps(arguments.update({'nb_of_cores_to_use' : 4}) or arguments,
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


# #### BAY OF BISCAY

 
# #### BAY OF SEINE
# myRIOMAR.plume_detection.apply_plume_mask(work_dir = work_dir, 
#                                           Zone = 'BAY_OF_SEINE',
#                                           Data_Source = 'SEXTANT',
#                                           Satellite_sensor = 'merged',
#                                           Atmospheric_correction = 'Standard',
#                                           Years = range(2018,2019),
#                                           Time_resolution = 'WEEKLY')


# #### BAY OF SEINE
# myRIOMAR.plume_detection.apply_plume_mask(work_dir = work_dir, 
#                                           Zone = 'BAY_OF_SEINE',
#                                           Data_Source = 'SEXTANT',
#                                           Satellite_sensor = 'modis',
#                                           Atmospheric_correction = 'Standard',
#                                           Years = range(2018,2019),
#                                           Time_resolution = 'WEEKLY')



# #### GULF OF LION
# myRIOMAR.plume_detection.apply_plume_mask(work_dir = work_dir,
#                 Zone = 'GULF_OF_LION',
#                 Data_Source = 'SEXTANT',
#                 Satellite_sensor = 'merged',
#                 Atmospheric_correction = 'Standard',
#                 Years = range(1998,1999),
#                 Time_resolution = 'WEEKLY')



# #### EASTERN CHANNEL
# myRIOMAR.plume_detection.apply_plume_mask(work_dir = work_dir, 
#                   Zone = 'EASTERN_CHANNEL',
#                   Data_Source = 'ODATIS',
#                   Satellite_sensor = 'modis',
#                   Atmospheric_correction = 'nirswir',
#                   Years = range(1998,2024),
#                   Time_resolution = 'WEEKLY')


# #### SOUTHERN BRITTANY
# myRIOMAR.plume_detection.apply_plume_mask(work_dir = work_dir, 
#                   Zone = 'SOUTHERN_BRITTANY',
#                   Data_Source = 'ODATIS',
#                   Satellite_sensor = 'modis',
#                   Atmospheric_correction = 'nirswir',
#                   Years = range(1998,2024),
#                   Time_resolution = 'WEEKLY')


# # #### ETANG DE BERRE
# # apply_plume_mask(work_dir = work_dir,
# #                   Zone = 'ETANG_DE_BERRE',
# #                   Data_Source = 'SAMUEL',
# #                   Satellite_sensor = '',
# #                   Atmospheric_correction = '',
# #                   Years = np.arange(2018,2019),
# #                   Time_resolution = 'WEEKLY',
# #                 lon_new_resolution = 0.005,
# #                 lat_new_resolution = 0.005,
# #                 searching_strategies = {'EDF' : {'grid' : np.array([  [False, False, False, False, False],
# #                                                                         [False, False, False, False, False],
# #                                                                         [False, False, True,  False, False],
# #                                                                         [True,  True,  True,  True,  False],
# #                                                                         [False, False, False, False, False],
# #                                                                         ]),
# #                                               'coordinates_of_center' : (2,2)}},

# #                 SPM_threshold = 5,
# #                 bathymetric_threshold = 1,
# #                 starting_points = {'EDF' : (43.527, 5.071)}, # lat, lon 
# #                 core_of_the_plumes = {'EDF' : (43.520, 5.06)},
# #                 lat_range_of_the_area_to_check_for_clouds = [43.484, 43.533],
# #                 lon_range_of_the_area_to_check_for_clouds = [5.038, 5.130],
# #                 threshold_of_cloud_coverage_in_percentage = 25,
# #                 lat_range_of_the_map_to_plot = [43.37, 43.57],
# #                 lon_range_of_the_map_to_plot = [4.97, 5.26],
# #                 lat_range_to_search_plume_area = [-90, 90],
# #                 lon_range_to_search_plume_area = [-180, 5.10],
# #                 use_L4_maps = False)
