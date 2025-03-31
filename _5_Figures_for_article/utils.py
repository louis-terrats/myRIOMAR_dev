#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:31:25 2025

@author: terrats
"""

import numpy as np
import xarray as xr
import pandas as pd

import glob, os

from utils import path_to_fill_to_where_to_save_satellite_files
from _3_plume_detection.utils import define_parameters
from _1_data_validation.utils import get_insitu_measurements
                   
def save_files_for_Figure_1(where_are_saved_satellite_data, where_to_save_the_figure, 
                            date_of_the_map, coordinates_of_the_map) :

    folder_where_to_save_Figure_1_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURES', 'FIGURE_1', 'DATA')
    os.makedirs(folder_where_to_save_Figure_1_data, exist_ok = True)
    
    path_to_nc_file = (path_to_fill_to_where_to_save_satellite_files(where_are_saved_satellite_data)
                       .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                                'SEXTANT/SPM/merged/Standard/DAILY')
                       .replace('[YEAR]/[MONTH]/[DAY]',
                                date_of_the_map))
        
    with xr.open_dataset( glob.glob(path_to_nc_file + "/*.nc")[0] ) as ds :
        SPM_map = (ds['analysed_spim']
                   .sel(lat=slice(coordinates_of_the_map['lat_min'], coordinates_of_the_map['lat_max']), 
                        lon=slice(coordinates_of_the_map['lon_min'], coordinates_of_the_map['lon_max'])) 
                   .to_dataframe()
                   .reset_index()
                   .drop(columns=["time"]))
        
        SPM_map.to_csv(folder_where_to_save_Figure_1_data + "/SPM_map.csv")
            
    coordinates_of_the_RIOMARS = { zone_name : {'lat_min' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][0],
                                                'lat_max' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][1],
                                                'lon_min' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][0],
                                                'lon_max' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][1]}
                                  for zone_name in ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE'] }
    
    pd.DataFrame.from_dict(coordinates_of_the_RIOMARS).to_csv(folder_where_to_save_Figure_1_data + "/RIOMAR_limits.csv")
    
    _, insitu_stations = get_insitu_measurements()
    station_LATITUDES = insitu_stations.LATITUDE.to_numpy(dtype=float)
    station_LONGITUDES = insitu_stations.LONGITUDE.to_numpy(dtype=float)
    
    index_of_stations_in_the_RIOMARS = [ np.where((station_LATITUDES >= coords['lat_min']) & 
                                                  (station_LATITUDES <= coords['lat_max']) & 
                                                  (station_LONGITUDES >= coords['lon_min']) & 
                                                  (station_LONGITUDES <= coords['lon_max']) )[0] 
                                        for coords in coordinates_of_the_RIOMARS.values() ]
    
    index_of_stations_in_the_RIOMARS = np.unique(np.concatenate(index_of_stations_in_the_RIOMARS))
    
    stations_in_the_RIOMARS = insitu_stations.iloc[index_of_stations_in_the_RIOMARS]
    
    stations_in_the_RIOMARS.to_csv(folder_where_to_save_Figure_1_data + "/Stations_position.csv")