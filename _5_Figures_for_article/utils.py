#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:31:25 2025

@author: terrats
"""

import numpy as np
import xarray as xr
import pandas as pd

import glob, os, pickle

from _99_common.utils import path_to_fill_to_where_to_save_satellite_files, define_parameters
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
          
    extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_1_data)
    
    
def load_the_regional_maps_and_save_them_for_plotting(where_are_saved_regional_maps, where_to_save_the_figure, dates_for_each_zone) :
        
    folder_where_to_save_Figure_2_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURES', 'FIGURE_2', 'DATA')
    os.makedirs(folder_where_to_save_Figure_2_data, exist_ok = True)
    
    path_to_regional_maps = {key : (path_to_fill_to_where_to_save_satellite_files( os.path.join(where_are_saved_regional_maps, 'RESULTS', key) )
                                       .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                                                'SEXTANT/SPM/merged/Standard/MAPS/DAILY')
                                       .replace('[YEAR]/[MONTH]/[DAY]', f'{date[:4]}/{date}.pkl')) 
                               for key, date in dates_for_each_zone.items()}
    
    for key, path_to_map in path_to_regional_maps.items() : 
    
        coordinates_of_the_map = define_parameters(key)    
    
        with open(path_to_map, 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data']  
            
            SPM_map = (ds
                       .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0], 
                                      coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]), 
                            lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                      coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])) 
                       .to_dataframe()
                       .reset_index())
            
            SPM_map.to_csv(folder_where_to_save_Figure_2_data + f"/{key}.csv")
            
    extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_2_data)
            
def extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_data) :     

    coordinates_of_the_RIOMARS = { zone_name : {'lat_min' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][0],
                                                'lat_max' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][1],
                                                'lon_min' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][0],
                                                'lon_max' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][1]}
                                  for zone_name in ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE'] }
    
    pd.DataFrame.from_dict(coordinates_of_the_RIOMARS).to_csv(folder_where_to_save_Figure_data + "/RIOMAR_limits.csv")
    
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
    
    stations_in_the_RIOMARS.to_csv(folder_where_to_save_Figure_data + "/Stations_position.csv")
    
    
def dates_for_each_zone() : 
    return {'GULF_OF_LION' : '2014-01-04',
            'BAY_OF_BISCAY' : '2009-04-22',
            'SOUTHERN_BRITTANY' : '2016-05-23',# '2022-01-21',
            'BAY_OF_SEINE' : '2018-02-25'}