#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:40:41 2025

@author: terrats
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from myRIOMAR_dev._4_X11_analysis.utils import (temporal_decomp_V2_7_x11, apply_X11_method_and_save_results)
from utils import (get_all_cases_to_process_for_regional_maps_or_plumes_or_X11, path_to_fill_to_where_to_save_satellite_files, 
                   fill_the_sat_paths, load_csv_files_in_the_package_folder)


### S'INSPIRER DE LA FONCTION DNAS LE DOSSIER SCRIPTS

def Apply_X11_method_on_time_series(core_arguments, Zones, nb_of_cores_to_use,
                                    on_which_temporal_resolution_the_plumes_have_been_detected,
                                    where_are_saved_plume_time_series,
                                    where_to_save_X11_results,
                                    include_river_flow = False) :

    core_arguments.update({'Zones' : Zones,
                           'Years' : "*",
                           'Satellite_variables': ['SPM'],
                           'Temporal_resolution' : ([on_which_temporal_resolution_the_plumes_have_been_detected] 
                                                    if isinstance(on_which_temporal_resolution_the_plumes_have_been_detected, str) 
                                                    else on_which_temporal_resolution_the_plumes_have_been_detected)})
    
    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)
    
    var_to_use = 'area_of_the_plume_mask_in_km2'
    
    for i, info in cases_to_process.iterrows() : 
        
        # info = cases_to_process.iloc[i]
            
        file_names_pattern = (fill_the_sat_paths(info, 
                                               path_to_fill_to_where_to_save_satellite_files(where_are_saved_plume_time_series + "/" + info.Zone),
                                               local_path = True)
                                  .replace(info.atmospheric_correction, f'{info.atmospheric_correction}/PLUME_DETECTION/')
                                  .replace('/*/*/*', '/Time_series_of_plume_area_and_SPM_threshold.csv'))
    
        if os.path.exists(file_names_pattern) == False : 
            print(f'File does not exists : {file_names_pattern}')
            continue
        
        ts_data = pd.read_csv(file_names_pattern)
        
        apply_X11_method_and_save_results(values = ts_data[var_to_use].tolist(), variable_name = var_to_use, 
                                          dates = ts_data.date, info = info, 
                                          where_to_save_X11_results = where_to_save_X11_results)
        
        if include_river_flow :
            
            river_flow_data = load_csv_files_in_the_package_folder(RIVER_FLOW = True, Zone_of_river_flow = info.Zone)
            
            apply_X11_method_and_save_results(values = river_flow_data.Values.tolist(), variable_name = 'river_flow', 
                                              dates = river_flow_data.Date, 
                                              info = info.replace({info.Satellite_variable : "River_flow",
                                                                   info.atmospheric_correction : "",
                                                                   info.sensor_name : "",
                                                                   info.Data_source : "River_flow"}), 
                                              where_to_save_X11_results = where_to_save_X11_results)
            

        
        