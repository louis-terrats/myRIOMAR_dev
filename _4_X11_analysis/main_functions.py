#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:40:41 2025

@author: terrats
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from myRIOMAR_dev._4_X11_analysis.utils import (temporal_decomp_V2_7_x11)
from utils import get_all_cases_to_process_for_regional_maps_or_plumes_or_X11, path_to_fill_to_where_to_save_satellite_files, fill_the_sat_paths


### S'INSPIRER DE LA FONCTION DNAS LE DOSSIER SCRIPTS

def Apply_X11_method_on_time_series(core_arguments, Zones, nb_of_cores_to_use,
                                    on_which_temporal_resolution_the_plumes_have_been_detected,
                                    where_are_saved_plume_time_series,
                                    where_to_save_X11_results) :

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
    
        results = temporal_decomp_V2_7_x11(values = ts_data[var_to_use].tolist(), dates = pd.to_datetime(ts_data.date), 
                                           time_frequency = info.Temporal_resolution,
                                            filter_outlier=False, overall_cutoff=50, 
                                            out_limit=3, perc_month_limit=50, 
                                            var_stationary=False, lin_interpol=False, 
                                            cutoff_fill=30, season_test=True)
        
        
        fig, axs = plt.subplots(4, 1, figsize=(24, 14))
        
        axs[0].plot(results['7_dates'], results['8_values_ini'])
        axs[0].set_title('Initial values')
        axs[1].plot(results['7_dates'], results['10_Interannual_signal'])
        axs[1].set_title(f'Inter-annual signal ({round(results["1_variance_due_to_Interannual_signal"], 1)}%)')
        axs[2].plot(results['7_dates'], results['9_Seasonal_signal'])
        axs[2].set_title(f'Seasonal signal ({round(results["0_variance_due_to_Seasonal_signal"], 1)}%)')
        axs[3].plot(results['7_dates'], results['11_Residual_signal'])
        axs[3].set_title(f'Residual signal ({round(results["2_variance_due_to_Residual_signal"], 1)}%)')
        
        # Adjust layout       
        plt.tight_layout()
    
        folder_where_to_store_the_plot = os.path.join(where_to_save_X11_results, 'X11_analysis', var_to_use) 
        os.makedirs(folder_where_to_store_the_plot, exist_ok=True)
        file_name = "_".join(info.drop(["Year", "Satellite_variable"]).astype(str).values)

        plt.savefig(folder_where_to_store_the_plot + "/" + file_name + '.png')
        
        plt.close(fig)
        
        pd.DataFrame({'dates' : results['7_dates'],
                        'Raw_signal' : results['8_values_ini'],
                        'Interannual_signal' : results['10_Interannual_signal'],
                        'Seasonal_signal' : results['9_Seasonal_signal'],
                        'Residual_signal' : results['11_Residual_signal'],
                        'Variation_coefficient' : results['5_var_coeff'],
                        'Variance_due_to_Interannual_signal' : results['1_variance_due_to_Interannual_signal'],
                        'Variance_due_to_Seasonal_signal' : results['0_variance_due_to_Seasonal_signal'],
                        'Variance_due_to_Residual_signal' : results['2_variance_due_to_Residual_signal'],
                        'Monotonic_change' : results['18_Kendall_Sen_analyses_on_Interannual_signal']['Is_the_change_of_annual_values_monotonic'],
                        'Rate_of_Change' : results['18_Kendall_Sen_analyses_on_Interannual_signal']['Rate_of_change_of_annual_values_in_percentage_per_time']
                        }).to_csv(folder_where_to_store_the_plot + "/" + file_name + '.csv', index = False)
    
