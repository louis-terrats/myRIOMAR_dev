import multiprocessing
import pandas as pd

from _2_regional_maps.utils import Create_and_save_the_maps, QC_maps, get_all_possibilities

from utils import (store_arguments, unique_years_between_two_dates)


def create_regional_maps(arguments) :
        
    (Data_sources, 
     Sensor_names, 
     Satellite_variables, 
     Atmospheric_corrections,
     Temporal_resolution, 
     start_day, end_day, 
     working_directory, 
     where_to_save_satellite_data, 
     overwrite_existing_satellite_files,
     redo_the_MU_database,
     path_to_SOMLIT_insitu_data,
     path_to_REPHY_insitu_data,
     Zones,
     overwrite_existing_regional_maps,
     plot_the_daily_regional_maps,
     nb_of_cores_to_use) = store_arguments(arguments, return_arguments = True)
    
    Years = unique_years_between_two_dates(start_day, end_day)
    
    cases_to_process = get_all_possibilities(Zones, Data_sources, Sensor_names, Atmospheric_corrections, Years, Satellite_variables)
    
    for i, info in cases_to_process.iterrows() : 
                
        # info = cases_to_process.iloc[i]
        info = pd.concat([info, pd.Series([start_day if info.Year == Years[0] else f'{info.Year}/01/01',
                                           end_day if info.Year == Years[-1] else f'{info.Year}/12/31'],
                                          index = ['date_min', 'date_max'])])
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        maps_creation = Create_and_save_the_maps(working_directory, where_to_save_satellite_data, info) 
        
        if ('map_files' not in vars(maps_creation)) or len(maps_creation.map_files) == 0 : 
            print(f"No satellite file here : {where_to_save_satellite_data}/{info.Data_source}/' - Switch to the next iterate")
            continue
        
        if maps_creation.are_the_maps_already_produced and (overwrite_existing_regional_maps == False) : 
            print("Maps already exist - Switch to the next iterate")
            continue       
        
        maps_creation._1_create_weekly_maps(nb_of_cores_to_use, plot_the_daily_regional_maps)
        
        maps_creation._2_create_monthly_maps(nb_of_cores_to_use)
        
        maps_creation._3_create_annual_maps()
        
    global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates()
        
    for i, info in global_cases_to_process.iterrows() :  
        
        # info = global_cases_to_process.iloc[i].copy()
        info['Year'] = 'MULTIYEAR'
        info = pd.concat([info, pd.Series([start_day, end_day], index = ['date_min', 'date_max'])])
        
        maps_creation = Create_and_save_the_maps(working_directory, where_to_save_satellite_data, info) 
        
        maps_creation._4_create_the_multiyear_map()
        

        
def QC_of_regional_maps(arguments) : 
    
    (Data_sources, 
     Sensor_names, 
     Satellite_variables, 
     Atmospheric_corrections,
     Temporal_resolution, 
     start_day, end_day, 
     working_directory, 
     where_to_save_satellite_data, 
     overwrite_existing_satellite_files,
     redo_the_MU_database,
     path_to_SOMLIT_insitu_data,
     path_to_REPHY_insitu_data,
     Zones,
     overwrite_existing_regional_maps,
     plot_the_daily_regional_maps,
     nb_of_cores_to_use) = store_arguments(arguments, return_arguments = True)
    
    Years = unique_years_between_two_dates(start_day, end_day)
                                            
    cases_to_process = get_all_possibilities(Zones, Data_sources, Sensor_names, Atmospheric_corrections, 
                                             Years, Satellite_variables)
    
    for i, info in cases_to_process.iterrows() :  
                
        # info = cases_to_process.iloc[i].copy()
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        QC_data_and_plots = QC_maps(info, working_directory)
        
        if len(QC_data_and_plots.map_files) == 0 : 
            print("No satellite regional maps - Switch to the next iterate")
            continue
                
        try : 
            QC_data_and_plots.compute_mask_for_coastal_waters(minimal_bathymetry_in_m = 1000,
                                                              minimal_distance_from_land_in_km = 20)
        except Exception : 
            print('End of the process - Switch to the next iterate')
            continue
        
        QC_data_and_plots.compute_QC_metrics(nb_of_cores_to_use, exclude_coastal_areas = True)
        
        
        
    global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates()
        
    for i in range(global_cases_to_process.shape[0]) : 
        
        info = global_cases_to_process.iloc[i].copy()
        info['Year'] = '*'
        
        # print(f'{i} over {global_cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
       
        QC_data_and_plots = QC_maps(info, working_directory)
               
        if len(QC_data_and_plots.map_files) == 0 : 
            continue
       
        try : 
            QC_data_and_plots.combine_QC_metrics()
        except Exception : 
            continue
       
