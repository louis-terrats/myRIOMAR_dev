import pandas as pd

from _2_regional_maps.utils import Create_and_save_the_maps, QC_maps

from utils import (unique_years_between_two_dates, get_all_cases_to_process_for_regional_maps_or_plumes_or_X11)


def create_regional_maps(core_arguments, Zones, overwrite_existing_regional_maps, save_map_plots_of_which_time_frequency, nb_of_cores_to_use,
                         where_are_saved_satellite_data, where_to_save_regional_maps) :
            
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones})
    
    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)
    
    for i, info in cases_to_process.iterrows() : 
                
        # info = cases_to_process.iloc[i]
        info = pd.concat([info, pd.Series([core_arguments['start_day'] if info.Year == core_arguments['Years'][0] else f'{info.Year}/01/01',
                                           core_arguments['end_day'] if info.Year == core_arguments['Years'][-1] else f'{info.Year}/12/31'],
                                          index = ['date_min', 'date_max'])])
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        maps_creation = Create_and_save_the_maps(where_to_save_regional_maps, where_are_saved_satellite_data, info) 
        
        if ('map_files' not in vars(maps_creation)) or len(maps_creation.map_files) == 0 : 
            print("Switch to the next iterate")
            continue
        
        if maps_creation.are_the_maps_already_produced and (overwrite_existing_regional_maps == False) : 
            print("Maps already exist - Switch to the next iterate")
            continue       
        
        maps_creation._1_create_weekly_maps(nb_of_cores_to_use, save_map_plots_of_which_time_frequency)
        
        maps_creation._2_create_monthly_maps(nb_of_cores_to_use, save_map_plots_of_which_time_frequency)
        
        maps_creation._3_create_annual_maps(save_map_plots_of_which_time_frequency)
        
    global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates()
        
    for i, info in global_cases_to_process.iterrows() :  
        
        # info = global_cases_to_process.iloc[i].copy()
        info['Year'] = 'MULTIYEAR'
        info = pd.concat([info, pd.Series([core_arguments['start_day'], core_arguments['end_day']], index = ['date_min', 'date_max'])])
        
        maps_creation = Create_and_save_the_maps(where_to_save_regional_maps, where_are_saved_satellite_data, info) 
        
        maps_creation._4_create_the_multiyear_map()
        

        
def QC_of_regional_maps(core_arguments, Zones, nb_of_cores_to_use, where_are_saved_regional_maps) : 
    
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones})
    
    cases_to_process = get_all_possibilities(core_arguments)
    
    for i, info in cases_to_process.iterrows() :  
                
        # info = cases_to_process.iloc[i].copy()
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        QC_data_and_plots = QC_maps(info, where_are_saved_regional_maps)
        
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
       
        QC_data_and_plots = QC_maps(info, where_are_saved_regional_maps)
               
        if len(QC_data_and_plots.map_files) == 0 : 
            continue
       
        try : 
            QC_data_and_plots.combine_QC_metrics()
        except Exception : 
            continue
       
