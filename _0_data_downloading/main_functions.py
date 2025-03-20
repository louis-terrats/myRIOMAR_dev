import pandas as pd
import multiprocessing

from _0_data_downloading.utils import (download_satellite_data, merge_and_save_the_download_report,
                                               remove_empty_folders, fill_the_sat_paths, plot_the_maps_in_the_folder)

from utils import (get_all_cases_to_process)

def Download_satellite_data(core_arguments, nb_of_cores_to_use, overwrite_existing_satellite_files, where_to_save_satellite_data) : 
    
    """
    Function to download satellite data from the ODATIS FTP server using rsync.
    """
            
    cases_to_process = get_all_cases_to_process(core_arguments)

    download_report = {}
    
    for i in range(cases_to_process.shape[0]) : 
                
        info = cases_to_process.iloc[i].copy()
        
        progress = f'{i} over {cases_to_process.shape[0]-1} ({info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable} / {info.Temporal_resolution})'
        
        print(progress)
                            
        satellite_data = download_satellite_data(info, 
                                                 core_arguments['start_day'], 
                                                 core_arguments['end_day'], 
                                                 where_to_save_satellite_data, 
                                                 nb_of_cores_to_use,
                                                 overwrite_existing_satellite_files) 
        
        if satellite_data.to_process == False : 
            download_report[progress] = satellite_data.download_report    
            continue
        
        satellite_data.download_missing_satellite_data()
        
        download_report[progress] = satellite_data.download_report
                
    remove_empty_folders(where_to_save_satellite_data)
    
    merge_and_save_the_download_report(download_report, where_to_save_satellite_data)
    
def Plot_and_Save_the_map(core_arguments,
                          nb_of_cores_to_use,
                          where_are_saved_satellite_data,
                          start_day_of_maps_to_plot,
                          end_day_of_maps_to_plot) : 
        
    cases_to_process = get_all_cases_to_process(core_arguments)

    with multiprocessing.Pool(nb_of_cores_to_use) as pool:

        for i in range(cases_to_process.shape[0]) : 
                    
            info = cases_to_process.iloc[i].copy()
            
            init = download_satellite_data(info, start_day_of_maps_to_plot, end_day_of_maps_to_plot, 
                                           where_are_saved_satellite_data, nb_of_cores_to_use) 
    
            paths_to_sat_data = fill_the_sat_paths(info, init.destination_path_to_fill, 
                                                   local_path = True, 
                                                   dates = pd.date_range(start=init.start_day, end=init.end_day, freq="D"))
            
            pool.map(plot_the_maps_in_the_folder, paths_to_sat_data)
