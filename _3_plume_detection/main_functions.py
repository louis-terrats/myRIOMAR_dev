import os, pickle, glob, multiprocessing, imageio, re, gc
import pandas as pd
import rpy2.robjects as robjects


from _99_common.utils import (align_bathymetry_to_resolution, 
                            unique_years_between_two_dates, load_shapefile_data,
                            path_to_fill_to_where_to_save_satellite_files,
                            fill_the_sat_paths, get_all_cases_to_process_for_regional_maps_or_plumes_or_X11,
                            define_parameters)

from _3_plume_detection.utils import (reduce_resolution, create_polygon_mask, preprocess_annual_dataset_and_compute_land_mask)

import _3_plume_detection.utils


def apply_plume_mask(core_arguments, Zones, detect_plumes_on_which_temporal_resolution_data,
                    nb_of_cores_to_use, use_dynamic_threshold, 
                    where_are_saved_regional_maps, where_to_save_plume_results) :
    
    """
    Apply a plume mask to satellite data.

    This function processes satellite-derived data to detect and analyze plumes 
    (e.g., sediment plumes). It creates maps, applies filters, and extracts plume 
    regions based on specified parameters. Results are saved in various formats, 
    including CSV files and images.

    Parameters
    ----------
    working_directory : str
        Base directory for data storage and results.
    Zone : str
        The region of interest for plume detection.
    Data_Source : str
        Source of satellite data (e.g., 'SEXTANT').
    Satellite_sensor : str
        The satellite sensor used for the data.
    Atmospheric_correction : str
        Type of atmospheric correction applied to the data.
    Years : list of int
        List of years for which the data will be processed.
    Temporal_resolution : str
        Temporal resolution of the data ('DAILY', 'MONTHLY', etc.).

    Returns
    -------
    None
        The results are saved as files in the specified `working_directory`.
    """
    
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones,
                           'Satellite_variables': ['SPM'],
                           'Temporal_resolution' : ([detect_plumes_on_which_temporal_resolution_data] 
                                                    if isinstance(detect_plumes_on_which_temporal_resolution_data, str) 
                                                    else detect_plumes_on_which_temporal_resolution_data)})
    
    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)
        
    france_shapefile = load_shapefile_data()
        
    for i, info in cases_to_process.iterrows() : 
            
        # info = cases_to_process.iloc[i]
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Temporal_resolution})')
            
        # Retrieve specific parameters based on the selected zone.
        parameters = define_parameters(info.Zone)
                    
        # Build the file pattern to locate the satellite data files.       
        file_names_pattern = fill_the_sat_paths(info, 
                                               path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + "/" + info.Zone).replace('[TIME_FREQUENCY]', ''),
                                               local_path = True).replace('/*/*/*', f'MAPS/{info.Temporal_resolution}/{info.Year}/*.pkl')
                
        # Check if the directory for the year's data exists; if not, skip processing.
        if not os.path.exists( os.path.dirname(file_names_pattern) ) : 
            print(f"Missing satellite data here : {file_names_pattern}")
            print("Skip to the next iterate")
            continue
                    
        # Find all files matching the specified pattern.
        file_names = glob.glob(file_names_pattern)
        
        # Open the first file to extract the dataset structure.
        with open(file_names[0], 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data'] 
                
        # Reduce the spatial resolution of the dataset to match the new resolution
        ds_reduced = ( reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                          if parameters['lat_new_resolution'] is not None 
                          else ds )
                
        # Align bathymetry data to the dataset resolution.
        bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced, f'{where_are_saved_regional_maps}/{info.Zone}/Bathy_data.pkl')
            
        # Create a mask that delineate the searching area.
        inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)
        
        # Create an annual dataset used for validating cloud-free areas, and generate a land mask
        (map_wo_clouds, land_mask) = preprocess_annual_dataset_and_compute_land_mask( (re.sub(r'/[0-9]{2,}', '', file_names_pattern) 
                                                                                        .replace(info.Temporal_resolution, 'MULTIYEAR')
                                                                                        .replace('*.pkl', 'Averaged_over_multi-years.pkl')),
                                                                                        parameters)
        
        # Process each file in parallel
        with multiprocessing.Pool(nb_of_cores_to_use) as pool:

            results = pool.starmap(_3_plume_detection.utils.main_process, 
                        [(    file_name, 
                              file_names_pattern, 
                              parameters,
                              bathymetry_data_aligned_to_reduced_map,
                              france_shapefile, 
                              map_wo_clouds,
                              land_mask,
                              inside_polygon_mask,
                              where_to_save_plume_results,
                              where_are_saved_regional_maps,
                              use_dynamic_threshold) for file_name in file_names ])
        
        # Create a DataFrame from the results, sort by date, and save it to a CSV
        statistics = pd.DataFrame([x for x in results if x is not None]).sort_values('date').reset_index(drop = True)
        folder_name = os.path.dirname(file_names[0]).replace(where_are_saved_regional_maps, where_to_save_plume_results).replace("MAPS", "PLUME_DETECTION")
        os.makedirs(folder_name, exist_ok=True)
        statistics.to_csv(f'{folder_name}/Results.csv', index=False)
         
        # Create a GIF from the saved maps by combining all PNG images
        saved_maps = sorted( glob.glob(f'{folder_name}/MAPS/*.png') )
        with imageio.get_writer(f'{folder_name}/GIF.gif', mode='I', fps= 1) as writer:
            for figure_file in saved_maps :
                image = imageio.imread(figure_file)
                writer.append_data(image)
                
        del results, ds_reduced, bathymetry_data_aligned_to_reduced_map, inside_polygon_mask, map_wo_clouds, land_mask
        gc.collect()
        
        # # For debugging
        # for file_name in file_names : 
        #     print(file_name)
        #     main_process( file_name, 
        #           file_names_pattern, 
        #           parameters,
        #           bathymetry_data_aligned_to_reduced_map,
        #           france_shapefile, 
        #           map_wo_clouds,
        #           land_mask,
        #           inside_polygon_mask,
        #           where_to_save_plume_results,
        #           where_are_saved_regional_maps)
            
    # global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates().reset_index(drop = True)
        
    # for i in range(global_cases_to_process.shape[0]) : 
        
    #     info = global_cases_to_process.iloc[i].copy()
        
    #     plot_time_series_of_plume_areas(working_directory = working_directory, 
    #                                     Zone = info.Zone, 
    #                                     Data_source = info.Data_source, 
    #                                     Satellite_sensor = info.Satellite_sensor, 
    #                                     atmospheric_correction = info.atmospheric_correction, 
    #                                     Temporal_resolution = info.Temporal_resolution, 
    #                                     Years = np.arange(1998,2025))
    
    
def make_and_plot_time_series_of_plume_areas(core_arguments, Zones, nb_of_cores_to_use, on_which_temporal_resolution_the_plumes_have_been_detected, 
                                    where_are_saved_plume_results, where_to_save_plume_time_series):
    
    """
    Calls the R function `plot_time_series_of_plume_area_and_thresholds` from Python.

    Args:
        working_directory (str): Working directory path.
        Zone (str): Zone name.
        Data_source (str): Data source name.
        Satellite_sensor (list): List of satellite sensors.
        atmospheric_correction (str): Atmospheric correction type.
        Temporal_resolution (str): Time resolution.
        Years (list): List of years.
    """
    
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones,
                           'Satellite_variables': ['SPM'],
                           'Temporal_resolution' : ([on_which_temporal_resolution_the_plumes_have_been_detected] 
                                                    if isinstance(on_which_temporal_resolution_the_plumes_have_been_detected, str) 
                                                    else on_which_temporal_resolution_the_plumes_have_been_detected)})
    
    Plumes_per_zone = { Zone : list( define_parameters(Zone)['core_of_the_plumes'].keys() ) for Zone in Zones}   
    
    # Source the R script
    robjects.r['source']("myRIOMAR_dev/_3_plume_detection/utils.R")

    r_function = robjects.r['plot_time_series_of_plume_area_and_thresholds']

    # Call the R function
    r_function(
        where_are_saved_plume_results = robjects.StrVector([where_are_saved_plume_results]),
        where_to_save_plume_time_series = robjects.StrVector([where_to_save_plume_time_series]),
        Zone = robjects.StrVector(core_arguments['Zones']),
        Data_source = robjects.StrVector(core_arguments['Data_sources']),
        Satellite_sensor = robjects.StrVector(core_arguments['Sensor_names']),
        atmospheric_correction = robjects.StrVector(core_arguments['Atmospheric_corrections']),
        Temporal_resolution = robjects.StrVector(core_arguments['Temporal_resolution']),
        Years = robjects.IntVector(core_arguments['Years']),
        Plumes = robjects.ListVector(Plumes_per_zone),
        nb_of_cores_to_use = robjects.IntVector([nb_of_cores_to_use])
    )
