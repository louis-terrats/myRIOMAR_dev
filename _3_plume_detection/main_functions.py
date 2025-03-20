import os, pickle, glob, multiprocessing, imageio, re

import numpy as np
import geopandas as gpd
import pandas as pd
import matplotlib as mpl
# import rpy2.robjects as robjects


from myRIOMAR.utils import (exit_program, expand_grid, align_bathymetry_to_resolution, 
                            unique_years_between_two_dates, store_arguments, load_shapefile_data,
                            path_to_fill_to_where_to_save_satellite_files,
                            fill_the_sat_paths)

from myRIOMAR._3_plume_detection.utils import (main_process, define_parameters, reduce_resolution, 
                                            create_polygon_mask, preprocess_annual_dataset_and_compute_land_mask,
                                            get_all_possibilities_for_plume_detection)


def apply_plume_mask(arguments) :
    
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
    Time_resolution : str
        Temporal resolution of the data ('DAILY', 'MONTHLY', etc.).

    Returns
    -------
    None
        The results are saved as files in the specified `working_directory`.
    """
    
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
    
    france_shapefile = load_shapefile_data()
    
    cases_to_process = get_all_possibilities_for_plume_detection(Zones, Data_sources, Sensor_names, 
                                                                 Atmospheric_corrections, Years, Temporal_resolution)
    
    for i in range(cases_to_process.shape[0]) : 
                
        info = cases_to_process.iloc[i].copy()
        info['Satellite_variable'] = 'SPM'
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.Satellite_sensor} / {info.atmospheric_correction} / {info.Year} / {info.Time_resolution})')

        # Check if the data source is from SEXTANT or similar databases.
        is_SEXTANT_file = np.isin(info.Data_source, ['SEXTANT', 'SAMUEL'])
            
        # Retrieve specific parameters based on the selected zone.
        parameters = define_parameters(info.Zone)
                    
        # Build the file pattern to locate the satellite data files.       
        file_names_pattern = fill_the_sat_paths(info, 
                                               path_to_fill_to_where_to_save_satellite_files(working_directory + 'RESULTS/' + info.Zone).replace('[TIME_FREQUENCY]', ''),
                                               local_path = True).replace('/*/*/*', f'MAPS/{info.Time_resolution}/{info.Year}/*.pkl')
                
        # Check if the directory for the year's data exists; if not, skip processing.
        if not os.path.exists( os.path.dirname(file_names_pattern) ) : 
            print(f"Missing satellite data here : {file_names_pattern}")
            print("Skip to the next iterate")
            continue
                    
        # Find all files matching the specified pattern.
        file_names = glob.glob(file_names_pattern)
        
        # Open the first file to extract the dataset structure.
        with open(file_names[0], 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data'] if info.Time_resolution != 'DAILY' else pickle.load(f)['map_data']                     
                
        # Reduce the spatial resolution of the dataset to match the new resolution
        ds_reduced = ( reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                          if parameters['lat_new_resolution'] is not None 
                          else ds )
                
        # Align bathymetry data to the dataset resolution.
        bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced, f'{working_directory}/RESULTS/{info.Zone}/Bathy_data.pkl')
            
        # Create a mask that delineate the searching area.
        inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)
        
        # Create an annual dataset used for validating cloud-free areas, and generate a land mask
        (map_wo_clouds, land_mask) = preprocess_annual_dataset_and_compute_land_mask( (re.sub(r'/[0-9]{2,}', '', file_names_pattern) 
                                                                                        .replace(info.Time_resolution, 'MULTIYEAR')
                                                                                        .replace('*.pkl', 'Averaged_over_multi-years.pkl')),
                                                                                        parameters)
        
        # Process each file in parallel
        with multiprocessing.Pool(nb_of_cores_to_use) as pool:

            results = pool.starmap(main_process, 
                        [(    file_name, 
                              file_names_pattern, 
                              parameters,
                              bathymetry_data_aligned_to_reduced_map,
                              is_SEXTANT_file,
                              france_shapefile, 
                              map_wo_clouds,
                              land_mask,
                              inside_polygon_mask) for file_name in file_names ])
        
        # Create a DataFrame from the results, sort by date, and save it to a CSV
        statistics = pd.DataFrame([x for x in results if x is not None]).sort_values('date').reset_index(drop = True)
        statistics.to_csv(f'{os.path.dirname(file_names[0]).replace("MAPS", "PLUME_DETECTION")}/Results.csv', index=False)
         
        # Create a GIF from the saved maps by combining all PNG images
        saved_maps = sorted( glob.glob(f'{os.path.dirname(file_names[0]).replace("MAPS", "PLUME_DETECTION")}/MAPS/*.png') )
        with imageio.get_writer(f'{os.path.dirname(file_names[0]).replace("MAPS", "PLUME_DETECTION")}/GIF.gif', mode='I', fps= 1) as writer:
            for figure_file in saved_maps :
                image = imageio.imread(figure_file)
                writer.append_data(image)
        
        # # For debugging
        # for file_name in file_names : 
        #     print(file_name)
        #     main_process( file_name, 
        #           file_names_pattern, 
        #           parameters,
        #           bathymetry_data_aligned_to_reduced_map,
        #           is_SEXTANT_file,
        #           france_shapefile, 
        #           map_wo_clouds,
        #           land_mask,
        #           inside_polygon_mask)
            
    # global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates().reset_index(drop = True)
        
    # for i in range(global_cases_to_process.shape[0]) : 
        
    #     info = global_cases_to_process.iloc[i].copy()
        
    #     plot_time_series_of_plume_areas(working_directory = working_directory, 
    #                                     Zone = info.Zone, 
    #                                     Data_source = info.Data_source, 
    #                                     Satellite_sensor = info.Satellite_sensor, 
    #                                     atmospheric_correction = info.atmospheric_correction, 
    #                                     Time_resolution = info.Time_resolution, 
    #                                     Years = np.arange(1998,2025))
    
    
def plot_time_series_of_plume_areas(working_directory, Zones, Data_sources, Satellite_sensors, Atmospheric_corrections, Time_resolutions, Years):
    
    """
    Calls the R function `plot_time_series_of_plume_area_and_thresholds` from Python.

    Args:
        working_directory (str): Working directory path.
        Zone (str): Zone name.
        Data_source (str): Data source name.
        Satellite_sensor (list): List of satellite sensors.
        atmospheric_correction (str): Atmospheric correction type.
        Time_resolution (str): Time resolution.
        Years (list): List of years.
    """
    
    # Source the R script
    robjects.r['source']("myRIOMAR/plume_detection/utils.R")

    r_function = robjects.r['plot_time_series_of_plume_area_and_thresholds']

    # Call the R function
    r_function(
        working_directory = robjects.StrVector([working_directory]),
        Zone = robjects.StrVector(Zones),
        Data_source = robjects.StrVector(Data_sources),
        Satellite_sensor = robjects.StrVector(Satellite_sensors),
        atmospheric_correction = robjects.StrVector(Atmospheric_corrections),
        Time_resolution = robjects.StrVector(Time_resolutions),
        Years = robjects.IntVector(list(Years))
    )
