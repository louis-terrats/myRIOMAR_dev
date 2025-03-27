#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:31:44 2025

@author: terrats
"""

import sys, pickle, os, bathyreq, glob, datetime, importlib.resources, tempfile, shutil
from itertools import product, chain
import pandas as pd
import xarray as xr
import numpy as np
import geopandas as gpd
import re
from functools import reduce
from collections.abc import Mapping, Iterable


def exit_program():
    print("Exiting the program...")
    sys.exit(0)
    
def load_file(file_name):
    
    with open(file_name, 'rb') as f:
        return pickle.load(f)
    
def expand_grid(**kwargs):
    
    """
    Create a DataFrame from the Cartesian product of input arrays.

    Parameters
    ----------
    **kwargs : dict
        Keyword arguments where keys are column names and values are arrays.

    Returns
    -------
    pandas.DataFrame
        DataFrame representing the Cartesian product of input arrays.
    """
    
    # Compute the Cartesian product of input values.
    rows = product(*kwargs.values())
    return pd.DataFrame(rows, columns=kwargs.keys())


def load_bathymetric_data(path_to_bathy_data, min_lon, max_lon, min_lat, max_lat) : 
    
    # If the bathymetric data doesn't exist, request it and save it
    if os.path.exists( path_to_bathy_data ) == False : 
        
        req = bathyreq.BathyRequest() # Create a bathymetric data request
        data, lonvec, latvec = req.get_area(longitude=[min_lon, max_lon], 
                                            latitude=[min_lat, max_lat])
        bathymetric_data = xr.DataArray(data, coords=[latvec[::-1], lonvec], dims=['lat', 'lon']) # Create a data array for bathymetry
        
        # Save the bathymetric data for future use
        with open(path_to_bathy_data, 'wb') as f:
            pickle.dump(bathymetric_data, f)
        
        # bathymetric_data.plot()
    else : 
        
        # Load the pre-saved bathymetric data
        with open(path_to_bathy_data, 'rb') as f:
            bathymetric_data = pickle.load(f)
            
    return bathymetric_data

def align_bathymetry_to_resolution(dataset, path_to_bathy_data) : 
               
    """
    Align bathymetric data to the resolution of the input dataset.

    Parameters
    ----------
    dataset : xarray.DataArray
        The input dataset to which the bathymetry should be aligned.
    parameters : dict
        Configuration parameters for plume detection.
    path_to_bathy_data : str
        Path to the raw bathymetric map of the Zone (e.g. f'{work_dir}/RESULTS/{Zone}/Bathy_data.pkl').
        

    Returns
    -------
    xarray.DataArray
        The bathymetric data aligned to the input dataset's resolution.
    """
        
    bathymetric_data = load_bathymetric_data(path_to_bathy_data, 
                                             min_lon = np.min(dataset.lon)-1, max_lon = np.max(dataset.lon)+1, 
                                             min_lat = np.min(dataset.lat)-1, max_lat = np.max(dataset.lat)+1)    
               
    # Align the bathymetric data to the reduced resolution dataset
    bathymetry_data_aligned_to_reduced_map = bathymetric_data.interp(lat= dataset.lat, lon= dataset.lon)
    
    return bathymetry_data_aligned_to_reduced_map


def degrees_to_km(lat_deg, lon_deg, latitude):
    """
    Convert distances in degrees to kilometers.
    
    Parameters:
    - lat_deg: Distance in degrees of latitude
    - lon_deg: Distance in degrees of longitude
    - latitude: Latitude where the conversion is needed
    
    Returns:
    - lat_km: Distance in kilometers for latitude
    - lon_km: Distance in kilometers for longitude
    """
    # Conversion factor for latitude
    lat_km = lat_deg * 111.32
    
    # Conversion factor for longitude, adjusted by latitude
    lon_km = lon_deg * 111.32 * np.cos(np.radians(latitude))
    
    return lat_km, lon_km


def km_to_degrees(lat_km, lon_km, latitude):
    """
    Convert a distance in meters to degrees of latitude and longitude.

    Args:
        meters (float): Distance in meters.
        latitude (float): Latitude in degrees where the conversion is applied.

    Returns:
        (float, float): (Latitude degrees, Longitude degrees)
    """
    # 1 degree latitude ≈ 111.32 km (constant)
    lat_deg = lat_km / 111.320  

    # 1 degree longitude ≈ 111.32 * cos(latitude) km
    lon_deg = lon_km / (111.320 * np.cos(np.radians(latitude)))

    return lat_deg, lon_deg

def generic_variable_names() : 
    
    return ['CHLA', 'SPM', 'SST']

def find_sat_data_files(info, path_to_sat_data) : 

    if isinstance(path_to_sat_data, str) : 
        path_to_sat_data = [path_to_sat_data]

    if ('Year' in info) and isinstance(info['Year'], str) and (info['Year'] == 'MULTIYEAR') : 
        file_pattern = '/*.pkl'
    elif np.isin(info.Satellite_variable, generic_variable_names()) : 
        file_pattern = '/*.nc'
    else : 
        file_pattern = f'/*{info.Satellite_variable}*.nc'    

    map_files = []
    
    if '[YEAR]' in path_to_sat_data[0] : 
        
        for Year in info['Year'] :    
    
            path_to_files = (path_to_sat_data 
                                 .replace('[YEAR]', str(Year))
                                 .replace('[MONTH]', '*')
                                 .replace('[DAY]', '*'))
                              
            files = glob.glob(path_to_sat_data + file_pattern)
            
            map_files.extend( files )
            
    else : 
        
        map_files = list(chain.from_iterable(glob.glob(path + file_pattern) for path in path_to_sat_data))
        
    return map_files

def store_arguments(arguments, locally = False, globally = False, return_arguments = False):
    
    arguments_to_return = []
    for key, value in arguments.items():
        if locally : 
            locals()[key] = value
            continue
        if globally : 
            globals()[key] = value
            continue
        if return_arguments : 
            arguments_to_return.append(value)
        
    if return_arguments : 
        return arguments_to_return
            
    # print(Data_sources)  # Works inside this function
    
def path_to_fill_to_where_to_save_satellite_files(where_to_save_files) : 
    
    path = f'{where_to_save_files}/[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]/[YEAR]/[MONTH]/[DAY]'

    return path

def fill_the_sat_paths(info, path_to_fill, local_path = False, dates = []) : 
    
    if len(dates) > 0 and isinstance(dates[0], str) : 
        dates = pd.to_datetime(dates)
    
    path_to_fill = ( path_to_fill  
                        .replace('[DATA_SOURCE]', info.Data_source)
                        .replace('[ATMOSPHERIC_CORRECTION]', info.atmospheric_correction)
                        .replace('[SENSOR]', info.sensor_name if 'sensor_name' in info else info.Satellite_sensor) )
                      
    if 'Temporal_resolution' in info.keys() or 'Temporal_resolution' in info.keys() : 
        
        Temporal_resolution = info.Temporal_resolution if 'Temporal_resolution' in info.keys() else info.Temporal_resolution
        if local_path == False : 
            Temporal_resolution = (Temporal_resolution
                                    .replace('DAILY', 'day')
                                    .replace('MONTHLY', 'month')
                                    .replace('WEEKLY', '8-day'))        
            
    elif isinstance(info.Year, str) and (info.Year == 'MULTIYEAR') :
        
        Temporal_resolution = 'ANNUAL'
        
    else :
        
        Temporal_resolution = 'DAILY'
            
    path_to_fill = path_to_fill.replace('[TIME_FREQUENCY]', Temporal_resolution)
        
    if local_path : 
        
        Folder_name_for_the_variable = ('CHLA' if 'CHL' in info.Satellite_variable.upper()
                                        else 'SPM' if 'SPM' in info.Satellite_variable.upper()
                                        else 'SST' if 'SST' in info.Satellite_variable.upper()
                                        else info.Satellite_variable)
        
        path_to_fill = (path_to_fill.replace('[PARAMETER]', Folder_name_for_the_variable))
        
    else : 
        
        path_to_fill = (path_to_fill.replace('[PARAMETER]', info.Satellite_variable_name_on_remote_folder))

    if len(dates) > 0 : 
        
        paths_to_sat_files = [ ( path_to_fill
                                      .replace('[YEAR]', str(date.year))
                                      .replace('[MONTH]', str(date.month).zfill(2))
                                      .replace('[DAY]', str(date.day).zfill(2))
                                      .replace("[DOY]", date.strftime("%j")) ) 
                                for date in dates ]
        
    else : 
        
        paths_to_sat_files = (path_to_fill
                                .replace('[YEAR]', '*')
                                .replace('[MONTH]', '*')
                                .replace('[DAY]', '*')
                                .replace("[DOY]", '*')) 
    
    return paths_to_sat_files

def get_all_cases_to_process(core_arguments) : 
        
    cases_to_process = expand_grid(Data_source = core_arguments['Data_sources'], 
                                    sensor_name = core_arguments['Sensor_names'], 
                                    atmospheric_correction = core_arguments['Atmospheric_corrections'],
                                    Satellite_variable = core_arguments['Satellite_variables'],
                                    Temporal_resolution = core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments else ['DAILY'])
    
    cases_to_process['atmospheric_correction'] = cases_to_process.apply(lambda row: 'Standard' 
                                                                        if np.isin(row['Data_source'], ['SEXTANT', 'EUMETSAT']) 
                                                                        else row['atmospheric_correction'], axis=1)
    cases_to_process['Satellite_variable'] = cases_to_process.apply(lambda row: 'SPM' 
                                                                        if row['Data_source'] == 'SEXTANT' and 'SPM' in row['Satellite_variable'] 
                                                                        else row['Satellite_variable'], axis=1)
    
    cases_to_process = cases_to_process.drop_duplicates().reset_index(drop = True)  
    
    return cases_to_process

def get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments) : 

    all_possibilities = expand_grid( Zone = core_arguments['Zones'],
                                    Data_source = core_arguments['Data_sources'], 
                                    sensor_name = core_arguments['Sensor_names'], 
                                    atmospheric_correction = core_arguments['Atmospheric_corrections'],
                                    Year = core_arguments['Years'],
                                    Satellite_variable = core_arguments['Satellite_variables'],
                                    Temporal_resolution = (core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments 
                                                       else core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments 
                                                       else '') )
    all_possibilities['atmospheric_correction'] = all_possibilities.apply(lambda row: 'Standard' 
                                                                        if row['Data_source'] == 'SEXTANT' 
                                                                        else row['atmospheric_correction'], axis=1)
    all_possibilities = all_possibilities.drop_duplicates()
    
    return all_possibilities

def create_arborescence(paths):
    arborescence = {}
    for path in paths:
        keys = path.split("/")
        current_level = arborescence
        for key in keys:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]  # Move deeper
    return arborescence

def return_the_parameter_name_based_on_file_name(file_name) : 
                        
    regular_expression = r'(?:(SPM-[G|R]|SPIM|suspended_matters|TSM_NN|CHL|CHL1|CHL-OC5|CHL-GONS|chlorophyll_a|POC|NRRS[0-9]*|RRS[0-9]*|DOC|CDOM|BBP|T-FNU|SST(?:-NIGHT|)))'
    match = re.search(regular_expression, file_name)
    
    if match : 
        
        return match.group(0) 
            
    else : 
        
        print('!!! Impossible to find the parameter name from the file name (see the function return_the_parameter_name_based_on_file_name) !!!')                        
                                                
def add_array_to_dict(dictionary, path, array):
    """
    Adds an array to a specific position in a nested dictionary.
    
    Parameters:
        dictionary (dict): The main dictionary to update.
        path (str): The path in the form "A/B/C/D/filename.nc".
        array (numpy.ndarray): The array to be stored.
    
    Returns:
        None (modifies the dictionary in place).
    """
    keys = path.split("/")
    filename = keys[-1]  # Extract the filename
    
    # Extract parameter name (e.g., SPM-G)
    param = return_the_parameter_name_based_on_file_name(filename)

    # Navigate the dictionary hierarchy
    current_level = dictionary
    for key in keys[:-1]:  # Exclude filename from navigation
        if key not in current_level:
            current_level[key] = {}
        current_level = current_level[key]

    # Add the array under the parameter name
    current_level[param] = array
    
def access_item_in_a_dictionnary(dictionary, path):
    
    keys = path.split("/")

    item = reduce(lambda dictionary, key: dictionary[key], keys[1:], dictionary)
    
    return item
        
    
def merge_dicts(dicts):
    """
    Merges multiple nested dictionaries into a single dictionary.
    
    Parameters:
        dicts (list): List of dictionaries to merge.
    
    Returns:
        dict: A merged dictionary.
    """
    def recursive_merge(d1, d2):
        """Recursively merges two dictionaries."""
        for key, value in d2.items():
            if key in d1 and isinstance(d1[key], dict) and isinstance(value, dict):
                recursive_merge(d1[key], value)  # Merge sub-dictionaries
            else:
                d1[key] = value  # Overwrite or add new keys
        return d1

    merged_dict = {}
    for d in dicts:
        merged_dict = recursive_merge(merged_dict, d)

    return merged_dict



def get_empty_paths(dictionary, prefix=""):
    
    paths = []
    
    if isinstance(dictionary, dict) and len(dictionary) > 0:  # If it's a non-empty dictionary
        for key, val in dictionary.items():
            paths.extend(get_empty_paths(val, f"{prefix}/{key}"))  # Recursive call
    elif len(dictionary) == 0:  # Exclude NaN values
        paths.append(prefix)  # Add the path if it's a valid value
        # paths.append( [f'{prefix}/{x}' for x in list(dictionary.keys())] )  # Add the path if it's a valid value
    
    return paths

def get_non_empty_paths(dictionary, prefix=""):
    
    paths = []
    
    if isinstance(dictionary, dict) and len(dictionary) > 0 :  # If it's a non-empty dictionary
        for key, val in dictionary.items():
            paths.extend(get_non_empty_paths(val, prefix = f"{prefix}/{key}"))  # Recursive call

    elif isinstance(dictionary, dict) == False:  # Exclude NaN values
        paths.append(prefix)  # Add the path if it's a valid value
    else :  # Exclude NaN values
        for key, val in dictionary.items():
            paths.append(f"{prefix}")  # Recursive call
            break
    
    # paths = np.unique( [x.replace('Sat_values', '').replace('Time', '') for x in paths] )
    
    return paths


def check_time_format(time_str):
    
    import re
    import numpy as np

    # Regular expression pattern for HH:MM:SS format
    time_pattern = r'^([0-1]?[0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9] UTC$'
    
    # Check if the time string matches the pattern
    if re.match(time_pattern, time_str):
        return time_str
    else:
        return np.nan
    
def extract_the_time_from_the_satellite_file(map_data) : 
        
    if 'image_reference_time' in map_data.attrs : # For SEXTANT products
        time = map_data._attrs['image_reference_time']
    elif 'DSD_entry_id' in map_data.attrs and 'L4' in map_data._attrs['DSD_entry_id'] : # For SEXTANT merged products
        time = ""
    elif 'start_time' in map_data.attrs :  # For ODATIS products
        time = pd.to_datetime(map_data.attrs['start_time']).strftime('%H:%M:%S UTC')
    elif 'time' in map_data.attrs :  # For ODATIS products    
        time = map_data.attrs['time']
        
    time = check_time_format(time)
    
    return time
       
def extract_dataframes_iterative(data):
    """Efficiently extract all DataFrames from a nested dictionary using an iterative approach."""
    stack = [data]  # Use a stack to avoid deep recursion

    while stack:
        current = stack.pop()

        if isinstance(current, pd.DataFrame):
            yield current  # Yield instead of appending to a list (memory-efficient)
        elif isinstance(current, Mapping):  # Check if it's a dictionary
            stack.extend(current.values())  # Add dictionary values to the stack
        elif isinstance(current, Iterable) and not isinstance(current, (str, bytes)):  
            stack.extend(current)  # Add list/tuple elements to the stack
            
def unique_years_between_two_dates(start_date: str, end_date: str):
    start_year = datetime.datetime.strptime(start_date, "%Y/%m/%d").year
    end_year = datetime.datetime.strptime(end_date, "%Y/%m/%d").year
    return list(range(start_year, end_year + 1))


def load_shapefile_data() : 
    
    france_shapefile = load_csv_files_in_the_package_folder(FRANCE_shapefile = True)
    
    return france_shapefile 
            
    # try : 
    #     france_shapefile = pygadm.Items(name="FRANCE", content_level=0)
    #     return france_shapefile
    # except Exception as e :
    #     print(f"The France shapefile can't be accessed through pygadm : {e}")
    #     print("The France shapefiles can be manually downloaded for free : e.g. https://gadm.org/download_country.html ")
 
    
def extract_and_format_date_from_path(path):
    match = re.search(r'/(\d{4})/(\d{2})/(\d{2})/', path)
    return ''.join(match.groups()) if match else None   

def load_csv_files_in_the_package_folder(SOMLIT = False, REPHY = False, FRANCE_shapefile = False):
    
    if SOMLIT : 
        with importlib.resources.open_text('myRIOMAR_dev._1_data_validation.INSITU_data.SOMLIT', 'Somlit.csv') as f:
            return (pd.read_csv(f, sep = ";", header = 2).iloc[1:]
                                .rename(columns = {'gpsLat*':'LATITUDE', 
                                                   'gpsLong*':'LONGITUDE',
                                                   'nomSite*':"Site"}))
        
    if REPHY : 
        with importlib.resources.open_binary('myRIOMAR_dev._1_data_validation.INSITU_data.REPHY', 'Table1_REPHY_hydro_RIOMAR.csv.gz') as f:
            return pd.read_csv(f, sep = ";", header = 0, encoding="ISO-8859-1", compression = {'method' : 'gzip'})
        
    if FRANCE_shapefile : 
        
        shp_folder = importlib.resources.files('myRIOMAR_dev._3_plume_detection.FRANCE_shapefile')  # Directly get the package folder path
    
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Extract all necessary shapefile components
            for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg']:
                shp_file = shp_folder / f'gadm41_FRA_0{ext}'
                if shp_file.exists():  # Ensure the file exists before copying
                    shutil.copy(shp_file, os.path.join(tmp_dir, f'gadm41_FRA_0{ext}'))
    
            # Read the shapefile from the temporary directory
            shapefile_path = os.path.join(tmp_dir, 'gadm41_FRA_0.shp')
            return gpd.read_file(shapefile_path)

    
    