import subprocess, re, multiprocessing, time, bz2, os, gzip, glob, pyinterp, gc, shutil
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
from ftplib import FTP
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    from ..shared_utils import (find_sat_data_files, km_to_degrees, path_to_fill_to_where_to_save_satellite_files,
                                fill_the_sat_paths, extract_the_time_from_the_satellite_file)
except ImportError:
    from myRIOMAR_dev.CODES.shared_utils import (find_sat_data_files, km_to_degrees, path_to_fill_to_where_to_save_satellite_files,
                                fill_the_sat_paths, extract_the_time_from_the_satellite_file)


def get_connexion_parameters(info) : 

    all_ftp_connexion_parameters = {
        'SEXTANT' : {
                'user' :"",
                'password' : "",
                'host' : "ftp.ifremer.fr",
                'extension' : ( '.bz2' if info.sensor_name == 'merged' else '.gz' ),
                'remote_path_to_fill' : ( 'ifremer/cersat/products/gridded/ocean-color/atlantic/EUR-L4-[PARAMETER]-ATL-v01/[YEAR]/[DOY]/' 
                                             if info.sensor_name == 'merged' 
                                             else 'ifremer/cersat/products/gridded/ocean-color/atlantic/[SENSOR]/[YEAR]/')
            }, 
        'ODATIS' : {
                'user' :"ftp_odatis-cnes",
                'password' : "Acri%cs4",
                'host' : "ftp.acrist-services.com",
                'extension' : '.nc',
                'remote_path_to_fill' : '[ATMOSPHERIC_CORRECTION]/[SENSOR]/[TIME_FREQUENCY]/[YEAR]/[MONTH]/[DAY]/' 
            },
        'EUMETSAT' : {
                'user' :"QBZG6Py6TmJ8_0xZWyBSStADPica",
                'password' : "mCT8oYZDZjh0Yb9XElKbjzhOheoa",
                'host' : "",
                'extension' : '.nc',
                'remote_path_to_fill' : '' 
            },
        
        }
    
    return all_ftp_connexion_parameters[info.Data_source]

def download_files_from_ftp(self):
    """
    Download an entire folder from an FTP server.
    
    Parameters:
        ftp_host (str): FTP server address.
        ftp_user (str): FTP username.
        ftp_pass (str): FTP password.
        remote_folder (str): Path to the remote folder on the FTP server.
        local_folder (str): Local destination folder.

    Example:
        download_ftp_folder("ftp.example.com", "myuser", "mypassword", "/remote/folder/", "./local_folder/")
    """
    
    remote_folders = fill_the_sat_paths(self.info, 
                                        path_to_fill = get_connexion_parameters(self.info)['remote_path_to_fill'], 
                                        local_path = False,
                                        dates = self.dates_to_download)
        
    local_folders = fill_the_sat_paths(self.info, 
                                        path_to_fill = (self.destination_path_to_fill.replace('[PARAMETER]', 'ALL') 
                                                        if (self.info.Data_source == 'SEXTANT') and (self.info.sensor_name != 'merged')
                                                        else self.destination_path_to_fill),
                                        local_path = True,
                                        dates = self.dates_to_download) 
    
    ftp_connexion_parameters = get_connexion_parameters(self.info)
               
    try : 
        ftp = FTP(ftp_connexion_parameters['host'])
        ftp.login(ftp_connexion_parameters['user'], ftp_connexion_parameters['password'])
    except Exception as e:
        self.download_report.loc[self.download_report.index.isin(self.dates_to_download), 'Message'] = f"❌ FTP Error when connecting to FTP address: {e}"
        return None
        
    downloaded_files = 0
    for i, (remote_folder, local_folder, date_folder) in enumerate(zip(remote_folders, local_folders, self.dates_to_download.astype(str))) :  
                    
        # print( f'{i} over {len(remote_folders)-1}' )
        # remote_folder = remote_folders[0]
        # local_folder = local_folders[0]
        # date_folder = self.dates_to_download.astype(str)[0]
        
        if downloaded_files > 0 and downloaded_files % 30 == 0 : # Pause the execution for 30'' every 20 files downloaded.
            print('30-second run stop')
            time.sleep(30) 
        
        os.makedirs(local_folder, exist_ok=True)  # Ensure local directory exists
    
        try:
            ftp.cwd(remote_folder)  
        except Exception as e:
            if self.info.Temporal_resolution == 'DAILY' : 
                self.download_report.loc[self.download_report.index == date_folder, 'Message'] = f"❌ This folder does not exist {remote_folder}: {e}"
            elif self.info.Temporal_resolution == 'WEEKLY' : 
                self.download_report.loc[self.download_report.index == date_folder, 'Message'] = np.nan # f"❌ No weekly data for {remote_folder}: {e}"
            continue 
        
        try:
            file_names = ftp.nlst()  # List all files in the directory
        except Exception as e:
            self.download_report.loc[self.download_report.index == date_folder, 'Message'] = f"❌ Error listing files in {remote_folder}: {e}"
            ftp.cwd("~")
            continue  # Skip to the next folder if listing fails
            
        # Pre-filter files based on extension and satellite variable
        valid_files = [
            file for file in file_names
            if file.endswith(ftp_connexion_parameters['extension']) and 
                ( "".join(re.search(r"(\d{4})/(\d{2})/(\d{2})$", local_folder).groups()) in file if self.info.Data_source == 'SEXTANT' and self.info.sensor_name != "merged"  
                 else self.info.Satellite_variable_name_on_remote_folder in file )
        ]
        
        if len(valid_files) == 0 : 
            self.download_report.loc[self.download_report.index == date_folder, 'Message'] = f"❌ No satellite file in {remote_folder}"
            ftp.cwd("~")
            continue  # Skip to the next folder if listing fails
            
        for file_name in valid_files:
            
            local_file = os.path.join(local_folder, file_name)
            
            try : 
                with open(local_file, "wb") as f:
                    ftp.retrbinary("RETR " + file_name, f.write)
                self.download_report.loc[self.download_report.index == date_folder, 'Message'] += f" ✅ Downloaded"
                downloaded_files =+ 1
                
                if 'global_files_already_downloaded' in vars(self) :
                    self.global_files_already_downloaded = np.append(self.global_files_already_downloaded, local_file.replace('.gz', '').replace('.bz2', ''))
                                
            except Exception as e:
                
                self.download_report.loc[self.download_report.index == date_folder, 'Message'] += f" ❌ FTP Error for {file_name}: {e}"
                continue
            
            if file_name.endswith('.bz2') : 
                unzip_bz2_file(local_file)
            
            if file_name.endswith('.gz') : 
                unzip_gz_file(local_file)
                
        ftp.cwd("~")

    ftp.quit()
    
    return None
        
        

        
def convert_dates_to_file_path_format(dates, destination_path_to_fill) : 
    
    file_path_format = ""
            
    file_path_format = re.search(r"\[YEAR\].*", destination_path_to_fill ).group(0)

    file_path_format = [ (file_path_format.replace("[YEAR]", str(date.year))
                                           .replace("[MONTH]", str(date.month).zfill(2))
                                           .replace("[DAY]", str(date.day).zfill(2))
                                           .replace("[DOY]", str(date.day_of_year).zfill(3))) 
                         for date in dates ]
            
    return file_path_format 
                

def which_dates_occured_in_the_sat_file_names(dates_to_look_for, sat_file_names, check_for_non_occurrence_instead = False) : 

    # Convert files_in_the_time_range to a set for O(1) lookups
    files_set = set(sat_file_names)
    
    # Generate the list of dates in the required format
    date_patterns = set(dates_to_look_for.strftime('%Y/%m/%d'))  # Convert to set for fast membership checking
    
    # Check which dates are NOT in files_set    
    if check_for_non_occurrence_instead : 
        matches = {date for date in date_patterns if any(date in file for file in files_set) == False}
    else : 
        matches = {date for date in date_patterns if any(date in file for file in files_set)}
        
    matches = pd.to_datetime( list(matches) )
    
    return(matches)
        

def format_variable_name_with_server_names( info ) :
     
    if info.Data_source == 'ODATIS' :
    
        to_return = ('CHL-OC5' if info.Satellite_variable == 'CHLA' 
                     else info.Satellite_variable ) 
        
    if info.Data_source == 'SEXTANT' :
        
        to_return = ('CHL' if info.Satellite_variable == 'CHLA' 
                     else 'SPIM' if 'SPM' in info.Satellite_variable
                     else info.Satellite_variable ) 
         
    if info.Data_source == 'EUMETSAT' :
        
        to_return = ('chl' if 'CHLA' == info.Satellite_variable
                     else 'chl_oc4me' if 'CHLA_OC4' == info.Satellite_variable
                     else 'chl_nn' if 'CHLA_NN' == info.Satellite_variable
                     else 'tsm' if 'SPM' in info.Satellite_variable.upper()
                     else 'tsm_nn' if 'SPM_NN' == info.Satellite_variable
                     else info.Satellite_variable ) 
         
    return to_return


def unzip_bz2_file(local_file):
                
    # Extract the filename without .bz2 extension
    output_file = local_file.replace(".bz2", "")

    # Decompress the file
    with bz2.BZ2File(local_file, 'rb') as f_in, open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        
    os.remove(local_file)
    
def unzip_gz_file(local_file):
                
    # Extract the filename without .bz2 extension
    output_file = local_file.replace(".gz", "")

    # Decompress the file
    with gzip.open(local_file, 'rb') as f_in, open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        
    os.remove(local_file)
    
    
def merge_and_save_the_download_report(download_report, where_to_save_satellite_data) : 

    # Convert dictionary to a single DataFrame with dates as index
    download_report = pd.concat(download_report, axis=1)
    
    # Rename columns for clarity
    download_report.columns = download_report.columns.droplevel(1)  # Remove unnecessary sub-column level
    download_report = ( download_report
                       .rename(columns=lambda x: f"{x}_Status")
                       .reset_index()
                       .rename(columns={'index': 'Dates'})) # Rename columns
        
    os.makedirs(f"{where_to_save_satellite_data}/DOWNLOAD_REPORTS/", exist_ok=True)  # Ensure local directory exists
    report_name = f"{where_to_save_satellite_data}/DOWNLOAD_REPORTS/{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.csv"
    download_report.to_csv(report_name, index=False, encoding="utf-8")

    print(f'Download report has been saved here : {report_name}')
   
                    
def remove_empty_folders(root_folder):
    """
    Efficiently removes empty folders recursively, starting from the deepest level.
    
    Args:
        root_folder (str): Path to the root folder to clean.
    """
    for dirpath, _, filenames in os.walk(root_folder, topdown=False):
        try:
            os.rmdir(dirpath)  # Removes the folder if it's empty
            # print(f"🗑️ Removed empty folder: {dirpath}")
        except OSError:
            pass  # Ignore errors if the folder is not empty
    
    # shutil.rmtree(root_folder + "SEXTANT/ALL", ignore_errors=True)
            
            
            
def download_L2_maps_from_EUMDAC(info, where_to_save_satellite_data, min_lon, max_lon, min_lat, max_lat, dates_to_download) : 
        
    # =============================================================================
    #% Pre-requisite
    # =============================================================================
    
    # Create an account here : https://eoportal.eumetsat.int/
    
    # Generate a token here : https://api.eumetsat.int/api-key/
    
    # Find your collection of interest : https://api.eumetsat.int/data/browse/1.0.0/collections?format=html # Each collection is a bunch of products
    
    # Install eumdac : 
        ### 'pip install eumdac' (more details here : https://user.eumetsat.int/resources/user-guides/eumetsat-data-access-client-eumdac-guide)
    
    # Modify the function cli.py to prevent any prompt during the downloading : 
        ### Search for the eumdac folder (e.g. /home/terrats/anaconda3/lib/python3.11/site-packages/eumdac/)
        ### In the folder, open cli.py
        ### Comment these lines : "     # if products_count >= 10 and not args.yes:
                                            # user_in = input("Do you want to continue (Y/n)? ")
                                            # if user_in.lower() == "n":
                                                # return "
    
    # =============================================================================
    #% Documents
    # =============================================================================
    
    # Help of the command eumdac : https://user.eumetsat.int/resources/user-guides/eumetsat-data-access-client-eumdac-guide
    
    # =============================================================================
    #% Check the permission
    # =============================================================================
    
    # Before using the CLI for the first time, ensure that you’re registered.
    # Run this bash command : 'eumdac set-credentials ConsumerKey ConsumerSecret'. ConsumerKey and ConsumerSecret are the token identifiers generated here : https://api.eumetsat.int/api-key/
    
    # =============================================================================
    #% Configuration 
    # =============================================================================
    
    # Set credentials 
    # ConsumerKey = "QBZG6Py6TmJ8_0xZWyBSStADPica"
    # ConsumerSecret = "mCT8oYZDZjh0Yb9XElKbjzhOheoa"
    connexion_parameters = get_connexion_parameters(info)
    
    # Name of the collection
    if 'SST' in info.Satellite_variable : 
        collection_names = {"Processed" : {"Dates" : [pd.to_datetime('2017-07-05'), pd.Timestamp.now()],
                                           "Collection_name" : "EO:EUM:DAT:0412"}, 
                            "Re-Processed" : {"Dates" : [ pd.to_datetime('2016-04-18'), pd.to_datetime('2018-04-04') ],
                                              "Collection_name" : "EO:EUM:DAT:0582"}} # list of collections here : https://api.eumetsat.int/data/browse/1.0.0/collections?format=html
    else : 
        collection_names = {"Processed" : {"Dates" : [pd.to_datetime('2017-07-05'), pd.Timestamp.now()],
                                           "Collection_name" : "EO:EUM:DAT:0407"}, 
                            "Re-Processed" : {"Dates" : [ pd.to_datetime('2016-04-25'), pd.to_datetime('2021-04-28') ],
                                              "Collection_name" : "EO:EUM:DAT:0556"}} # list of collections here : https://api.eumetsat.int/data/browse/1.0.0/collections?format=html

    collection_names = adjust_collection_names_to_start_and_end_dates(collection_names, 
                                                                      np.min(dates_to_download), 
                                                                      pd.to_datetime(f'{np.max(dates_to_download)} 23:59:59'))
                
    # Name of the platform
    Sensor_name = info.sensor_name.replace('olcia', 'S3A').replace('olcib', 'S3B')
    
    # Parameters    
    # Variables_to_download = [ var.replace('CHLA', 'chl_nn').replace('SPM', 'tsm_nn') for var in info.Satellite_variable ]
    Variables_to_download = [ info.Satellite_variable.replace('CHLA', 'chl_nn').replace('SPM', 'tsm_nn') , "geo_coordinates" ]
    if 'chl_nn' in Variables_to_download : 
        Variables_to_download.extend(['chl_oc4'])        
        
    # Type
    if 'SST' in info.Satellite_variable : 
        Type = "WST" # Water Full Resolution
    else :
        Type = "WFR" # Water Full Resolution
    
    # =============================================================================
    #% Execute the Bash command
    # =============================================================================
    
    merged_command = []
    
    for collection_name in collection_names.values() :
        
        # print(collection_name['Dates'])
        command_to_download_data = f'eumdac download -c {collection_name["Collection_name"]} --yes --start {collection_name["Dates"][0].strftime("%Y-%m-%dT%H:%M")} --end {collection_name["Dates"][1].strftime("%Y-%m-%dT%H:%M")} --bbox {min_lon} {min_lat} {max_lon} {max_lat} --filename "{Sensor_name}*2*{Type}*" --entry ' + f'{" ".join([ f"*{var}*" for var in Variables_to_download])}'
            
        # command_to_activate_the_permissions_to_download = 'curl -k -d "grant_type=client_credentials" \
        #                                                     -H "Authorization: Basic UUJaRzZQeTZUbUo4XzB4Wld5QlNTdEFEUGljYTptQ1Q4b1laRFpqaDBZYjlYRWxLYmp6aE9oZW9h" \
        #                                                     https://api.eumetsat.int/token'
                                                            
        command_to_activate_the_permissions_to_download = f'eumdac set-credentials --yes {connexion_parameters["user"]} {connexion_parameters["password"]}'
                        
        where_to_save_data_of_the_sat = f'"{where_to_save_satellite_data}/EUMETSAT/L2/{info.sensor_name}/"'
        command_to_create_the_destination_folder_if_it_does_not_exist = f'mkdir -p {where_to_save_data_of_the_sat}'            
        command_to_set_the_destination_folder = f'-o {where_to_save_data_of_the_sat}'
                    
        command_to_execute = f'{command_to_create_the_destination_folder_if_it_does_not_exist} && {command_to_activate_the_permissions_to_download} && {command_to_download_data + " " + command_to_set_the_destination_folder}'
        
        merged_command.append(command_to_execute)   
              
    the_merged_command = " && ".join(merged_command)
    
    for the_command in the_merged_command.split(' && ') :
        # print(the_command)
        subprocess.run(the_command.replace('"', '').split(' '), stdout=subprocess.DEVNULL)
                
    
    
def adjust_collection_names_to_start_and_end_dates(collection_names, start_date, end_date) : 
    
    if end_date <= collection_names["Re-Processed"]['Dates'][1] :
        
        del collection_names['Processed']
        collection_names['Re-Processed']['Dates'] = [start_date, end_date]
        
    elif start_date >= collection_names["Re-Processed"]['Dates'][0] :
        
        del collection_names['Re-Processed']
        collection_names['Processed']['Dates'] = [start_date, end_date]
        
    else : 
        
        collection_names['Re-Processed']['Dates'][0] = start_date
        collection_names['Processed']['Dates'][0] = collection_names['Re-Processed']['Dates'][1] + datetime.timedelta(days=1)
        collection_names['Processed']['Dates'][1] = end_date
        
    return collection_names


def regrid_L2_sat_data_to_regular_grid(L2_files,
                                       L2_coordinate_ds,
                                       new_grid) : 
    
    regridded_map_files = {}
    
    time_value = extract_the_time_from_the_satellite_file(L2_coordinate_ds)
    
    for L2_map_file in [ x for x in L2_files if x.endswith('geo_coordinates.nc') == False ] : 
                
        L2_map_ds = xr.open_dataset( L2_map_file )

        if 'sea_surface_temperature' in L2_map_ds.data_vars :
            var_name = 'sea_surface_temperature'
        else :
            var_name = list(L2_map_ds.data_vars)[0]

        new_grid.clear()
        new_grid.push(L2_coordinate_ds['longitude'].values, 
                     L2_coordinate_ds['latitude'].values, 
                     L2_map_ds[var_name].values, 
                     simple=False)
        nearest = new_grid.variable('mean')
        
        L2_map_ds.close()

        lon, lat = np.meshgrid(new_grid.x, new_grid.y, indexing='ij')
        ds_regular = xr.DataArray(data = nearest,
                                  dims=('lon', 'lat'),
                                  coords={'lon': lon[:,0], 'lat': lat[0]},
                                  attrs={'time': time_value}).T

        if 'sea_surface_temperature' in L2_map_ds.data_vars :
            ds_regular.values = ds_regular.values - 273.15 # Conversion from Kelvin to Celsius

        # ds_regular.plot()
    
        regridded_map_files[var_name] = ds_regular
        
    return regridded_map_files
    

def regrid_and_save_EUMDAC_maps_of_one_day(info, where_to_save_satellite_data, destination_path_to_fill, 
                                           the_day, new_grid, min_lon, max_lon, min_lat, max_lat) :
              
    name_pattern_of_L2_map_folders = f'{where_to_save_satellite_data}/EUMETSAT/L2/{info.sensor_name}/*____{the_day.strftime("%Y%m%dT")}*'
    
    L2_map_folders = glob.glob(name_pattern_of_L2_map_folders)

    if len(L2_map_folders) == 0 : 
        downloading_message = f" ❌ No satellite data for {the_day}"
        return {the_day : downloading_message}

    all_regridded_maps = {}
    downloading_message = ""
    for L2_map_folder in L2_map_folders :     

        L2_files = glob.glob(L2_map_folder + "/*.nc")            

        if len(L2_files) == 0 : 
            downloading_message = downloading_message + f" ❌ No satellite data for the scene {L2_map_folder.rsplit('/', 1)[-1]}"
            continue
        
        try : 
            
            if 'SST' in info.Satellite_variable :
                L2_coordinate_ds = xr.open_dataset( L2_files[0] ).rename({'lat': 'latitude', 'lon': 'longitude'}) 
            else :
                L2_coordinate_ds = xr.open_dataset( [ x for x in L2_files if x.endswith('/geo_coordinates.nc') ][0] ) 
        
            regridded_maps = regrid_L2_sat_data_to_regular_grid(L2_files = L2_files, 
                                                                L2_coordinate_ds = L2_coordinate_ds,
                                                                new_grid = new_grid)
                        
            all_regridded_maps[L2_map_folder] = regridded_maps
            
            L2_coordinate_ds.close()
            shutil.rmtree(L2_map_folder)
            
        except Exception as e:
            downloading_message = downloading_message + f" ❌ Error regridding satellite data for the scene {L2_map_folder.rsplit('/', 1)[-1]} : {e}"
            continue
            
    for parameter in np.unique( np.concatenate([ list(x.keys()) for _ , x in all_regridded_maps.items() ]) ) : 
        
        regridded_maps_without_the_parameter = [key for key, da in all_regridded_maps.items() if parameter not in da.keys()]
        for key_of_the_map in regridded_maps_without_the_parameter :
            downloading_message = downloading_message + f" ❌ No {parameter} data for the scene {key_of_the_map.rsplit('/', 1)[-1]}"
        
        where_to_save_regridded_map = fill_the_sat_paths(info = info, path_to_fill = destination_path_to_fill, 
                                                         local_path = True, dates = [the_day])[0]
        
        try : 
            
            regridded_map = ( xr.concat( [x[parameter] for _, x in all_regridded_maps.items()] , dim = 'stacked' )
                                 .mean(dim="stacked", skipna=True) )
            regridded_map.name = parameter
            
            list_of_times = [x[parameter].attrs['time'] for _, x in all_regridded_maps.items()]
            mean_time = pd.to_datetime( [t for t in list_of_times] ).mean().time().strftime('%H:%M:%S UTC')
            
            regridded_map = regridded_map.to_dataset().assign_attrs({'time': mean_time})
            
            os.makedirs(where_to_save_regridded_map, exist_ok=True)
            regridded_map.to_netcdf(f"{where_to_save_regridded_map}/{parameter}.nc",
                                    encoding={regridded_map[parameter].name: {"zlib": True, "complevel": 4, "dtype": "float32"}})
            
            downloading_message = downloading_message + f" ✅ Downloaded and Regridded: {parameter}.nc"
            
            regridded_map.close()
            del regridded_map
            gc.collect()
            
        except Exception as e:
            downloading_message = downloading_message + f" ❌ Error creating the satellite map of {parameter} : {e}"
            continue
            
    del all_regridded_maps
    gc.collect()
        
    return {the_day : downloading_message}


def download_EUMDAC_maps(self, min_lon, max_lon, min_lat, max_lat, new_map_resolution) : 
            
    # for dates_to_download in self.splitted_dates_to_download : 
    #     download_L2_maps_from_EUMDAC(self.info, self.where_to_save_satellite_data,
    #                                     min_lon, max_lon, 
    #                                     min_lat, max_lat, dates_to_download) 
    
    # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
    with multiprocessing.Pool(self.nb_of_cores_to_use) as pool:

        download_reports = pool.starmap(download_L2_maps_from_EUMDAC, [(self.info, self.where_to_save_satellite_data,
                                                       min_lon, max_lon, 
                                                       min_lat, max_lat, dates_to_download) 
                                                   for dates_to_download in self.splitted_dates_to_download ])
        
        new_grid = define_the_new_grid_for_map_regridding(resolution_in_m_of_the_new_grid = new_map_resolution, 
                                                          mean_latitude_of_the_area = 46.23, 
                                                          min_lon_new_grid = min_lon, 
                                                          max_lon_new_grid = max_lon, 
                                                          min_lat_new_grid = min_lat, 
                                                          max_lat_new_grid = max_lat) 
            
        # Use multiprocessing to process each week
        download_reports = pool.starmap(regrid_and_save_EUMDAC_maps_of_one_day, [(self.info, 
                                                                                  self.where_to_save_satellite_data, 
                                                                                  self.destination_path_to_fill, 
                                                       the_day, new_grid, 
                                                       min_lon, max_lon, 
                                                       min_lat, max_lat) 
                                                   for the_day in self.dates_to_download ])
        
        # For debugging
        # download_reports = []
        # for the_day in pd.date_range(start=self.start_day, end=self.end_day, freq="D") :
            
        #     download_report = regrid_and_save_EUMDAC_maps_of_one_day(self, the_day, new_grid, 
        #                                            min_lon, max_lon, 
        #                                            min_lat, max_lat) 
        #     download_reports.append(download_report)
        
    for date, download_message in {k.strftime('%Y-%m-%d'): v for d in download_reports for k, v in d.items()}.items() : 
        
        self.download_report.loc[self.download_report.index == date, 'Message'] = download_message

            

def define_the_new_grid_for_map_regridding(resolution_in_m_of_the_new_grid, 
                                           mean_latitude_of_the_area,
                                           min_lon_new_grid,
                                           max_lon_new_grid,
                                           min_lat_new_grid,
                                           max_lat_new_grid) : 
    
    diff_lat_in_degrees = km_to_degrees(lat_km = resolution_in_m_of_the_new_grid / 1000, lon_km = 0, latitude = mean_latitude_of_the_area)[0] # latitude of France
    diff_lon_in_degrees = km_to_degrees(lat_km = 0, lon_km = resolution_in_m_of_the_new_grid / 1000, latitude = mean_latitude_of_the_area)[1] # latitude of France

    # Compute bin edges and centers
    new_lons_edges = np.arange(min_lon_new_grid, max_lon_new_grid + diff_lon_in_degrees, diff_lon_in_degrees)
    new_lats_edges = np.arange(min_lat_new_grid, max_lat_new_grid + diff_lat_in_degrees, diff_lat_in_degrees)

    # Compute bin centers
    new_lons = (new_lons_edges[:-1] + new_lons_edges[1:]) / 2
    new_lats = (new_lats_edges[:-1] + new_lats_edges[1:]) / 2

    binning = pyinterp.Binning2D(
        pyinterp.Axis(new_lons, is_circle=False),
        pyinterp.Axis(new_lats))
    
    return binning
    
def test_if_satellite_data_DOES_NOT_exist_in_the_Data_source(info) : 
                
    test = (
                
        (info.Data_source == 'ODATIS' and 
         ((info.atmospheric_correction == 'Standard') or 
          (np.isin(info.sensor_name, ['modis', 'meris', 'olcia', 'olcib']) == False) or
          (info.atmospheric_correction == 'acolite' and info.sensor_name == 'modis') or
          ('SST' in info.Satellite_variable and info.sensor_name != 'modis') or
          (info.atmospheric_correction == 'nirswir' and info.sensor_name != 'modis') )) or
                        
        (info.Data_source == 'SEXTANT' and 
         ((info.atmospheric_correction != 'Standard') or 
          (info.Temporal_resolution != 'DAILY') or
          (np.isin(info.sensor_name, ['modis', 'meris', 'seawifs', 'viirs', 'merged']) == False) or
          all(param not in info.Satellite_variable for param in ['CHL', 'SPM']) )) or
                        
        (info.Data_source == 'EUMETSAT' and 
         ((info.atmospheric_correction != 'Standard') or 
          (info.Temporal_resolution != 'DAILY') or
          (np.isin(info.sensor_name, ['olcia', 'olcib']) == False) or
          all(param not in info.Satellite_variable for param in ['CHL', 'SPM', 'SST']) ))
        
    )
    
    return test


def split_consecutive_dates(dates):
    """
    Splits a DatetimeIndex into groups of consecutive dates.

    Parameters:
        dates (pd.DatetimeIndex): The dates to split.

    Returns:
        list of pd.DatetimeIndex: List of consecutive date groups.
    """
    if not isinstance(dates, pd.DatetimeIndex):
        dates = pd.to_datetime(dates)  # Ensure it's a DatetimeIndex

    dates = dates.sort_values()  # Sort dates

    # Find where gaps occur (difference is not 1 day)
    breaks = np.where(np.diff(dates) > pd.Timedelta(days=1))[0]
    
    splitted_dates = [dates[i:j] for i, j in zip(np.r_[0, breaks + 1], np.r_[breaks + 1, len(dates)])]

    # Efficiently split using numpy
    return splitted_dates


def plot_the_maps_in_the_folder(path) : 
    
    # mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM memory !)
    mpl.use('agg') # Prevent showing plot in the Plot panel (this saves RAM memory)
    
    if os.path.isdir(path) == False : 
        return
    
    nc_files = [f for f in os.listdir(path) if f.endswith(".nc")]
    
    for nc_file in nc_files:
        
        file_path = os.path.join(path, nc_file)
        
        with xr.open_dataset(file_path, chunks={}) as ds:
        
            var_name = list(ds.data_vars)[0]
            
            min_limit = ds.quantile(0.1, skipna=True)[var_name].values
            max_limit = ds.quantile(0.9, skipna=True)[var_name].values
            
            plt.figure(figsize=(12, 9))
            ds[var_name].plot(vmin = min_limit, vmax = max_limit)
            plt.title(var_name)
            
            output_file = os.path.join(path, f"{nc_file.replace('.nc', '.png')}")
            plt.savefig(output_file, dpi=150)
            plt.close()
        
            print(f"Plot saved here : {output_file}")
            
        gc.collect()


def find_files_in_the_time_range(info, destination_path_to_fill, all_dates) : 
                
    path_to_sat_data = fill_the_sat_paths(info = info, path_to_fill = destination_path_to_fill, local_path = True, dates = all_dates)
    
    all_files_already_downloaded = find_sat_data_files(info, path_to_sat_data = path_to_sat_data)
    
    pattern_to_match = f"(^.*/{'.*$)|(^.*/'.join( convert_dates_to_file_path_format(all_dates, destination_path_to_fill) )}.*$)"

    files_in_the_time_range = np.array(all_files_already_downloaded)[  np.where( pd.Series(all_files_already_downloaded).str.match(pattern_to_match) )[0] ]
    
    return files_in_the_time_range



def extract_values_from_the_global_file(self) : 
    
    where_to_save_values = fill_the_sat_paths(self.info, 
                                        path_to_fill = self.destination_path_to_fill,
                                        local_path = True,
                                        dates = self.dates_to_extract) 
    
    var_name_map = {"SPM": "suspended_matters", "CHL": "chlorophyll_a"}
    var_name = next((v for k, v in var_name_map.items() if k in self.info.Satellite_variable), None)
    
    if var_name is None : 
        self.download_report.loc[self.dates_to_extract.strftime("%Y-%m-%d")] = '❌ Unknown parameter name to extract'
        return
    
    formatted_dates = self.dates_to_extract.strftime("%Y%m%d")
    report_index = self.dates_to_extract.strftime("%Y-%m-%d")
    
    # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
    with multiprocessing.Pool(self.nb_of_cores_to_use) as pool:
        
        download_reports = pool.starmap(extract_values_from_one_global_file, 
                                        [(formatted_dates, report_index, date_to_extract,
                                          self.global_files_already_downloaded,
                                          self.download_report, var_name,
                                          where_to_save_values) 
                                                   for date_to_extract in formatted_dates ])
        
    # for date_to_extract in formatted_dates : 
    #     extract_values_from_one_global_file(formatted_dates, report_index, date_to_extract,
    #                                       self.global_files_already_downloaded,
    #                                       self.download_report, var_name,
    #                                       where_to_save_values) 
        
    self.download_report = pd.concat(download_reports, axis = 1).T.sort_index()
    
    
    
def extract_values_from_one_global_file(formatted_dates, report_index, date_to_extract,
                                        global_files_already_downloaded,
                                        download_report,var_name,
                                        where_to_save_values) :
        
    # matched_files = global_files_already_downloaded[np.char.endswith(global_files_already_downloaded, f'{date_to_extract}.nc')]
    matched_files = [f for f in global_files_already_downloaded if f.endswith(f'{date_to_extract}.nc')]
          
    if len(matched_files) == 0:
        download_report.loc[date_to_extract, 'Message'] = '❌ Satellite global file not found'
        return download_report.loc[date_to_extract]
    
    global_file = matched_files[0]
                    
    try : 
        with xr.open_dataset(global_file) as ds :
            time = ds.attrs['image_reference_time']
            map_ini = ds.squeeze('time')[var_name].drop_vars('time').to_dataset().assign_attrs({'time' : time})
                                                                
    except Exception : 
        download_report.loc[date_to_extract, 'Message'] = '❌ Error while extracting the values : {e}'
        return download_report.loc[date_to_extract]
                        
    try : 
        
        where_to_save_values_of_the_file = where_to_save_values[ np.where( date_to_extract == formatted_dates )[0][0] ]
        os.makedirs(where_to_save_values_of_the_file, exist_ok=True)
        map_ini.to_netcdf(f"{where_to_save_values_of_the_file}/{var_name}.nc",
                                encoding={map_ini[var_name].name: {"zlib": True, "complevel": 4, "dtype": "float32"}},
                                engine = 'netcdf4')
        download_report.loc[date_to_extract,'Message'] = f"✅ Downloaded and Extracted: {var_name}.nc"
        return download_report.loc[date_to_extract]
        
    except Exception : 
        download_report.loc[date_to_extract,'Message'] = '❌ Error while saving the file : {e}'
        return download_report.loc[date_to_extract]




# =============================================================================
#### Classes 
# =============================================================================

class download_satellite_data : 
        
    def __init__(self, info, start_day, end_day, where_to_save_satellite_data, nb_of_cores_to_use, overwrite_existing_satellite_files = False)  :
        
        info['Satellite_variable_name_on_remote_folder'] = format_variable_name_with_server_names( info )
                
        destination_path_to_fill = path_to_fill_to_where_to_save_satellite_files(where_to_save_satellite_data)
                    
        start_day_formatted = datetime.strptime(start_day, "%Y/%m/%d")
        end_day_formatted = datetime.strptime(end_day, "%Y/%m/%d")
        
        if np.isin(info.Temporal_resolution, ['DAILY', 'MONTHLY']) : 
        
            all_dates = pd.date_range((start_day_formatted if info.Temporal_resolution == 'DAILY' 
                                       else start_day_formatted.replace(day = 1) if info.Temporal_resolution == 'MONTHLY'
                                       else None), 
                                      (end_day_formatted if info.Temporal_resolution == 'DAILY' 
                                        else end_day_formatted.replace(day = 1) if info.Temporal_resolution == 'MONTHLY'
                                        else None),
                                      freq={'DAILY' : 'D',
                                            'MONTHLY' : 'MS'}[info.Temporal_resolution])
            
        elif info.Temporal_resolution == 'WEEKLY' : 
            
            all_dates = [pd.date_range(start=f"{year}-01-01", end=f"{year}-12-31", freq="8D") 
                         for year in np.arange(start_day_formatted.year, end_day_formatted.year +1)]

            all_dates = pd.to_datetime( np.concatenate(all_dates) )
            
            all_dates = all_dates[(all_dates >= (start_day_formatted - timedelta(days = 7))) & (all_dates <= end_day_formatted)]
            
        Download_messages = pd.DataFrame({'Message' : ['' for _ in all_dates]}, index = all_dates) 
                                        
        if test_if_satellite_data_DOES_NOT_exist_in_the_Data_source(info) :
            self.to_process = False
            Download_messages.loc[:,'Message'] = "The Data source does not distribute such data"
            self.download_report = Download_messages
            self.destination_path_to_fill = destination_path_to_fill
            return 
        
        info['Year'] = np.unique(all_dates.year)
                   
        if overwrite_existing_satellite_files == False : 
                        
            files_in_the_time_range = find_files_in_the_time_range(info, destination_path_to_fill, all_dates)

        else : 
            
            files_in_the_time_range = np.array([])
       
        dates_to_download = which_dates_occured_in_the_sat_file_names(all_dates, files_in_the_time_range, 
                                                                      check_for_non_occurrence_instead = True)
        
        dates_to_extract = dates_to_download.copy()
        
        if (info.Data_source == 'SEXTANT') and (info.sensor_name != 'merged') :
            
            global_files_in_the_time_range = find_files_in_the_time_range(info, destination_path_to_fill.replace('[PARAMETER]', 'ALL'), all_dates)
            
            dates_to_download = which_dates_occured_in_the_sat_file_names(all_dates, global_files_in_the_time_range, 
                                                                          check_for_non_occurrence_instead = True)
        
        Download_messages.loc[Download_messages.index.isin(dates_to_extract) == False, 'Message'] = "✅ Already downloaded"
                
        self.info = info
        self.start_day = start_day
        self.end_day = end_day
        self.where_to_save_satellite_data = where_to_save_satellite_data
        self.all_dates = all_dates
        self.files_already_downloaded = files_in_the_time_range
        self.dates_to_download = dates_to_download
        self.destination_path_to_fill = destination_path_to_fill
        self.download_report = Download_messages
        self.nb_of_cores_to_use = nb_of_cores_to_use
        
        if len(dates_to_extract) == 0 :
            self.to_process = False
        else : 
            self.to_process = True  
            
        if (info.Data_source == 'SEXTANT') and (info.sensor_name != 'merged') :
            self.global_files_already_downloaded = global_files_in_the_time_range
            self.dates_to_extract = dates_to_extract
            

        
    def download_missing_satellite_data(self)  :
        
        ### Create path for each file to download. 
        # Create the remote folder

        if self.info.Data_source == 'EUMETSAT' : 
            
            self.splitted_dates_to_download = split_consecutive_dates(self.dates_to_download)
            
            if "SST" in self.info.Satellite_variable : 
                new_map_resolution = 600
            else :
                new_map_resolution = 300
            
            download_EUMDAC_maps(self, 
                                 min_lon = -8, 
                                 max_lon = 11, 
                                 min_lat = 41, 
                                 max_lat = 52,
                                 new_map_resolution = new_map_resolution)
            
        else : 
            
            download_files_from_ftp(self)
            
            if (self.info.Data_source == 'SEXTANT') and (self.info.sensor_name != 'merged') :
                extract_values_from_the_global_file(self)
                                
        return None
    
