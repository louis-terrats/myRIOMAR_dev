import gc, os

import rpy2.robjects as robjects

from _1_data_validation.utils import MU_database_processing

from utils import (store_arguments, get_all_cases_to_process)

import pandas as pd
import numpy as np


def Match_up_with_insitu_measurements(core_arguments, redo_the_MU_database, nb_of_cores_to_use,
                                      where_are_saved_satellite_data, where_to_save_Match_Up_data) : 
             
    cases_to_process = get_all_cases_to_process(core_arguments)
               
    MU_database = MU_database_processing(where_to_save_Match_Up_data = where_to_save_Match_Up_data, 
                              cases_to_process = cases_to_process,
                              redo_the_MU_database = redo_the_MU_database,
                              nb_of_cores_to_use = nb_of_cores_to_use)
    
    if len(MU_database.MU_database) == 0 :
        MU_database.Create_the_MU_database(where_are_saved_satellite_data)
    else : 
        MU_database.Complete_the_MU_database(where_are_saved_satellite_data)
                
    MU_database.Summarize_MU_with_statistics()
    
    MU_database.Compile_and_Save_MU_summary( core_arguments['Data_sources'] )
    
    del MU_database
    gc.collect()


def Make_validation_plots_and_statistics(core_arguments, where_are_saved_Match_Up_data, MU_criteria, nb_of_cores_to_use):
    
    cases_to_process = get_all_cases_to_process(core_arguments)
    
    path_to_summary_df = where_are_saved_Match_Up_data + '/MU_summary.csv'
    
    # Load the summary_df
    if os.path.exists(path_to_summary_df) : 
        MU_summary_df = pd.read_csv(path_to_summary_df)
    else : 
        return f"No Summary file here : {self.path_to_the_MU_database.replace('database.joblib', '') + 'summary.csv'}"
    
    # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
    with multiprocessing.Pool(nb_of_cores_to_use) as pool:

        def make_scatterplot_and_save_statistics(info) : 
    
            index_to_keep = np.where((MU_summary_df.Data_source == info.Data_source) & 
                                     (MU_summary_df.sensor_name == info.sensor_name) & 
                                     (MU_summary_df.atmospheric_correction == info.atmospheric_correction) & 
                                     (MU_summary_df.Satellite_variable == info.Satellite_variable))
            
            MU_summary_df_of_the_case = MU_summary_df.iloc[index_to_keep]
            
            MU_criteria = get_MU_criteria_for_each_product(info)
            
            # Source the R script
            robjects.r['source']("myRIOMAR_dev/_1_data_validation/utils.R")
        
            r_function = robjects.r['Main_function']
            
            for satellite_algorithm in np.unique(MU_summary_df_of_the_case.Satellite_algorithm) : 
                
                MU_summary_df_of_the_sat_algo = MU_summary_df_of_the_case[ MU_summary_df_of_the_case.Satellite_algorithm == satellite_algorithm ]
            
                MU_summary_df_of_the_sat_algo = MU_summary_df_of_the_sat_algo.iloc[:10]
                    
                # Call the R function
                r_function(
                    
                    satellite_median = robjects.FloatVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_mean'].to_list()),
                    satellite_n = robjects.IntVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_n'].to_list()),
                    satellite_sd = robjects.FloatVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_std'].to_list()),
                    satellite_times = robjects.StrVector(MU_summary_df_of_the_sat_algo.Satellite_time.astype(str).to_list()),
                    insitu_value = robjects.FloatVector(MU_summary_df_of_the_sat_algo.Insitu_value.to_list()),
                    insitu_time = robjects.StrVector(MU_summary_df_of_the_sat_algo.Insitu_time.to_list()),
                    site_name = robjects.StrVector(MU_summary_df_of_the_sat_algo.SITE.to_list()),
                    min_n = robjects.IntVector([MU_criteria["min_n"]]),
                    max_CV = robjects.IntVector([MU_criteria["max_CV"]]),
                    max_hour_diff = robjects.IntVector([MU_criteria["max_hour_diff_between_insitu_and_satellite_measurement"]]),
                    grid_size = robjects.IntVector([MU_criteria["grid_size"]]),
                    date = robjects.StrVector(MU_summary_df_of_the_sat_algo.DATE.to_list()),
                    satellite_source = robjects.StrVector([info.Data_source]),
                    satellite_sensor = robjects.StrVector([info.sensor_name]),
                    satellite_atm_corr = robjects.StrVector([info.atmospheric_correction]),
                    satellite_algorithm = robjects.StrVector([satellite_algorithm]),
                    where_to_save_MU_results = robjects.StrVector([where_are_saved_Match_Up_data])
                    
                )
                
                    
            

    
    
