from utils import (get_all_cases_to_process)
from _1_data_validation.utils import MU_database_processing


def Match_up_with_insitu_measurements(core_arguments, zones, redo_the_MU_database, nb_of_cores_to_use,
                                      where_are_saved_satellite_data, where_to_save_Match_Up_data) : 
             
    core_arguments.update({'Temporal_resolution' : ['DAILY']})
    
    cases_to_process = get_all_cases_to_process(core_arguments)
               
    MU_database = MU_database_processing(where_to_save_Match_Up_data = where_to_save_Match_Up_data, 
                              cases_to_process = cases_to_process,
                              zones = zones,
                              redo_the_MU_database = redo_the_MU_database,
                              nb_of_cores_to_use = nb_of_cores_to_use)
    
    if len(MU_database.MU_database) == 0 :
        MU_database.Create_the_MU_database(where_are_saved_satellite_data)
    else : 
        MU_database.Complete_the_MU_database(where_are_saved_satellite_data)
      
    MU_database.Summarize_MU_with_statistics()
    
    MU_database.Compile_and_Save_MU_summary( core_arguments['Data_sources'] )
            
    MU_database.Make_scatterplot_and_save_statistics()
    
    # MU_database.Make_time_series()
