import gc

from myRIOMAR._1_data_validation.utils import MU_database_processing

from myRIOMAR.utils import (store_arguments, get_all_cases_to_process)


def Match_up_with_insitu_measurements(core_arguments, redo_the_MU_database, nb_of_cores_to_use,
                                      where_to_save_Match_Up_data) : 
             
    cases_to_process = get_all_cases_to_process(core_arguments)
               
    MU_database = MU_database_processing(where_to_save_Match_Up_data = where_to_save_Match_Up_data, 
                              cases_to_process = cases_to_process,
                              redo_the_MU_database = redo_the_MU_database,
                              nb_of_cores_to_use = nb_of_cores_to_use)
    
    if len(MU_database.MU_database) == 0 :
        MU_database.Create_the_MU_database(path_to_SOMLIT_insitu_data, path_to_REPHY_insitu_data, where_to_save_satellite_data)
    else : 
        MU_database.Complete_the_MU_database(where_to_save_satellite_data)
                
    MU_database.Summarize_MU_with_statistics()
    
    MU_database.Compile_and_Save_MU_summary(Data_sources)
    
    del MU_database
    gc.collect()