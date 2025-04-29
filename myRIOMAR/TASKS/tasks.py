#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 11:30:15 2025

@author: loulou
"""

import os
import glob
import yaml
import sys

# --- Define PACKAGE_ROOT and add to sys.path ---
PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PACKAGE_ROOT)

# --- Import your package ---
import myRIOMAR_dev as myRIOMAR

# --- Load YAML configuration ---
yaml_files = glob.glob(os.path.join(PACKAGE_ROOT, "CONFIG", "*.yaml"))

core_arguments = {}
for yaml_file in yaml_files:
    with open(yaml_file, "r") as f:
        content = yaml.safe_load(f)
        if content:
            # Merge dictionaries
            core_arguments.update(content)

# --- Define common variables ---
RESULTS_DIR = core_arguments['Paths']['Where_to_save_results']
GLOBAL_MAPS_dir = core_arguments['Paths']['Where_to_save_satellite_global_map']
nb_cores = core_arguments['Multithreading']['nb_of_cores_to_use']
use_dynamic_threshold_for_plume_detection = core_arguments['Plume_detection']['Use_dynamic_threshold']


# --- TASK FUNCTIONS ---

def download_data(core_arguments, 
                  overwrite=False):
    myRIOMAR._0_data_downloading.Download_satellite_data(
        core_arguments,
        nb_of_cores_to_use=nb_cores,
        overwrite_existing_satellite_files=overwrite,
        where_to_save_satellite_data=GLOBAL_MAPS_dir,
    )

def plot_global_maps(core_arguments, 
                     start_day=None, end_day=None):
    myRIOMAR._0_data_downloading.Plot_and_Save_the_map(
        core_arguments,
        nb_of_cores_to_use=nb_cores,
        where_are_saved_satellite_data=GLOBAL_MAPS_dir,
        start_day_of_maps_to_plot=start_day,
        end_day_of_maps_to_plot=end_day,
    )

def match_up_with_insitu(core_arguments,
                         redo_the_MU_database=False):
    myRIOMAR._1_data_validation.Match_up_with_insitu_measurements(
        core_arguments,
        zones=["FRANCE"],
        redo_the_MU_database=redo_the_MU_database,
        nb_of_cores_to_use=nb_cores,
        where_are_saved_satellite_data=RESULTS_DIR,
        where_to_save_Match_Up_data=os.path.join(RESULTS_DIR, "MU_data"),
    )

def create_regional_maps(core_arguments, 
                         overwrite = True,
                         time_frequencies_to_map = {'DAILY': False, 'WEEKLY': False, 'MONTHLY': True, 'ANNUAL': True}):
    
    myRIOMAR._2_regional_maps.create_regional_maps(
        core_arguments,
        Zones=["SOUTHERN_BRITTANY"],
        overwrite_existing_regional_maps=overwrite,
        save_map_plots_of_which_time_frequency=time_frequencies_to_map,
        nb_of_cores_to_use=nb_cores,
        where_are_saved_satellite_data=GLOBAL_MAPS_dir,
        where_to_save_regional_maps=RESULTS_DIR,
    )
    
def detect_plumes(core_arguments):
    
    myRIOMAR._3_plume_detection.apply_plume_mask(
        core_arguments,
        Zones = ['SOUTHERN_BRITTANY'],
        detect_plumes_on_which_temporal_resolution_data = 'WEEKLY',
        nb_of_cores_to_use = nb_cores,
        use_dynamic_threshold = use_dynamic_threshold_for_plume_detection,
        where_are_saved_regional_maps = RESULTS_DIR,
        where_to_save_plume_results = os.path.join(RESULTS_DIR, ('DYNAMIC' if use_dynamic_threshold_for_plume_detection else 'FIXED') + '_THRESHOLD'),
    )


def plot_timeseries_of_plume_surface(core_arguments,
                                     temporal_resolution_to_use = 'WEEKLY'):
    
    myRIOMAR._3_plume_detection.make_and_plot_time_series_of_plume_areas(
        core_arguments,
        Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
        nb_of_cores_to_use = nb_cores,
        on_which_temporal_resolution_the_plumes_have_been_detected = temporal_resolution_to_use,
        where_are_saved_plume_results = os.path.join(RESULTS_DIR, ('DYNAMIC' if use_dynamic_threshold_for_plume_detection else 'FIXED') + '_THRESHOLD'),
        where_to_save_plume_time_series = os.path.join(RESULTS_DIR, ('DYNAMIC' if use_dynamic_threshold_for_plume_detection else 'FIXED') + '_THRESHOLD'),
    )

def X11_analysis_on_time_series(core_args):

    myRIOMAR._4_X11_analysis.Apply_X11_method_on_time_series(
        core_arguments,
        Zones = core_arguments['Product_selection']['Regions'],
        nb_of_cores_to_use = nb_cores,
        on_which_temporal_resolution_the_plumes_have_been_detected = core_arguments['Temporal_resolution']['Temporal_resolution_of_plume_detection'],
        where_are_saved_plume_time_series = os.path.join(RESULTS_DIR, ('DYNAMIC' if use_dynamic_threshold_for_plume_detection else 'FIXED') + '_THRESHOLD'),
        where_to_save_X11_results = os.path.join(RESULTS_DIR, ('DYNAMIC' if use_dynamic_threshold_for_plume_detection else 'FIXED') + '_THRESHOLD'),
        include_river_flow = True,
    )
