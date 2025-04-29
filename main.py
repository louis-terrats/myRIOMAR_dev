#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:41:48 2025

@author: terrats
"""

import argparse
from myRIOMAR_dev.CODES.TASKS import tasks
import matplotlib

# Set matplotlib backend to prevent plots from displaying
# matplotlib.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM memory !)
matplotlib.use('agg') # Prevent showing plot in the Plot panel (this saves RAM memory)

def main():
    
    parser = argparse.ArgumentParser(
        description="myRIOMAR: Satellite river plume data analysis tools"
    )

    subparsers = parser.add_subparsers(dest='command')

    # --- download ---
    parser_download = subparsers.add_parser('download', help='Download satellite data')
    parser_download.add_argument('--overwrite', action='store_true', help='Overwrite existing data')

    # --- plot ---
    parser_plot = subparsers.add_parser('plot', help='Plot global maps')
    parser_plot.add_argument('--start_day', type=str, default=None, help='Start date (YYYY-MM-DD)')
    parser_plot.add_argument('--end_day', type=str, default=None, help='End date (YYYY-MM-DD)')

    # --- match-up ---
    parser_matchup = subparsers.add_parser('matchup', help='Match satellite with in-situ data')
    parser_matchup.add_argument('--redo', action='store_true', help='Redo the Match-Up database')

    # --- regional maps ---
    parser_maps = subparsers.add_parser('maps', help='Create regional maps')
    parser_maps.add_argument('--overwrite', action='store_true', help='Overwrite regional maps')

    # --- plume detection ---
    parser_plume = subparsers.add_parser('plumes', help='Detect river plumes')
    parser_plume.add_argument('--dynamic', action='store_true', help='Use dynamic threshold')

    # --- timeseries ---
    parser_timeseries = subparsers.add_parser('timeseries', help='Plot time series of plume surface')

    # --- X11 ---
    parser_x11 = subparsers.add_parser('X11', help='Apply X11 decomposition on time series')

    # === Parse arguments ===
    args = parser.parse_args()

    # === Call the corresponding tasks ===
    if args.command == 'download':
        tasks.download_data(tasks.core_arguments, overwrite=args.overwrite)
    elif args.command == 'plot':
        tasks.plot_global_maps(tasks.core_arguments, start_day=args.start_day, end_day=args.end_day)
    elif args.command == 'matchup':
        tasks.match_up_with_insitu(tasks.core_arguments, redo_the_MU_database=args.redo)
    elif args.command == 'maps':
        tasks.create_regional_maps(tasks.core_arguments, overwrite=args.overwrite)
    elif args.command == 'plumes':
        tasks.detect_plumes(tasks.core_arguments, use_dynamic_threshold=args.dynamic)
    elif args.command == 'timeseries':
        tasks.plot_timeseries_of_plume_surface(tasks.core_arguments)
    elif args.command == 'X11':
        tasks.X11_analysis_on_time_series(tasks.core_arguments)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()



# #%% Figures_for_the_article

# myRIOMAR._5_Figures_for_article.Figure_1(where_are_saved_satellite_data = "/media/terrats/My Book/LOUIS_TERRATS/RIOMAR/DATA/OCEAN_COLOR/",
#                                          where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

# myRIOMAR._5_Figures_for_article.Figure_2(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
#                                          where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

# myRIOMAR._5_Figures_for_article.Figure_4(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
#                                          where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

# myRIOMAR._5_Figures_for_article.Figure_5(where_are_saved_regional_maps = "/home/terrats/Desktop/RIOMAR/TEST/",
#                                          where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

# myRIOMAR._5_Figures_for_article.Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD/",
#                                            where_are_saved_plume_results_with_fixed_threshold = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/FIXED_THRESHOLD/",
#                                            where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")

# myRIOMAR._5_Figures_for_article.Figure_8_9_10(
#                                            where_are_saved_X11_results = "/home/terrats/Desktop/RIOMAR/TEST/RESULTS/DYNAMIC_THRESHOLD/",
#                                            where_to_save_the_figure = "/home/terrats/Desktop/RIOMAR/TEST/")




