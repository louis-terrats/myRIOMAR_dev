#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:40:41 2025

@author: terrats
"""

from _4_X11_analysis.utils_functions import (temporal_decomp_V2_7_x11)

### S'INSPIRER DE LA FONCTION DNAS LE DOSSIER SCRIPTS

def Apply_X11_method_on_time_series(core_arguments, apply_X11_to_which_time_series)

results = temporal_decomp_V2_7_x11(values = ts_data[var_to_use].tolist(), dates = ts_data.index, 
                                   time_frequency = Time_resolution,
                                    filter_outlier=False, overall_cutoff=50, 
                                    out_limit=3, perc_month_limit=50, 
                                    var_stationary=False, lin_interpol=False, 
                                    cutoff_fill=30, season_test=True)


fig, axs = plt.subplots(4, 1, figsize=(24, 14))

axs[0].plot(results['7_dates'], results['8_values_ini'])
axs[0].set_title('Initial values')
axs[1].plot(results['7_dates'], results['10_Interannual_signal'])
axs[1].set_title(f'Inter-annual signal ({round(results["1_variance_due_to_Interannual_signal"], 1)}%)')
axs[2].plot(results['7_dates'], results['9_Seasonal_signal'])
axs[2].set_title(f'Seasonal signal ({round(results["0_variance_due_to_Seasonal_signal"], 1)}%)')
axs[3].plot(results['7_dates'], results['11_Residual_signal'])
axs[3].set_title(f'Residual signal ({round(results["2_variance_due_to_Residual_signal"], 1)}%)')

# Adjust layout
plt.tight_layout()

plt.savefig(where_to_save_results + apply_X11_to_which_time_series + "_" + var_to_use.replace(' ', '_') + '.png')

plt.close(fig)

pd.DataFrame({'dates' : results['7_dates'],
                'Raw_signal' : results['8_values_ini'],
                'Interannual_signal' : results['10_Interannual_signal'],
                'Seasonal_signal' : results['9_Seasonal_signal'],
                'Residual_signal' : results['11_Residual_signal'],
                'Variation_coefficient' : results['5_var_coeff'],
                'Variance_due_to_Interannual_signal' : results['1_variance_due_to_Interannual_signal'],
                'Variance_due_to_Seasonal_signal' : results['0_variance_due_to_Seasonal_signal'],
                'Variance_due_to_Residual_signal' : results['2_variance_due_to_Residual_signal'],
                'Monotonic_change' : results['18_Kendall_Sen_analyses_on_Interannual_signal']['Is_the_change_of_annual_values_monotonic'],
                'Rate_of_Change' : results['18_Kendall_Sen_analyses_on_Interannual_signal']['Rate_of_change_of_annual_values_in_percentage_per_time']
                }).to_csv(where_to_save_results + apply_X11_to_which_time_series + "_" + var_to_use.replace(' ', '_') + '.csv', index = False)

