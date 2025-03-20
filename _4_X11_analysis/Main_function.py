#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:40:41 2025

@author: terrats
"""

from myRIOMAR._4_X11_analysis.utils_functions import (F_EVF_V1_2, F_census_X_11_pezzulli_V1_2, f_test_season_v1_1, f_kendall_sen_v1_0, 
                                                   f_seasonal_kendall_v2_0, linfit, F_test_pente)
import pandas as pd
import numpy as np
import sys
from scipy.stats import f, kruskal, norm, chi2, t
from scipy import stats

def temporal_decomp_V2_7_x11(values, dates, time_frequency,
                             filter_outlier=False, overall_cutoff=50.0, 
                             out_limit=3.0, perc_month_limit=50.0, 
                             var_stationary=False, lin_interpol=False, 
                             cutoff_fill=30.0, season_test=False):
    
    """
    PURPOSE:
    Decompose a time series var into 3 components: T: TREND, S: SEASON, Y: Irregular
    
    ; PURPOSE:
    ; Decompose a time series var into 3 components: T: TREND, S: SEASON, Y: Irregular
    ;
    ; CATEGORY:
    ;
    ; CALLING SEQUENCE:
    ; temporal_decomp_V2_6_x11, Var, month, hbad
    ;
    ; INPUTS:
    ; Var:  vector representing the ANNUAL TIME SERIES  to be decomposed (12 months cycle including bad values)
    ; month: vector of month
    ; hbad: bad value
    ; data: ouput structure
    
    ; KEYWORD PARAMETERS:
    ;
    ;filter_outlier : if set filter of the outlier assuming a normal distribution (mean +/- (out_lim*std))
    ;out_limit : limit considered for outlier detection ...default: 3
    ;overall_cutoff: cutoff value for the test on the intial % of valid data in the TS. If %age of valid data < overall cutoff -> flagged. Default 50%
    ;perc_month_limit : cutoff value for selecting the valid months for the "shortened years". in % defualt 50%
    ;X11_pezzulli: if set use the X11 method defined in Pezzulli et al., J of Climate, 2005
    ;X11_user : if set use the X11 method defined in the X11 user guide manual
    ;var_stationary : if stat compute the contribution of the X11 components to the variance of the stationary part of the initial TS
    ;lin_interpol : if set to 1 fill the gaps in the TS using linear interpolation instead of the evf method (if set to 0)
    ;cutoff_fill : cutoff value for the maximal %age of missing data acceptable for performing the gap filling procedure. Default 30% of missing data max
    ;season_test : if set then the test on the presence of sesoanlity in the data is performed
    ;
    
    ;========================================================================================
    ;OUTPUT STRUCTURE DEFINITION
    
    ;PtStr = {	var_S:fbad , $ ; Percentage of variance of X due to the Seasonal component
    ;            var_T:fbad , $; Percentage of variance of X due to the Trend component
    ;            var_Y:fbad ,  $; Percentage of variance of X due to the Irregular component
    ;			covar: fbad, $; Percentage of variance of X related to the covariance terms
    ;
    ;			var_X: fbad, $;  variance of the original time series
    ;			var_coeff: fbad, $;  variation coefficient of the original time series
    ;
    ;			code_pres_month: 0, $ ; Integer coding the m_test_year vector of monthly data validity
    ;
    ;			S : make_array(N_months, value =fbad), $ ; vector of N months containing the Seasonal component
    ;			T : make_array(N_months, value =fbad), $ ; vector of N months containing the Trend component
    ;			Y : make_array(N_months, value =fbad), $  ; vector of N months containing the Irregular component
    ;
    ;			year_n_months :fbad ,$ ; Number of months in 1 synthetic year
    ;			N_outlier: 0. ,$ ; Number of outlier
    ;			perc_valid_months: fbad,  $ ; Percentage of valid initial data (Number of valid data in the initial time series/ 87) *100.
    ;			nb_valid_month:  fbad,  $ ;Total number of valid data considered in the analysis (after short_year and gap_filling)
    ;			N_missing_data: 0. ,$  ; Number of data which have been filled (including missing values, outlier...)
    ;			test_season: make_array(2, value = fbad) ,$ ; test the homogeneity of the season
    ;
    ;			Kendall_sen: make_array(5, value = fbad), $
    ;			Seasonal_Kendall_sen: make_array(5, value = fbad), $
    ;
    ;			slope_trend: fbad, $
    ;            intercept_trend : fbad, $               .
    ;            prob_Test_trend:fbad , $
    ;
    ;			flag: fbad  $; TEST the efficiency of the gap filling method
    ;		}
    
    ; SIDE EFFECTS:
    ; none
    ;
    ; COMMON BLOCKS:
    ; none
    ;
    ; MODIFICATION HISTORY:
    ;
    ; Vantrepotte V.
    ;Written  20.10.07
    ;Last Update: 26-09-2008: RC rate of Change + comment
    
    """
    
    # missing_values_are = np.nan # -99999
    # fbad = np.nan
        
    var = pd.Series(values, index= dates )

    if time_frequency == "ANNUAL" : 
        full_date_range = [pd.Timestamp(f'{date.year}-07-01') for date in 
                           pd.date_range(start=dates.min(), end=dates.max() + pd.offsets.MonthEnd(0), freq = '1Y')]
    
    if time_frequency == "MONTHLY" : 
        full_date_range = [pd.Timestamp(f'{date.year}-{date.month:02d}-15') for date in 
                           pd.date_range(start=dates.min(), end=dates.max() + pd.offsets.MonthEnd(0), freq = '1M')]
        
    if time_frequency == "WEEKLY" : 
        full_date_range = [pd.Timestamp(f'{date.year}-{date.month:02d}-{the_day:02d}') 
                             for date in pd.to_datetime( np.unique( pd.date_range(start=dates.min(), end=dates.max(), freq = '1D').strftime('%Y-%m') ) )
                             for the_day in [ 4, 12, 20, 28]]
        
    if time_frequency == "INITIAL" : 
        full_date_range = dates
        
    var_reindexed = var.reindex(full_date_range, fill_value=np.nan)  
    # var_reindexed[np.isnan(var_reindexed)] = missing_values_are
    
    N_values = len(var_reindexed)
    data = {
            "0_variance_due_to_Seasonal_signal": np.nan,  # Percentage of variance of X due to the Seasonal component
            "1_variance_due_to_Interannual_signal": np.nan,  # Percentage of variance of X due to the Trend component
            "2_variance_due_to_Residual_signal": np.nan,  # Percentage of variance of X due to the Irregular component
            "3_covar": np.nan,  # Percentage of variance of X related to the covariance terms
    
            "4_var_X": np.nan,  # Variance of the original time series
            "5_var_coeff": np.nan,  # Variation coefficient of the original time series
    
            "6_code_pres_month": 0,  # Integer coding the m_test_year vector of monthly data validity
    
            "7_dates": np.full(N_values, np.nan).astype(float),
            "8_values_ini": np.full(N_values, np.nan).astype(float),
            "8_values_ini_interpolated": np.full(N_values, np.nan).astype(float),
            
            "9_Seasonal_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Seasonal component
            "10_Interannual_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Trend component
            "11_Residual_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Irregular component
    
            "12_year_n_months": np.nan,  # Number of months in 1 synthetic year
            "13_N_outlier": 0.0,  # Number of outliers
            "14_perc_valid_data": np.nan,  # Percentage of valid initial data (Number of valid data in the initial time series / 87) * 100
            "15_nb_valid_month": np.nan,  # Total number of valid data considered in the analysis (after short_year and gap_filling)
            "16_N_missing_data": 0.0,  # Number of data which have been filled (including missing values, outliers, etc.)
            "17_Homogeneity_of_the_Seasonal_signal": np.full(2, np.nan),  # Test the homogeneity of the season
    
            # "Kendall_sen": np.full(5, np.nan),  # Kendall's tau statistics
            "18_Kendall_Sen_analyses_on_Interannual_signal": np.nan,  # Kendall's tau statistics
            "19_Kendall_Sen_analyses_on_Seasonal_signal": np.nan,  # Seasonal Kendall's tau statistics
    
            "20_slope_trend_per_time_step": np.nan,  # Slope of the linear trend derived from Tt
            
            "22_intercept_trend": np.nan,  # Intercept of the linear trend derived from Tt
            "23_prob_Test_trend": np.nan,  # Probability associated with the linear trend derived from Tt
    
            "24_flag": np.nan,  # Test the efficiency of the gap filling method
            
            "25_Is_the_slope_trend_per_time_step_significant": np.nan,
        }
    
    data['7_dates'] = var_reindexed.index
    data['8_values_ini'] = var_reindexed.values
    
    # Calculate the percentage of valid months
    perc_valid_data = 100 - np.sum(np.isnan(var_reindexed)) / len(var_reindexed) * 100.0
    data['14_perc_valid_data'] = perc_valid_data

    if perc_valid_data < overall_cutoff :
        
        data['24_flag'] = -100
        return data

    if time_frequency in ["MONTHLY", "ANNUAL"] : 
        months = np.array([x.month for x in var_reindexed.index])
    if time_frequency == "WEEKLY" :     
        months = np.array([x.strftime('%m-%d') for x in var_reindexed.index])
        
    values_to_use = np.zeros(len(var_reindexed)).astype(bool)
    
    n_months_per_year = 12
    if time_frequency == "WEEKLY" : 
        n_months_per_year = int(n_months_per_year * 4)
    if time_frequency == "ANNUAL" : 
        n_months_per_year = int(n_months_per_year / 12)
        
    index_month_to_use = np.zeros(n_months_per_year).astype(bool)
    
    for i_m in np.unique(months) :
        
        ind_month_i = np.where(months == i_m)[0]
        ind_bad_val_month_i = np.where( np.isnan(var_reindexed[ind_month_i]) )[0]
        n_bad_val_month_i = len(ind_bad_val_month_i)
        count_n_month_i = len(ind_month_i)

        if (100 * float(n_bad_val_month_i) / float(count_n_month_i)) < perc_month_limit :
            index_month_to_use[ind_month_i[0]] = True
            values_to_use[ind_month_i] = True
    
    if all(index_month_to_use == False) : 
        data['24_flag'] = -101
        return data
    
    data['12_year_n_months'] = np.sum(index_month_to_use)
    data['15_nb_valid_month'] = np.sum(values_to_use)
    if time_frequency == "MONTHLY" : 
        data['14_perc_valid_month'] = 100 * data['15_nb_valid_month'] / len(np.unique([x.strftime('%Y-%m') for x in var_reindexed.index]))
    if time_frequency == "WEEKLY" :     
        data['14_perc_valid_month'] = 100 * data['15_nb_valid_month'] / len(np.unique([x.strftime('%Y-%m-%d') for x in var_reindexed.index]))

    index_var_to_use = np.where(values_to_use)[0]
    var_to_use = var_reindexed.iloc[ index_var_to_use ]
    
    if time_frequency in ["MONTHLY", "ANNUAL"] : 
        months = np.array([x.month for x in var_reindexed.index])
    if time_frequency == "WEEKLY" :     
        months = np.array([x.strftime('%m-%d') for x in var_reindexed.index])
    months = np.unique(months)

    n_outlier = 0
    
    if filter_outlier : 
    
        for month in months :
            
            var_month = var_to_use.iloc[ np.where([x == month for x in months])[0] ]
            
            # var_month_to_use_for_mean_std = var_month[np.isfinite(var_month)]
            mean_month = np.nanmean(var_month)
            std_month = np.nanstd(var_month)
            
            index_outlier = np.where((var_month > mean_month + out_limit * std_month) | 
                                     (var_month < mean_month - out_limit * std_month))[0]
    
            n_outlier += len(index_outlier)
    
            if len(index_outlier) != 0:                
                dates_outliers = var_month.iloc[index_outlier].index
                var_to_use.loc[dates_outliers] = np.nan
    
    # GAPS FILLING using EVF METHOD
    test_missing = np.where( np.isnan(var_to_use) )[0]
    count_missing = len(test_missing)

    var_interpolated = var_to_use
    if count_missing != 0:
        
        if (100 * float(count_missing) / len(var_to_use)) < cutoff_fill :
            
            if lin_interpol == False :
                
                # Use EVF method
                max_iter = 20
                level_SD_perc = 10
                ACP_cutoff = 80

                var_interpolated = F_EVF_V1_2(X = var_to_use, 
                                              level_SD_perc = level_SD_perc, ACP_cutoff = ACP_cutoff, max_iter = max_iter)
                data['24_flag'] = 1

                if (count_missing > 1) and (var_interpolated[test_missing[0]] == var_interpolated[test_missing[1]]) :
                                                
                    data['24_flag'] = 2
                    level_SD_perc = 10
                    ACP_cutoff = 50
                    var_interpolated[test_missing] = np.nan
                    var_interpolated = F_EVF_V1_2(X = var_interpolated,
                                                  level_SD_perc = level_SD_perc, ACP_cutoff = ACP_cutoff, max_iter = max_iter)

                    if var_interpolated[test_missing[0]] - var_interpolated[test_missing[1]] == 0:
                        data['24_flag'] = 3

            else:
                # Use linear interpolation
                var_interpolated = np.interp(np.arange(len(var_to_use)), 
                                    np.arange(len(var_to_use))[ np.isfinite(var_to_use) ], 
                                    var_to_use[ np.isfinite(var_to_use) ])
                data['24_flag'] = 4
        else:
            data['24_flag'] = -10
    else:
        data['24_flag'] = 0

    data['8_values_ini_interpolated'] = var_interpolated

    test_NAN = np.isnan(var_interpolated).sum()
    if test_NAN != 0:
        data['24_flag'] = -9

    # TIME SERIES analyses
    if test_NAN != 0 or data['24_flag'] == -10 or (time_frequency != "ANNUAL" and data['12_year_n_months'] <= 3) :
        
        return data
                
    Census_X11 = F_census_X_11_pezzulli_V1_2(Xt = var_interpolated, period = data['12_year_n_months'])

    data['9_Seasonal_signal'][index_var_to_use] = Census_X11['St']
    data['11_Residual_signal'][index_var_to_use] = Census_X11['Yt']
    data['10_Interannual_signal'][index_var_to_use] = Census_X11['Tt']

    if season_test:
        data['17_Homogeneity_of_the_Seasonal_signal'] = f_test_season_v1_1(
            X = var_interpolated[:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            T = Census_X11['Tt'][:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            S = Census_X11['St'][:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            period = data['12_year_n_months'])

    data['4_var_X'] = np.var(var_interpolated)
    data['5_var_coeff'] = np.std(var_interpolated) / np.mean(var_interpolated) * 100.0

    if var_stationary:
        
        fit_trend_LINEAR = np.polyfit(np.arange(len(Census_X11['Tt'])), Census_X11['Tt'], 1)
        yfit = np.polyval(fit_trend_LINEAR, np.arange(len(Census_X11['Tt'])))
        data['0_variance_due_to_Seasonal_signal'] = np.var(Census_X11['St']) / np.var(var_interpolated - yfit) * 100.0
        data['1_variance_due_to_Interannual_signal'] = np.var(Census_X11['Tt'] - yfit) / np.var(var_interpolated - yfit) * 100.0
        data['2_variance_due_to_Residual_signal'] = np.var(Census_X11['Yt']) / np.var(var_interpolated - yfit) * 100.0
        data['3_covar'] = 100.0 * (1.0 - np.sum([np.var(Census_X11['St']) / np.var(var_interpolated - yfit), 
                                                np.var(Census_X11['Tt'] - yfit) / np.var(var_interpolated - yfit), 
                                                np.var(Census_X11['Yt']) / np.var(var_interpolated - yfit)]))
        
    else:
        
        var_St = np.var(Census_X11['St'])
        var_Tt = np.var(Census_X11['Tt'])
        var_Yt =np.var(Census_X11['Yt'])
        sum_vars = var_St + var_Tt + var_Yt
        
        data['0_variance_due_to_Seasonal_signal'] = var_St / sum_vars * 100.0
        data['1_variance_due_to_Interannual_signal'] = var_Tt / sum_vars * 100.0
        data['2_variance_due_to_Residual_signal'] = var_Yt / sum_vars * 100.0
        data['3_covar'] = 100.0 * (1.0 - np.sum([np.var(Census_X11['St']) / np.var(var_interpolated), 
                                                np.var(Census_X11['Tt']) / np.var(var_interpolated), 
                                                np.var(Census_X11['Yt']) / np.var(var_interpolated)]))
       
    if time_frequency != 'ANNUAL' :
        v_mean_year = [np.mean(var_interpolated[ int(i * data['12_year_n_months']) : int( (i + 1) * data['12_year_n_months'] -1) ]) 
                       for i in np.arange(len(var_interpolated) / data['12_year_n_months'])
                       if len(var_interpolated[ int(i * data['12_year_n_months']) : int( (i + 1) * data['12_year_n_months'] ) ]) == data['12_year_n_months']]
    else : 
        v_mean_year = var_interpolated
    
    if len( np.where(np.isfinite(v_mean_year))[0] ) > 2 :
        data['18_Kendall_Sen_analyses_on_Interannual_signal'] = f_kendall_sen_v1_0(X = v_mean_year, alpha=0.05, correction_tie = True)
        
    data['19_Kendall_Sen_analyses_on_Seasonal_signal'] = f_seasonal_kendall_v2_0(var_interpolated, data['12_year_n_months'], alpha=0.05)
    
    # Generate x-values corresponding to the length of census_X11.Tt
    x_values = np.arange(len(Census_X11['Tt']))
    
    # Perform a linear fit
    fit_trend_test = linfit(x_values, Census_X11['Tt'])
    
    # Test the significance of the slope
    test_slope = F_test_pente(x_values, Census_X11['Tt'], fit_trend_test[1])
    
    # Store the slope and intercept in a data object (assuming 'data' is a pre-defined object with the necessary attributes)
    data['20_slope_trend_per_time_step'] = fit_trend_test[1]
    data['22_intercept_trend'] = fit_trend_test[0]
    
    # Calculate the significance level
    data['15_nb_valid_month'] = len(Census_X11['Tt'])  # Assuming this represents the number of valid months
    data['21_Significance_level_of_the_slope'] = 1.0 - abs(1.0 - 2.0 * (1.0 - t.cdf(abs(test_slope[1]), data['15_nb_valid_month'] - 2)))
        
    data['25_Is_the_slope_trend_per_time_step_significant'] = data['21_Significance_level_of_the_slope'] < 0.05
    
    return data

