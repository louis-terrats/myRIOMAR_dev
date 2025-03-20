import numpy as np
import pandas as pd
import sys
from scipy.stats import f, kruskal, norm, chi2, t
from scipy import stats



def generate_lagged_matrix(series, p):
    n = len(series)
    lagged_matrix = np.zeros((n - p, p))
    for i in range(p):
        lagged_matrix[:, i] = series[p - i - 1:n - i - 1]
    return lagged_matrix


def autocorrelation(Z, lags, Wt0 = None):
   
    Autocov = np.zeros_like(lags, dtype=float)
    # rk = np.zeros_like(lags, dtype=float) 
    if isinstance(Wt0, np.ndarray) == False : 
        Wt0 = np.ones(len(Z))
    
    # Autocovariance calculation for the FIRST ITERATION
    for i, lag in enumerate(lags[:-1]) :
        
        if i == 0 : 
            Zt = Z
            Zt_K = Z
            Wt_lag = Wt0
            Wt_lag_k = Wt0
        else : 
            Zt = Z[:-lag]
            Zt_K = Z[lag:]
            Wt_lag = Wt0[:-lag]
            Wt_lag_k = Wt0[lag:]

        Wt_global = Wt_lag + Wt_lag_k -1
        # Z_mean_k = np.nansum(Zt * Wt_global) / np.nansum(Wt_lag)
        Z_mean = np.nansum(Zt_K * Wt_global) / np.nansum(Wt_lag_k)

        Autocov[i] = (1.0 / np.nansum(Wt_lag)) * np.nansum(Wt_lag * (Zt - Z_mean) * (Zt_K - Z_mean))
        # rk[i] = 1.96 / np.sqrt(N - lag - 1)
        # Autocov[i] = 1.0 / (N - 1) * np.sum(Wt_lag * (Zt - Z_mean) * (Zt_K - Z_mean_k))


    # Auto Correlation Calculation (considering missing values)
    Autocor = Autocov / Autocov[0]
    
    return Autocor



def perform_EVF_computation(Z_ini, Z_act, lags, ind_MV, Wt0, ACP_cutoff, level_SD_perc, max_iter, it, test_end, first_iteration) :
        
    N = len(Z_ini)
    dt = 1.0
    ind_Z = np.arange(N)+1
    
    if first_iteration == False :
        Wt0 = np.ones(N)
    
    Autocor = autocorrelation(Z_act, lags, Wt0)
    
    # Confidence interval of the autocorrelation function
    rk = 1.96 / np.sqrt(N - lags - 1)
    ind_p = np.where(abs(Autocor) < rk)[0]

    # Number of lagged series for PCA
    p = ind_p[0] + 1 if len(ind_p) > 0 else 0

    # Exception to ignore the case of P = 0
    if p == 0 : 
        
        if first_iteration : 
            
            Z_ini[ind_MV] = np.nan
            return Z_ini
        
        else : 
            
            Z_it = np.copy(Z_ini)
            Z_it[ind_MV] = Z_act[ind_MV]
            
            test_val = np.sqrt(np.sum((Z_it[ind_MV] - Z_act[ind_MV]) ** 2)) / len(ind_MV)
            
            if test_val < np.std(Z_ini[np.isfinite(Z_ini)]) / level_SD_perc:
               test_end = True

            if it >= max_iter:
                Z_act[ind_MV] = np.nan  # NaN for values where convergence did not occur
                test_end = True
            
            return [Z_it, it, test_end]
            
        
    # P-lagged Matrix construction
    matrix = np.full((p, N - (p-1) * int(dt)), np.nan)
    matrix_ind_Z = np.copy(matrix)          

    # Fill lagged matrices
    index_max = N - (p-1) * int(dt)
    matrix[0, :] = Z_act[:index_max]
    matrix_ind_Z[0, :] = ind_Z[:index_max]
    
    M_weight = np.copy(matrix)
    M_weight[0, :] = Wt0[:index_max]

    for i in np.arange(1, p):
        
        index_max = N - (p -1 - i) * int(dt)
        matrix[i, :] = Z_act[i : index_max]
        matrix_ind_Z[i, :] = ind_Z[i : index_max]
        M_weight[i, :] = Wt0[i : index_max]

    if (first_iteration == False) and np.isnan(matrix).any() :
        sys.exit('There should be no missing data in the second iteration')
        
    matrix[ np.isnan(matrix) ] = np.nanmean(Z_act)
        
    matrix = np.flip(matrix, axis=1)
    matrix_ind_Z = np.flip(matrix_ind_Z, axis=1)
    M_weight = np.flip(M_weight, axis=1)

    # Calculation of the centered matrix
    matrix_std = np.full_like(matrix, np.nan)
    corr_mean = np.full(p, np.nan)

    for i in range(p):
        corr_mean[i] = np.nansum(matrix[i, :] * M_weight[i, :]) / np.nansum(M_weight[i, :])
        matrix_std[i, :] = matrix[i, :] - corr_mean[i]

    # Covariance matrix for PCA
    covar = np.full((p, p), np.nan)

    for i in range(p):
        for j in range(p):
            x1 = matrix_std[i, :]
            x2 = matrix_std[j, :]
            valid_indices = np.where(( np.isfinite(x1) ) & ( np.isfinite(x2) ))[0]
            covar[i, j] = np.cov(x1[valid_indices], x2[valid_indices])[0, 1]

    # Perform PCA
    eigen_val, eigen_vector = np.linalg.eigh(covar)
    eigen_val = np.flip(eigen_val)
    eigen_vector = np.flip(np.transpose(eigen_vector), axis = 1)
    exp_var = eigen_val / np.sum(eigen_val) * 100.0

    # (np.dot(covar, eigen_vector[:, 0]) - eigen_val[0] * eigen_vector[:, 0]) # Check if eigen_val and eigen_vectors objects are correct. The results hsould be zero.

    # Cumulative variance calculation
    cumul = np.cumsum(exp_var)
    Test_axes = np.where(cumul > ACP_cutoff)[0]
    N_axes = Test_axes[0] + 1 if len(Test_axes) > 0 else 1

    # Components calculation
    
    # Ctk = np.full((matrix_std.shape[1], len(eigen_vector)), -9999.0)
    # for t in np.arange(matrix_std.shape[1]) : 
    #     for k in np.arange(len(eigen_val)) : 
    #         a = []        
    #         for the_p in np.arange(p) : 
    #             a.append( M_weight[the_p, t] * eigen_vector[the_p, k] * matrix_std[the_p, t] )
    #         Ctk[t,k] = np.sum(a)    
                
    Ctk = np.dot( np.transpose(matrix_std * M_weight) , eigen_vector)
    Y_p_n_axes = np.full((p, Ctk.shape[1], N_axes), np.nan)

    for i in range(N_axes):
        Y_p_n_axes[:, :, i] = np.transpose(  np.dot( Ctk[i, :].reshape(-1, 1), eigen_vector[i, :].reshape(1, -1) ) )

    Y_p = np.nansum(Y_p_n_axes, axis=2) if N_axes > 1 else Y_p_n_axes[:, :, 0]
    Y_p += np.repeat(corr_mean.reshape(-1, 1), Y_p.shape[1], axis=1)

    if first_iteration : 
        # New vector creation after the first iteration
        Z_i = np.full(N, np.nan)
    else : 
        Z_i = np.arange(N).astype(float)
    
    for i in range(N):
        
        ind = np.where(matrix_ind_Z == (ind_Z)[i])[0]
        
        if first_iteration == False : 
            
            if len(ind) <= 0:
                sys.exit('Error: no valid indices found.')
                continue
                
            if np.isnan(np.mean(Y_p[ind])):
                sys.exit('Error: NaN values present.')
        
            Z_i[i] = np.nanmean(Y_p[ind])
        
        else : 
        
            Z_i[i] = np.nanmean(Y_p[ind]) if len(ind) > 0 else np.nan

    Z_it = np.copy(Z_ini)
    Z_it[ind_MV] = Z_i[ind_MV]
    
    if first_iteration : 
    
        return Z_it    
    
    else : 
        
        it += 1
        test_val = np.sqrt(np.nansum((Z_it[ind_MV] - Z_act[ind_MV]) ** 2)) / len(ind_MV)
        
        if test_val < np.nanstd(Z_ini) / level_SD_perc:
           test_end = True

        if it >= max_iter:
            Z_it[ind_MV] = np.nan  # NaN for values where convergence did not occur
            test_end = True
            
        return [Z_it, it, test_end]



def F_EVF_V1_2(X, level_SD_perc, ACP_cutoff, max_iter):
    
    """
    PURPOSE:
    Fill the gaps present in a time series using the method defined in Ibanez and Conversi (2002).
    
    INPUTS:
    X: Original time series with gaps
    level_SD_perc: Threshold decision value for deciding the end of the iteration
    ACP_cutoff: Cutoff value of the cumulated variance used for selecting the number of axes of the PCA
    max_iter: Maximum number of iterations
    
    OUTPUTS:
    The original time series with gaps replaced by estimated values
    """
    
    Z = np.copy(X)
    ind_MV = np.where( np.isnan(Z) )[0]  # Indices of missing values
    
    if len(ind_MV) == 0 :
        
        return Z
    
    # Initial weight setup: 1 where not missing, 0 where missing
    Wt0 = np.ones(len(Z))
    Wt0[ind_MV] = 0
        
    # Definition of the lag for autocorrelation
    lags = np.arange(len(Z) // 3)
    
    Z_act = perform_EVF_computation(Z_ini = Z, Z_act = Z, 
                                     lags = lags, ind_MV = ind_MV, Wt0 = Wt0, 
                                     ACP_cutoff = ACP_cutoff, level_SD_perc = level_SD_perc, max_iter = max_iter, 
                                     it = 0, test_end = False, first_iteration = True)
        
    
    test_end = False
    it = 0

    while test_end == False:
        
        [Z_act, it, test_end] = perform_EVF_computation(Z_ini = Z, Z_act = Z_act, 
                                                         lags = lags, ind_MV = ind_MV, Wt0 = np.ones(len(Z)), 
                                                         ACP_cutoff = ACP_cutoff, level_SD_perc = level_SD_perc, max_iter = max_iter, 
                                                         it = it, test_end = test_end, first_iteration = False)
    
    return Z_act

    
def F_census_X_11_pezzulli_V1_2(Xt, period):
        
    # Test if the period is even or odd
    if period % 2 != 0:
        test_period = 0  # period is odd
    else:
        test_period = 1  # period is even

    Xt = np.array(Xt)  # Reform the input into a NumPy array
    N = len(Xt)

    # Step 1
    if test_period == 1:
        Tt_1 = F_moving_average_v1_0(Xt, 2, period, period)
    else:
        Tt_1 = F_moving_average_v1_0(Xt, 1, period, period)

    Zt_1 = Xt - Tt_1

    SMA_22_1 = np.full(N, np.nan)

    for i in range(period):
        t_SMA_22_1 = F_moving_average_v1_0(Zt_1[i::period], 2, 2, 1)
        SMA_22_1[ np.arange(len(t_SMA_22_1)) * period + i ] = t_SMA_22_1

    # Normalisation
    if test_period == 1:
        St_1 = SMA_22_1 - F_moving_average_v1_0(SMA_22_1, 2, period, period)
    else:
        St_1 = SMA_22_1 - F_moving_average_v1_0(SMA_22_1, 1, period, period)

    # Step 2
    Yt_2 = Xt - St_1
    Tt_2 = F_henderson_v1_0(X= Yt_2, order = 2 * period - 1, period = period)

    Zt_2 = Xt - Tt_2

    SMA_22_2 = np.full(N, np.nan)

    for i in np.arange(period) :
        t_SMA_22_2 = F_moving_average_v1_0(Zt_2[i::period], 2, 2, 1)
        # index_to_fill = np.array( [x for x in np.arange(len(t_SMA_22_2)) * period + i if x <= len(SMA_22_2)] )
        index_to_fill = np.array( [x for x in np.arange(len(t_SMA_22_2)) * period + i] )
        SMA_22_2[ index_to_fill ] = t_SMA_22_2

    # Normalisation
    if test_period == 1:
        St_2 = SMA_22_2 - F_moving_average_v1_0(SMA_22_2, 2, period, period)
    else:
        St_2 = SMA_22_2 - F_moving_average_v1_0(SMA_22_2, 1, period, period)

    # Step 3
    Yt_3 = Xt - St_2
    Tt_3 = F_henderson_v1_0(Yt_3, 2 * period - 1, period)

    Yt = Xt - Tt_3 - St_2

    Tt = Tt_3
    St = St_2

    decomposed = {'Tt': Tt, 'St': St, 'Yt': Yt}

    return decomposed


def F_moving_average_v1_0(X, m, p, period):
        
    period_calc = period
    p_calc = p

    if p_calc % 2 == 0:
        p_calc += 1
    half_window = (p_calc-1)/2

    X = np.array(X)
    N = len(X)

    if period_calc == 1:
        if p == 2:
            period_calc = 1
        elif p == 3:
            period_calc = 2
        elif p == 5:
            period_calc = 3

    # Calculate indices for slicing
    start1 = int((period_calc - 1) - half_window + 1)
    end1 = int(period_calc)
    
    start2 = int(N - 1 - period_calc + 1)
    end2 = int(start2 + half_window)
    
    # Slice the array and concatenate
    X_calc = np.concatenate([
        X[start1:end1],
        X,
        X[start2:end2]
    ])

    N_calc = len(X_calc)
    MA_X = np.arange(N_calc).astype(float)

    if m % 2 == 0 :
        
        for i in np.arange(half_window, N_calc - (half_window) ) :
            
            MA_X[int(i)] = (m / 2 * X_calc[int(i - half_window)] + 
                       np.nansum(m * X_calc[ int( (i - half_window) + 1 ) : int(i - 1 + half_window +1)]) + 
                       m / 2 * X_calc[int(i + half_window)]) / ((p_calc - 1) * m)
            
        MA_X = MA_X[ int(half_window) : int(N + half_window) ]
        
    else :
        
        if m == 1 :
            
            MA_X = pd.Series(X_calc).rolling(window=p_calc, center=True, min_periods = 1).mean().to_numpy()
            MA_X = MA_X[ int(half_window) : int(N + half_window)]
            
        else:
            
            weight = np.concatenate([np.arange(2) + 1, np.full(p_calc - 2, m), np.arange(2)[::-1] + 1])
            X_calc = np.concatenate([X[int( (period_calc - 1) - (p_calc + 2 - 1) / 2 + 1) : int((period_calc - 1) +1)], 
                                     X, 
                                     X[ int(N - 1 - period_calc + 1) : int( (N - 1 - period_calc + 1) + ((p_calc + 2 - 1) / 2 - 1)) +1]])
            N_calc = len(X_calc)
            MA_X = np.zeros(N_calc)
            for i in range((p_calc + 1) / 2, N_calc - ((p_calc + 1) / 2)):
                MA_X[i] = np.sum(X_calc[ int(i - (p_calc + 1) / 2) : int( i + (p_calc + 1) / 2 +1 ) ] * weight) / np.sum(weight)
            MA_X = MA_X[ int(half_window +1) : int(N + half_window +1) ]

    return MA_X


def F_henderson_v1_0(X, order, period):
        
    half_order = (order-1)/2

    N = len(X)
    X_calc = np.concatenate([X[ int((period - 1) - half_order + 1) : int(period) ], 
                             X, 
                             X[ int(N - 1 - period + 1) : int((N -1 - period +1) + (half_order - 1) +1)]])
    N_calc = len(X_calc)
    MA_X = np.zeros(N_calc)

    if order % 2 == 0 :
        raise ValueError(f"ORDER MUST BE ODD {order}")

    i_H = np.concatenate([np.arange(-half_order, 0), np.array([0]), np.arange(half_order) + 1])
    n_H = float((len(i_H) - 1) / 2 + 2)
    h = (315 * ((n_H - 1) ** 2 - i_H ** 2) * (n_H ** 2 - i_H ** 2) * ((n_H + 1) ** 2 - i_H ** 2) * 
         (3 * n_H ** 2 - 16 - 11 * i_H ** 2) / (8 * n_H * (n_H ** 2 - 1) * (4 * n_H ** 2 - 1) * 
                                                (4 * n_H ** 2 - 9) * (4 * n_H ** 2 - 25)))

    for i in range( int(half_order), int(N_calc - (half_order) )):
        MA_X[i] = np.sum( X_calc[ int(i - half_order) : int(i + half_order + 1)] * h) / np.sum(h)
    
    MA_X = MA_X[ int(half_order) : int(N + half_order) ]

    return MA_X







def f_test_season_v1_1(X, T, S, period, test_on_X=False, test_on_S=False, n_years=None):
        
    N_years = len(X) // period

    RES_TEST = np.nan
    S_evolutive = 0.0

    S_I = X - T
    if test_on_S:
        S_I = S

    # Create arrays for the calculations
    M_S_I = np.full((period, N_years), np.nan)
    M_X = np.full((period, N_years), np.nan)

    # Populate the arrays
    for i in range(period):
        M_S_I[i, :len(X[i::period])] = S_I[i::period]
        M_X[i, :len(X[i::period])] = X[i::period]

    ind_bad = np.where( np.isnan(M_S_I) )
    ind_valid = np.where( np.isfinite(M_S_I) )

    N_col_M_S_I, N_lign_M_S_I = M_S_I.shape

    N_col_M_X, N_lign_M_X = M_X.shape

    if test_on_X:
        
        # Test on X
        S_A2_X_i = np.full(N_col_M_X, -999.0)
        Mean_X_i = np.full(N_col_M_X, -999.0)
        
        for i in np.arange(N_col_M_X) :
            count_n_i = np.sum( np.isfinite( M_X[i, :]) )
            if count_n_i > 0:
                Mean_X_i[i] = np.nanmean(M_X[i, :])
                S_A2_X_i[i] = count_n_i * (Mean_X_i[i] - np.nanmean(M_X[ind_valid]))**2

        S_R2_X = np.nansum((M_X - np.tile(Mean_X_i[i], (N_lign_M_X, 1)).T)**2)
        S_A2_X = np.nansum(S_A2_X_i)

        Fs = (S_A2_X / (N_col_M_X - 1)) / (S_R2_X / (len(ind_valid[0]) - N_col_M_X))
        
        dfn = N_col_M_X - 1
        dfd = len(ind_valid[0]) - N_col_M_X

        Limit_FS = 0.1 / 100
        test_Fs = 1 - f.cdf(Fs, dfn, dfd)

    else:
        
        # Test on M_S_I
        S_A2_X_i = np.full(N_col_M_S_I, -999.0)
        Mean_MSI_i = np.full(N_col_M_S_I, -999.0)

        for i in range(N_col_M_S_I):
            
            count_n_i = np.where( np.isfinite(M_S_I[i, :]) )[0]
            if len(count_n_i) > 0:
                mean_MSI_i = np.nanmean(M_S_I[i, :])
                S_A2_X_i[i] = len(count_n_i) * (mean_MSI_i - np.nanmean(M_S_I[ind_valid]))**2

        S_R2_X = np.nansum((M_S_I - np.tile(Mean_MSI_i, (N_lign_M_S_I, 1)).T)**2)
        S_A2_X = np.nansum(S_A2_X_i)

        Fs = (S_A2_X / (N_col_M_S_I - 1)) / (S_R2_X / (len(ind_valid[0]) - N_col_M_S_I))
        dfn = N_col_M_S_I - 1
        dfd = len(ind_valid[0]) - N_col_M_S_I

        Limit_FS = 0.1 / 100
        test_Fs = 1 - f.cdf(Fs, dfn, dfd)

    if test_Fs < Limit_FS :
        
        abs_M_S_I = np.abs(M_S_I - np.nanmean(M_S_I[ind_valid]))

        if len(ind_bad[0]) != 0:
            abs_M_S_I[ind_bad] = np.nan

        # S_2 = np.nanvar(abs_M_S_I) * len(ind_valid[0])
        X_p_p = np.nanmean(abs_M_S_I)

        X_i_p = np.arange(N_col_M_S_I)
        X_p_j = np.arange(N_lign_M_S_I)

        for i in range(N_col_M_S_I):
            valid_indices = np.isfinite(abs_M_S_I[i, :])
            if np.any(valid_indices):
                X_i_p[i] = np.nanmean(abs_M_S_I[i, valid_indices])

        for j in range(N_lign_M_S_I):
            valid_indices = np.isfinite(abs_M_S_I[:, j])
            if np.any(valid_indices):
                X_p_j[j] = np.nanmean(abs_M_S_I[valid_indices, j])

        # S_A_2 = N_lign_M_S_I * np.sum((X_i_p - X_p_p)**2)
        S_B_2 = N_col_M_S_I * np.sum((X_p_j - X_p_p)**2)
        S_R_2 = np.sum((abs_M_S_I - np.tile(X_i_p, (N_lign_M_S_I, 1)).T - np.tile(X_p_j, (N_col_M_S_I, 1)) + X_p_p)**2)

        FM = (S_B_2 / (N_lign_M_S_I - 1)) / (S_R_2 / ((N_lign_M_S_I - 1) * (N_col_M_S_I - 1)))
        dfn = N_lign_M_S_I - 1
        dfd = (N_col_M_S_I - 1) * (N_lign_M_S_I - 1)

        Limit_Fm = 5.0 / 100
        test_Fm = 1 - f.cdf(FM, dfn, dfd)

        T1 = 7.0 / Fs
        T2 = 3.0 * FM / Fs
        T = np.sqrt(0.5 * (T1 + T2))

        if test_Fm > Limit_Fm:
            if T1 >= 1 or T2 >= 1:
                RES_TEST = 1.0
            else:
                Limit_KW = 0.001
                test_KW_M_X = kw_TEST(M_S_I, df=period - 1, MISSING=-1.0)

                if test_KW_M_X[1] > Limit_KW:
                    RES_TEST = 1.0
                else:
                    RES_TEST = 2.0
        else:
            if T >= 1:
                RES_TEST = 0.0
            else:
                S_evolutive = 1.0

                if T1 >= 1 or T2 >= 1:
                    RES_TEST = 1.0
                else:
                    Limit_KW = 0.1 / 100
                    test_KW_M_X = kw_TEST(M_S_I, df=period - 1, MISSING=-1.0)

                    if test_KW_M_X[1] > Limit_KW:
                        RES_TEST = 1.0
                    else:
                        RES_TEST = 2.0

    else:
        RES_TEST = 0.0

    if RES_TEST != 2.0:
        res_test_season = RES_TEST
    else:
        if S_evolutive == 0:
            res_test_season = 2.0
        else:
            res_test_season = 3.0

    # return res_test_season, S_evolutive
    return res_test_season == 0.0


def kw_TEST(data, df, MISSING):
        
    """
    Perform the Kruskal-Wallis H-test for independent samples.

    Parameters:
    - data: 2D array-like, where each row represents a group (season).
    - df: Degrees of freedom (number of groups - 1).
    - MISSING: Value used to denote missing data.

    Returns:
    - H: Kruskal-Wallis H statistic.
    - p_value: p-value for the test.
    """
        
    
    # Replace missing values with NaN
    data = np.where(data == MISSING, np.nan, data)

    # Split data into groups
    groups = [data[i, ~np.isnan(data[i, :])] for i in range(data.shape[0])]

    # Perform Kruskal-Wallis test
    H, p_value = kruskal(*groups)

    return [H, p_value]





def f_kendall_sen_v1_0(X, alpha=0.05, correction_tie=False):
    """
    Test for monotonic change in a time series using the non-parametric Kendall statistics
    and compute the Sen slope estimator.
    
    Parameters:
        X (array-like): Input time series vector
        hbad (float, optional): Bad value definition, default is -9999
        alpha (float, optional): Significance level, default is 0.05
        correction_tie (int, optional): If set, correct the calculation for the presence of tied data, default is 0
    
    Returns:
        tuple: (sen_slope, prob, res_kendall_test, b, RC)
    """
    
    output = {
        'Is_the_change_of_annual_values_monotonic': np.nan,
        'Sen_slope': np.nan,
        'Sen_intercept': np.nan,
        'Rate_of_change_of_annual_values_in_percentage_per_time': np.nan
    }
    
    # Prepare data
    dat_in = np.array(X)
    v_indice = np.arange(1, len(dat_in) + 1)
    ind_dat_in_valid = np.where( np.isfinite(dat_in) )[0]
    # count_dat_in_valid = len(ind_dat_in_valid)
    
    n = len(dat_in[ind_dat_in_valid])
    m_dat_in = np.full((n, n), np.nan).astype(float)
    
    # Compute the matrix of slopes
    for i in range(n):
        index_i = ind_dat_in_valid[i]
        for j in range(n):
            index_j = ind_dat_in_valid[j]
            m_dat_in[i, j] = (dat_in[index_i] - dat_in[ind_dat_in_valid[index_j]]) / \
                             (v_indice[index_i] - v_indice[ind_dat_in_valid[index_j]])
    
    # Handle invalid data by setting diagonals and upper triangle to hbad
    for i in range(n):
        m_dat_in[i, i:n] = np.nan
    
    r_m_dat_in = m_dat_in[1:, :-1]
    
    # Identify valid data
    valid_data = np.where( np.isfinite(r_m_dat_in) )
    sign_m_dat_in = np.sign(r_m_dat_in)
        
    # Handle ties
    if correction_tie:
        test_tie = f_count_tie(dat_in[ind_dat_in_valid])
        if len(test_tie) > 1:
            corr_factor_tie = 0
            for i in range(len(test_tie[0])):
                corr_factor_tie += test_tie[1][i] * (test_tie[0][i] * (test_tie[0][i] - 1) * (2 * test_tie[0][i] + 5))
        else:
            corr_factor_tie = 0
    else:
        corr_factor_tie = 0
    
    # Sum of signs
    Sum_sign = np.nansum(sign_m_dat_in[valid_data])
    var_S = 1. / 18. * ((n * (n - 1) * (2 * n + 5)) - corr_factor_tie)
    
    # Compute Sen's slope
    sen_slope = np.nanmedian(r_m_dat_in[valid_data])
    
    # Compute intercept
    b_mat = dat_in[ind_dat_in_valid] - (sen_slope * v_indice[ind_dat_in_valid])
    b = np.nanmedian(b_mat)
    
    # Rate of Change (%)
    RC = (sen_slope / abs(b)) * 100
    
    # Z-values calculation
    # Z_tab = norm.ppf(1 - alpha / 2)
    
    if Sum_sign > 0:
        Z_calc = (Sum_sign - 1) / np.sqrt(var_S)
    elif Sum_sign < 0:
        Z_calc = (Sum_sign + 1) / np.sqrt(var_S)
    else:
        Z_calc = 0
    
    # C_alpha = Z_tab * np.sqrt(var_S)
    
    # M1 = (len(valid_data[0]) - C_alpha) / 2
    # M2 = (len(valid_data[0]) + C_alpha) / 2
    
    # sorted_slope = np.sort(r_m_dat_in[valid_data])
    
    # if M1 < 1 :
    #     Lower_limit = sorted_slope[0]
    # else :
    #     Lower_limit = sorted_slope[int(np.round(M1 -1))]
    
    # if M2 > (len(sorted_slope)-1) :
    #     Upper_limit = sorted_slope[-1]
    # else : 
    #     Upper_limit = sorted_slope[int(np.round(M2))]
    
    prob = abs((norm.cdf(abs(Z_calc)) - 1) * 2)
    
    if prob < alpha:
        res_kendall_test = 1
    else:
        res_kendall_test = 0
        
    # return sen_slope, prob, res_kendall_test, b, RC
    
    output['Is_the_change_of_annual_values_monotonic'] = res_kendall_test == 1
    output['Sen_slope'] = sen_slope
    output['Sen_intercept'] = b
    output['Rate_of_change_of_annual_values_in_percentage_per_time'] = RC
    
    return output



def f_count_tie(Z):
    """
    Counts the number of ties in the input series.

    Parameters:
        Z (array-like): Input vector

    Returns:
        numpy.ndarray: If tie values exist, returns a matrix with types of ties (double, triple, etc.) and their occurrences.
                       Else, returns 0.
    """
    
    X = np.array(Z)
    N = len(X)
    tie_count = np.zeros(N + 1, dtype=int)

    # Loop through each element in the array
    for i in range(N):
        # Find indices where elements are equal to X[i]
        i_tie = np.where(X == X[i])[0]
        count_a = len(i_tie)
        
        # If there are ties (more than one occurrence)
        if count_a > 1:
            if np.isfinite(X[i]) :
                tie_count[count_a] += 1
            
            # Mark the counted ties with -9999
            X[i_tie] = np.nan

    # If there are any ties, prepare the result array
    if np.sum(tie_count) != 0:
        result_tie = np.transpose([
            np.where(tie_count >= 1)[0],  # Types of ties
            tie_count[tie_count >= 1]     # Corresponding occurrences
        ])
    else:
        result_tie = np.array([0])

    return result_tie



def f_seasonal_kendall_v2_0(X, period, alpha=0.05, correction_tie=0):
    """
    Performs the non-parametrical Seasonal Kendall test for trends.
    Computes the Sen slope for the whole time series and the monthly slope.
    Tests the homogeneity of the seasonal trends.
    Provides the significance level of the seasonal and global trends.

    Parameters:
        X (array): Input time series.
        period (int): The periodicity of the time series.
        alpha (float, optional): Significance level to be used. Default is 0.05.
        hbad (float, optional): Bad value. Default is -9999.
        correction_tie (int, optional): Correct the calculation for ties in the data. Default is 0.

    Returns:
        dict: Results containing the seasonal Kendall test information.
    """
        
    # Initialize the output structure
    # v_bad = -9999.
    # out_seasonal_kendall = {
    #     'sen_slope': np.nan,
    #     'prob_sen_slope': np.nan,
    #     'prob_chi_homog': np.nan,
    #     'prob_chi_trend': np.nan,
    #     'conf_limit': [np.nan, np.nan],
    #     'monthly_sen_slope': np.full(period, np.nan),
    #     'prob_monthly': np.full(period, np.nan),
    #     'intercept': np.nan,
    #     'RC': np.nan
    # }
    
    out_seasonal_kendall = {
        'Sen_slope': np.nan,
        'prob_Sen_slope' : np.nan,
        'Sen_intercept': np.nan,
        
        'Sen_slope_confidence_limit': [np.nan, np.nan],
        'Is_Sen_slope_significant': np.nan,
        
        'Sen_monthly_slopes': np.full(period, np.nan),
        'prob_Sen_monthly_slopes': np.full(period, np.nan),
        'prob_test_Homogeneity_of_the_seasonal_signal_across_years': np.nan,
        'prob_test_Trend_in_the_seasonal_signal_across_years': np.nan,
        'Homogeneity_of_the_seasonal_signal_across_years': True,
        'Significant_Trend_in_the_seasonal_signal_across_years' : False,

        'Rate_of_change_in_percentage_per_time': np.nan
    }
    
    # Prepare data
    X = np.array(X)
    ind_dat_in_valid = np.where( np.isfinite(X) )[0]
    n_val = len(ind_dat_in_valid)
    dat_in = np.copy(X)
    
    # Handle missing values
    ind_bad_val_in = np.where( np.isnan(X) )[0]
    if len(ind_bad_val_in) > 0:
        dat_in[ind_bad_val_in] = np.nan
    
    # Reshape data
    n_lign = int(np.ceil(len(dat_in) / period))
    if n_lign * period != len(dat_in):
        complete_y = np.full((n_lign * period) - len(dat_in), np.nan)
        m_dat_in = np.reshape(np.concatenate((dat_in, complete_y)), (period, n_lign), order='F')
    else:
        m_dat_in = np.reshape(dat_in, (period, n_lign), order='F')

    # Initialize arrays
    S_var = np.full(period, np.nan)
    n_month = np.full(period, np.nan)
    Sign_Diff = np.full(period, np.nan)
    sen_slope_m = np.full((period, int(np.sum(np.arange(n_lign)))), np.nan)
    
    # Calculate Sen's slope for each month
    for i_month in range(period):
        
        i_valid = np.where( np.isfinite(m_dat_in[i_month, :]) )[0]
        count_valid = len(i_valid)
        
        if count_valid > 0:
            
            ind = 0
            sub_diff_k = np.full(int(np.sum(np.arange(count_valid))), np.nan)
            sub_slope_s = np.full(int(np.sum(np.arange(count_valid))), np.nan)
            
            for i in range(count_valid - 1):
                for j in range(i + 1, count_valid):
                    sub_diff_k[ind] = (m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]]) / \
                                      abs(m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]])
                    sub_slope_s[ind] = (m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]]) / (j - i)
                    if sub_slope_s[ind] == 0:
                        sub_diff_k[ind] = 0  # Correct for differences of 0
                    ind += 1
            
            sen_slope_m[i_month, :len(sub_slope_s)] = sub_slope_s
            
            # Correction for ties in the data
            if correction_tie == 1:
                test_tie = np.unique(m_dat_in[i_month, i_valid], return_counts=True)
                if len(test_tie[0]) > 1:
                    corr_factor_tie = 0
                    for i in range(len(test_tie[0])):
                        corr_factor_tie += test_tie[1][i] * (test_tie[0][i] * (test_tie[0][i] - 1) * (2 * test_tie[0][i] + 5))
                    S_var[i_month] = (count_valid * (count_valid - 1) * (2 * count_valid + 5) - corr_factor_tie) / 18
                else:
                    S_var[i_month] = count_valid * (count_valid - 1) * (2 * count_valid + 5) / 18
            else:
                S_var[i_month] = count_valid * (count_valid - 1) * (2 * count_valid + 5) / 18
            
            Sign_Diff[i_month] = np.nansum(sub_diff_k)
            n_month[i_month] = count_valid

    # Test trend homogeneity over seasons
    Z_calc_i_month = np.full(period, np.nan)
    
    for i in range(period):
        if np.isfinite(Sign_Diff[i]) :
            if Sign_Diff[i] > 0:
                Z_calc_i_month[i] = Sign_Diff[i] / np.sqrt(S_var[i])
            elif Sign_Diff[i] < 0:
                Z_calc_i_month[i] = Sign_Diff[i] / np.sqrt(S_var[i])
            else:
                Z_calc_i_month[i] = 0
    
    i_m_valid = np.where( np.isfinite( Z_calc_i_month ) )[0]
    count_n_valid_months = len(i_m_valid)
    
    Chi_total = np.nansum(Z_calc_i_month ** 2)
    Chi_trend = count_n_valid_months * (np.nanmean(Z_calc_i_month)) ** 2
    Chi_homo = Chi_total - Chi_trend
    Chi_theo_homog = chi2.ppf(0.95, count_n_valid_months - 1)
    
    out_seasonal_kendall['prob_test_Homogeneity_of_the_seasonal_signal_across_years'] = 1 - chi2.cdf(Chi_homo, count_n_valid_months - 1)
    
    if abs(Chi_homo) < Chi_theo_homog:
        # homog = 1
        out_seasonal_kendall['Homogeneity_of_the_seasonal_signal_across_years'] = False
                
        if abs(Chi_trend) < chi2.ppf(0.95, 1):
            # overall_chi_test = 1
            out_seasonal_kendall['Significant_Trend_in_the_seasonal_signal_across_years'] = True

        out_seasonal_kendall['prob_test_Trend_in_the_seasonal_signal_across_years'] = 1 - chi2.cdf(Chi_trend, 1)
    
    # Monthly trend analyses
    out_seasonal_kendall['Sen_monthly_slopes'] = np.full(period, np.nan)
    
    for i in range(period):
        i_val_i_month = np.where( np.isfinite(sen_slope_m[i, :]) )[0]
        i_count_valid_slope = len(i_val_i_month)
        
        if i_count_valid_slope > 0:
            i_sel_slope = sen_slope_m[i, i_val_i_month]
            if i_count_valid_slope / 2 - round(i_count_valid_slope / 2) == 0:
                out_seasonal_kendall['Sen_monthly_slopes'][i] = np.nanmedian(i_sel_slope)
            else:
                out_seasonal_kendall['Sen_monthly_slopes'][i] = np.nanmedian(i_sel_slope)
    
    Z_calc_monthly = np.full(period, np.nan)
    out_seasonal_kendall['prob_Sen_monthly_slopes'] = np.full(period, np.nan)
    
    for i in range(period):
        if np.isfinite( Sign_Diff[i] ) :
            if Sign_Diff[i] > 0:
                Z_calc_monthly[i] = (Sign_Diff[i] - 1) / np.sqrt(S_var[i])
            elif Sign_Diff[i] < 0:
                Z_calc_monthly[i] = (Sign_Diff[i] + 1) / np.sqrt(S_var[i])
            else:
                Z_calc_monthly[i] = 0
    
    out_seasonal_kendall['prob_Sen_monthly_slopes'][i_m_valid] = np.abs((norm.cdf(np.abs(Z_calc_monthly[i_m_valid])) - 1) * 2)
    
    # Sen's slope
    valid_sen_slope = np.where( np.isfinite(sen_slope_m) )[0]
    count_valid_slope = len(valid_sen_slope)
    selec_slope = sen_slope_m.flatten()[valid_sen_slope]
    
    if count_valid_slope / 2 - round(count_valid_slope / 2) == 0:
        sen_slope = np.nanmedian(selec_slope)
    else:
        sen_slope = np.nanmedian(selec_slope)
    
    out_seasonal_kendall['Sen_slope'] = sen_slope / period
    
    # Intercept
    b_mat = dat_in[ind_dat_in_valid] - (sen_slope / period * ind_dat_in_valid)
    if n_val / 2 - round(n_val / 2) == 0:
        out_seasonal_kendall['Sen_intercept'] = np.nanmedian(b_mat)
    else:
        out_seasonal_kendall['Sen_intercept'] = np.nanmedian(b_mat)
    
    # Rate of Change (%)
    out_seasonal_kendall['Rate_of_change_in_percentage_per_time'] = (sen_slope / np.abs(out_seasonal_kendall['Sen_intercept'])) * 100
    
    # Significance of overall trend
    total_sign_diff = np.sum(Sign_Diff[i_m_valid])
    total_s_var = np.sum(S_var[i_m_valid])
    
    if total_sign_diff > 0:
        Z_calc = (total_sign_diff - 1) / np.sqrt(total_s_var)
    elif total_sign_diff < 0:
        Z_calc = (total_sign_diff + 1) / np.sqrt(total_s_var)
    else:
        Z_calc = 0
    
    Z_tab = norm.ppf(1 - alpha / 2)
    C_alpha = Z_tab * np.sqrt(total_s_var)
    
    M1 = (count_valid_slope - C_alpha) / 2
    M2 = (count_valid_slope + C_alpha) / 2
    
    sorted_slope = np.sort(selec_slope)
    
    Lower_limit = sorted_slope[int(np.round(M1 - 1))]
    Upper_limit = sorted_slope[int(np.round(M2))]
    
    out_seasonal_kendall['Sen_slope_confidence_limit'] = [Lower_limit, Upper_limit]
    out_seasonal_kendall['prob_Sen_slope'] = np.abs((norm.cdf(np.abs(Z_calc)) - 1) * 2)
    
    # Kendall test result
    if np.abs(Z_calc) > np.abs(Z_tab):
        season_kendall_test = 1
    else:
        season_kendall_test = 0
        
    out_seasonal_kendall['Is_Sen_slope_significant'] = season_kendall_test == 1
    
    return out_seasonal_kendall



def F_test_pente(x, y, slope):
    """
    Function to test if the slope is significantly different from zero.

    Parameters:
    x : array-like
        Input array for the x-values.
    y : array-like
        Input array for the y-values.
    slope : float
        The slope to test.

    Returns:
    tuple
        decision : int
            1 if slope is significantly different from zero, otherwise 0.
        t_a : float
            The t-value calculated for the slope.
    """

    N = len(x)  # Number of elements in x

    # Calculate the variance of the slope
    s_2_b = ((np.std(y) / np.std(x))**2 - slope**2) / (N - 2)

    if s_2_b >= 0:
        # Calculate the t-value
        t_a = slope / np.sqrt(s_2_b)

        # Threshold value from t-distribution for alpha = 0.05 (two-tailed test)
        thresh = t.ppf(1 - 0.025, N - 2)

        if abs(t_a) >= thresh:
            # Slope is significantly different from 0
            decision = 1
        else:
            # Slope is not significantly different from 0
            decision = 0
    else:
        # If s_2_b is negative, set default values
        decision = 0
        t_a = 0

    return decision, t_a

def linfit(x, y):
    """
    Perform a linear fit.
    
    Parameters:
    x : array-like
        Independent variable.
    y : array-like
        Dependent variable.
    
    Returns:
    tuple
        Slope and intercept of the linear fit.
    """
    
    slope, intercept, _, _, _ = stats.linregress(x, y)
    return intercept, slope

