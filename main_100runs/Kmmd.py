# !/usr/bin/env python
# coding: utf-8
from __future__ import division
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_kernels
import math
import warnings
import numpy as np
from sys import stdout
from scipy.interpolate import interp1d
from sklearn.metrics import pairwise_distances
warnings.filterwarnings("ignore")
from joblib import Parallel, delayed

path = '/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/'

def MMD2u(K, m, n):
    """The MMD^2_u unbiased statistic.
    """
    Kx = K[:m, :m]
    Ky = K[m:, m:]
    Kxy = K[:m, m:]
    return 1.0 / (m * (m - 1.0)) * (Kx.sum() - Kx.diagonal().sum()) + \
           1.0 / (n * (n - 1.0)) * (Ky.sum() - Ky.diagonal().sum()) - \
           2.0 / (m * n) * Kxy.sum()



def compute_null_distribution(K, m, n, iterations=1000, verbose=False,
                              random_state=None, marker_interval=1000):
    
    if type(random_state) == type(np.random.RandomState()):
        rng = random_state
    else:
        rng = np.random.RandomState(random_state)

    mmd2u_null = np.zeros(iterations)
    for i in range(iterations):
        if verbose and (i % marker_interval) == 0:
            # print(i),
            stdout.flush()
        idx = rng.permutation(m + n)
        K_i = K[idx, idx[:, None]]
        mmd2u_null[i] = MMD2u(K_i, m, n)
        # print(mmd2u_null[i])

    if verbose:
        print("")

    return mmd2u_null

def kernel_two_sample_test(X, Y, kernel_function, iterations=1000,
                           verbose=False, random_state=None, **kwargs):
    
    m = len(X)
    n = len(Y)
    XY = np.vstack([X, Y])
    K = pairwise_kernels(XY, metric=kernel_function, **kwargs)
    mmd2u = MMD2u(K, m, n)
    if verbose:
        print("MMD^2_u = %s" % mmd2u)
        print("Computing the null distribution.")

    mmd2u_null = compute_null_distribution(K, m, n, iterations,
                                           verbose=verbose,
                                           random_state=random_state)
    p_value = max(1.0 / iterations, (mmd2u_null > mmd2u).sum() /
                  float(iterations))

    if verbose:
        print("p-value ~= %s \t (resolution : %s)" % (p_value, 1.0 / iterations))

    return mmd2u, mmd2u_null, p_value


def calculate_seq(x_length):
    if x_length <= 4:
        return 0.4
    elif x_length <= 10:
        return 0.6 - 0.1 * math.ceil(x_length / 2)
    elif x_length <= 26:
        return 0.09 - 0.01 *  math.ceil((x_length - 10) / 2)
    elif x_length <= 40:
        return 0.009 - 0.001 * math.ceil((x_length - 26) / 2) 
    else:
        return 0.002



def process_row( row, tp, rep, x_fine,x, kernel_function, metric):
    y1 = row[:(tp * rep)].values
    y2 = row[(-tp * rep):].values
    interp_func1 = interp1d(x, y1, kind='linear')
    interp_func2 = interp1d(x, y2, kind='linear')

    y_fine1 = interp_func1(x_fine)
    y_fine2 = interp_func2(x_fine)

    X = np.column_stack((x_fine, y_fine1))
    Y = np.column_stack((x_fine, y_fine2))

    sigma2 = np.median(pairwise_distances(X, Y, metric=metric)) ** 2
    gamma = 1.0 / sigma2

    mmd2u2, mmd2u_null2, p_value2 = kernel_two_sample_test(X, Y, kernel_function=kernel_function, gamma=gamma, verbose=False)
    return ( mmd2u2, p_value2)

def kd_func(tp, rep,  normdata,  metric):
    x = np.arange(1, tp * rep + 1)
    seq = calculate_seq(len(x))
    x_fine = np.round(np.arange(1, tp * rep + 0.001, seq), 2)
    kernel_function = 'rbf'  # 

    results = Parallel(n_jobs=-1)(delayed(process_row)(row, tp, rep, x_fine,x, kernel_function, metric) for index, row in normdata.iterrows())
    # results.sort(key=lambda x: x[0])  
    MMD_list2, p_value_list2 = zip(*results)

    rank_df = pd.DataFrame({'gene': ["gene_" + str(i + 1) for i in range(normdata.shape[0])], 'MMD_list2': MMD_list2, 'p_value_list': p_value_list2})

    print(rank_df.head())
    print('ks func complete')
    return rank_df
