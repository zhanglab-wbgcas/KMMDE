# !/usr/bin/env python
# coding: utf-8
from __future__ import division

import os
import pandas as pd
import time
import sys
import numpy as np

from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.gaussian_process.kernels import RBF

import warnings
import numpy as np
from sys import stdout
import tst_est
from scipy.interpolate import interp1d
import statsmodels.stats.multitest as smm
from sklearn.metrics import pairwise_distances
warnings.filterwarnings('ignore')
from statsmodels.stats.multitest import multipletests

path = '/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/'

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
    """Compute the bootstrap null-distribution of MMD2u.
    """
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
    """Compute MMD^2_u, its null distribution and the p-value of the
    kernel two-sample test.

    Note that extra parameters captured by **kwargs will be passed to
    pairwise_kernels() as kernel parameters. E.g. if
    kernel_two_sample_test(..., kernel_function='rbf', gamma=0.1),
    then this will result in getting the kernel through
    kernel_function(metric='rbf', gamma=0.1).
    """
    m = len(X)
    n = len(Y)
    XY = np.vstack([X, Y])
    K = pairwise_kernels(XY, metric=kernel_function, **kwargs)
    # print(K)
    mmd2u = MMD2u(K, m, n)
    # print('mmd2u',mmd2u)
    if verbose:
        print("MMD^2_u = %s" % mmd2u)
        print("Computing the null distribution.")

    mmd2u_null = compute_null_distribution(K, m, n, iterations,
                                           verbose=verbose,
                                           random_state=random_state)
    p_value = max(1.0 / iterations, (mmd2u_null > mmd2u).sum() /
                  float(iterations))

    # print('num',(mmd2u_null > mmd2u).sum())
    if verbose:
        print("p-value ~= %s \t (resolution : %s)" % (p_value, 1.0 / iterations))

    return mmd2u, mmd2u_null, p_value





dir = path + "norm_grape_both.csv"
print('input data', "norm_grape_both.csv")

data = pd.read_csv(dir, index_col=0)
#data = data.head(100)
print(data.head());print(data.shape)
# print(data.index)
gene_list = data.index

tp=6
rep=2
kernel = 'rbf'

x = np.arange(1, tp*rep + 1).tolist()  # 
x_fine = np.round(np.arange(1, tp * rep + 0.001, 0.08), 2)
# print(x);
# print(x_fine)

p_value_list = []
MMD_list = []
for index, row in data.iterrows():
    # 
    # print(index)
    gene = index
    y1 = np.array((row[:(tp * rep)]).tolist())
    y2 = np.array((row[(-tp * rep):]).tolist())
    interp_func1 = interp1d(x, y1, kind='linear')
    interp_func2 = interp1d(x, y2, kind='linear')

    y_fine1 = interp_func1(x_fine)
    y_fine2 = interp_func2(x_fine)
    # print(y_fine1)

    # print(x)
    X = np.column_stack((x_fine, y_fine1))
    Y = np.column_stack((x_fine, y_fine2))
    
    sigma2 = np.median(pairwise_distances(X, Y, metric='euclidean')) ** 2
    gamma = 1.0 / sigma2
  
  
    mmd2u, mmd2u_null, p_value= kernel_two_sample_test(X, Y,kernel_function=kernel,
                                                          gamma=gamma,
                                                             verbose=False)


    MMD_list.append(mmd2u)
    p_value_list.append(p_value)


rank_df = pd.DataFrame({'gene': gene_list,'p_value_list':p_value_list,'MMD_list':MMD_list})

print(rank_df.head())
rank_df.to_csv(path + "25kmmde_grape_res.csv", index=False)

print('grape by 25kmmde func complete')





