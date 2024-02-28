# !/usr/bin/env python
# coding: utf-8
import os
import pandas as pd
import time
import sys
import numpy as np

#读取txt文件 并去除重复
# txt_dir = '/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/blast_coldgo_gene/gene_ids.txt'
txt_dir = '/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/blast_coldgo_gene/01allcold_1e-10_50gene_ids.txt'
#读取kmmde的结果 以及特异结果
kmmde_dir = '/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/25kmmde_grape_res.csv'
kmmde_spec_dir ='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/KMMDE_spec.csv'
edger_dir ='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/grape_edger_res.csv'
edger_spec_dir='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/edgeR_spec.csv'
deseq2_dir ='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/grape_deseq2_res.csv'
deseq2_spec_dir='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/DEseq_spec.csv'
funpat_dir ='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/grape_funpat_res.csv'
funpat_spec_dir='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/Funpat_spec.csv'
impulsede2_dir ='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/grape_impulsede2_res.csv'
impulsede2_spec_dir='/home/liukc/project/sdc/02func_algori/revise_dre/grapedata_pro/Impls_spec.csv'

with open(txt_dir, "r") as f:
    txt_dat = f.readlines()
txt_dat = pd.read_csv(txt_dir, sep='\t', header=None)
print(len)

txt_dat.columns= ['gene']
print(txt_dat)
#去除重复
txt_list = txt_dat['gene']
org_len = len(txt_list)
print(org_len)
# using naive method
# to remove duplicated
# from list
txt_rm = []
for i in txt_list:
    if i not in txt_rm:
        txt_rm.append(i)

after_len = len(txt_rm)
print(after_len)

#读取kmmde的结果 以及特异结果
kmmde_dat = pd.read_csv(kmmde_dir)
# print(kmmde_dat.head())

kmmde_spec_dat =pd.read_csv(kmmde_spec_dir,encoding='gb18030',index_col = 0)
# print(kmmde_spec_dat.head())
kmmde_spec_dat.columns= ['gene']

#看交集
kmmde_cold = list(set(kmmde_spec_dat['gene']) & set(txt_dat['gene']))
print(len(kmmde_cold))

edger_spec_dat= pd.read_csv(edger_spec_dir,encoding='gb18030',index_col = 0)
edger_spec_dat.columns= ['gene']
edger_cold = list(set(edger_spec_dat['gene']) & set(txt_dat['gene']))
print(len(edger_cold))

deseq2_spec_dat= pd.read_csv(deseq2_spec_dir,encoding='gb18030',index_col = 0)
deseq2_spec_dat.columns= ['gene']
deseq2_cold = list(set(deseq2_spec_dat['gene']) & set(txt_dat['gene']))
print(len(deseq2_cold))

funpat_spec_dat= pd.read_csv(funpat_spec_dir,encoding='gb18030',index_col = 0)
funpat_spec_dat.columns= ['gene']
funpat_cold = list(set(funpat_spec_dat['gene']) & set(txt_dat['gene']))
print(len(funpat_cold))

impulsede2_spec_dat = pd.read_csv(impulsede2_spec_dir,encoding='gb18030',index_col = 0)
impulsede2_spec_dat.columns= ['gene']

impulsede2_cold = list(set(impulsede2_spec_dat['gene']) & set(txt_dat['gene']))
print(len(impulsede2_cold))

# edger_dat =
# deseq2_dat =
# funpat_dat =
# impulsede2_dat =


