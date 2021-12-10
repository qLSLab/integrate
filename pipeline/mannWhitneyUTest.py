#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
import pandas as pd
import genericLib as gL
import itertools as itt
from scipy.stats import mannwhitneyu

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
nSamples = int(sys.argv[1])
nBatches = int(sys.argv[2])

modelId = 'ENGRO2'
lcellLines =  ['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']
allCombs = list(itt.combinations_with_replacement(lcellLines, 2))
allCombs=[(test[0],test[1]) for test in allCombs if test[0]!=test[1]]

## load input random sampling datasets
dOFdf = {}

for cellLine in lcellLines:
    df=pd.DataFrame()
    for i in range(nBatches):
        df2 =pd.read_csv(os.path.join("outputs\\randomSampling_ENGRO2_nSol_"+ str(nSamples)+"_nBatch_"+str(i) + '_'+cellLine+'.csv'), sep = '\t')
        df=df.append(df2.round(10))
    df=df.reset_index()
    
    dOFdf[cellLine] = df
    
## compute mann-whitney U test
valore_nan=np.nan

for comb in allCombs:
    print(comb)
    df1=dOFdf[comb[0]]
    df2=dOFdf[comb[1]]   

    lCols = dOFdf[comb[0]].columns.tolist()[1:]

    df_tests=pd.DataFrame(index=lCols,columns=[ 'statistic_less', 'pvalue_less', 'statistic_greater', 'pvalue_greater', 'mean_' + comb[0], 'median_' + comb[0], 'std_' + comb[0], 'mean_' + comb[1],'median_' + comb[1], 'std_' + comb[1]])
    
    lines=list()
    for col in lCols:
        df1_col=df1[col]
        df2_col=df2[col]
    
        if ((df1_col == 0).all() == True) and ((df2_col == 0).all() == True):

            lines.append([ valore_nan, valore_nan, valore_nan, valore_nan, df1_col.mean(), df1_col.median(), df1_col.std(),df2_col.mean(), df2_col.median(), df2_col.std()])
        elif df1_col.equals(df2_col) == True:

            lines.append([ valore_nan, valore_nan, valore_nan, valore_nan, df1_col.mean(), df1_col.median(), df1_col.std(),df2_col.mean(), df2_col.median(), df2_col.std()])
        else:
            resultsLess = mannwhitneyu(df1_col, df2_col, alternative = 'less')
            resultsGreater = mannwhitneyu(df1_col, df2_col, alternative = 'greater')
            lines.append([ resultsLess.statistic, resultsLess.pvalue, resultsGreater.statistic, resultsGreater.pvalue, df1_col.mean(), df1_col.median(), df1_col.std(),df2_col.mean(), df2_col.median(), df2_col.std()])
    
    
    df_tests.iloc[:,:]=lines

    df_tests.to_csv(os.path.join(OUTDIR, 'mwuTest_' + comb[0] +'_vs_' + comb[1]  +'.csv'),sep="\t")
