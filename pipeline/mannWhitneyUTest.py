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
modelId = 'ENGRO2'
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']

## load input random sampling datasets
dOFdf = {}
for cellLine in lcellLines:
    df = pd.read_csv(os.path.join(OUTDIR, 'randomSampling_' + modelId + '_nSol_' + str(nSamples) + '_' + cellLine + '.csv'), sep = '\t')
    dOFdf[cellLine] = df

## compute mann-whitney U test
allCombs = list(itt.combinations(lcellLines, 2))

for comb in allCombs:
    lCols = dOFdf[comb[0]].columns.tolist()
    outFile = open(os.path.join(OUTDIR, 'mwuTest_' + comb[0] +'_vs_' + comb[1] + '.csv'), mode='w')
    gL.writeLineByLineToFile(outFile, ['Rxn', 'statistic_less', 'pvalue_less', 'statistic_greater', 'pvalue_greater', 'mean_' + comb[0], 'median_' + comb[0], 'std_' + comb[0], 'mean_' + comb[1],'median_' + comb[1], 'std_' + comb[1]], '\t')
    for col in lCols:
        if ((dOFdf[comb[0]][col].round(5) == 0).all() == True) and ((dOFdf[comb[1]][col].round(5) == 0).all() == True):
            gL.writeLineByLineToFile(outFile, [col, np.nan, np.nan, np.nan, np.nan, dOFdf[comb[0]][col].round(5).mean(), dOFdf[comb[0]][col].round(5).median(), dOFdf[comb[0]][col].round(5).std(),dOFdf[comb[1]][col].round(5).mean(), dOFdf[comb[1]][col].round(5).median(), dOFdf[comb[1]][col].round(5).std()], '\t')
        elif dOFdf[comb[0]][col].round(5).equals(dOFdf[comb[1]][col].round(5)) == True:
            gL.writeLineByLineToFile(outFile, [col, np.nan, np.nan, np.nan, np.nan, dOFdf[comb[0]][col].round(5).mean(), dOFdf[comb[0]][col].round(5).median(), dOFdf[comb[0]][col].round(5).std(),dOFdf[comb[1]][col].round(5).mean(), dOFdf[comb[1]][col].round(5).median(), dOFdf[comb[1]][col].round(5).std()], '\t')
        else:
            resultsLess = mannwhitneyu(dOFdf[comb[0]][col].round(5), dOFdf[comb[1]][col].round(5), alternative = 'less')
            resultsGreater = mannwhitneyu(dOFdf[comb[0]][col].round(5), dOFdf[comb[1]][col].round(5), alternative = 'greater')
            gL.writeLineByLineToFile(outFile, [col, resultsLess.statistic, resultsLess.pvalue, resultsGreater.statistic, resultsGreater.pvalue, dOFdf[comb[0]][col].round(5).mean(), dOFdf[comb[0]][col].round(5).median(), dOFdf[comb[0]][col].round(5).std(),dOFdf[comb[1]][col].round(5).mean(), dOFdf[comb[1]][col].round(5).median(), dOFdf[comb[1]][col].round(5).std()], '\t')
    outFile.close()
