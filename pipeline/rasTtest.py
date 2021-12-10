#!/usr/bin/env python
# coding: utf-8

import os
import cobra as cb
import pandas as pd
import genericLib as gL
from scipy.stats import ttest_ind
import itertools as itt
import numpy as np

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
lcellLines = ['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']
rasScoreFile = 'ENGRO2_RAS'

## According to the splitted reactions in the irreversible models, the corresponding RAS scores are duplicated to consider both the forward and the backward reaction
RAS = pd.read_csv(os.path.join(OUTDIR, rasScoreFile + '.csv'), sep="\t", index_col = 'Rxn')
RAS = RAS.fillna(1)

dRxn2SpltRxn = {}
model = cb.io.sbml.read_sbml_model(os.path.join(MODELDIR, 'ENGRO2_MCF102A_wIrrRxns.xml'))

for rxn in model.reactions:
    if rxn.id.endswith('_f') is True or rxn.id.endswith('_b') is True or rxn.id.endswith('_r') is True:
        if rxn.id[:-2] not in dRxn2SpltRxn:
            dRxn2SpltRxn[rxn.id[:-2]] = [rxn.id]
        else:
            dRxn2SpltRxn[rxn.id[:-2]] += [rxn.id]
    else:
        if rxn.id not in dRxn2SpltRxn:
            dRxn2SpltRxn[rxn.id] = [rxn.id]
        else:
            dRxn2SpltRxn[rxn.id] += [rxn.id]

dfRxn2SpltRxns = pd.DataFrame(dRxn2SpltRxn.items(), columns=['Rxn', 'RxnSplt'])

dfmerge = pd.merge(RAS, dfRxn2SpltRxns, on='Rxn')
dfmerge = dfmerge.explode('RxnSplt')
dfmerge = dfmerge.drop(columns=['Rxn'])
dfmerge = dfmerge.rename(columns={"RxnSplt": "Rxn"})
dfmerge = dfmerge.set_index('Rxn')

dfmerge.to_csv(os.path.join(OUTDIR, rasScoreFile + '_splt.csv'), sep = '\t', index = True)

## Load computed RAS scores
dfmerge = dfmerge.dropna(how='all')
dfmerge = dfmerge.rename(columns={'Cellule361_Rep1': 'MDAMB361_Rep1', 'Cellule361_Rep2': 'MDAMB361_Rep2', 'Cellule361_Rep3': 'MDAMB361_Rep3', 'MDA-MB-231_Rep1': 'MDAMB231_Rep1', 'MDA-MB-231_Rep2': 'MDAMB231_Rep2', 'MDA-MB-231_Rep3': 'MDAMB231_Rep3',
                    'MCF-7_Rep1': 'MCF7_Rep1', 'MCF-7_Rep2': 'MCF7_Rep2', 'MCF-7_Rep3': 'MCF7_Rep3', 'SK-BR-3_Rep1': 'SKBR3_Rep1', 'SK-BR-3_Rep2': 'SKBR3_Rep2', 'SK-BR-3_Rep3': 'SKBR3_Rep3'})
dfmerge = dfmerge.round(5)

## Compute t-test for all the pairs of cell lines
allCombs = list(itt.combinations(lcellLines, 2))
for comb in allCombs:
    outFile = open(os.path.join(OUTDIR, 'ras_' + comb[0] +'_vs_' + comb[1] + '.csv'), mode='w')
    gL.writeLineByLineToFile(outFile, ['Rxn', 'statistic_less', 'pvalue_less', 'statistic_greater', 'pvalue_greater', 'mean_' + comb[0], 'std_' + comb[0], 'mean_' + comb[1], 'std_' + comb[1]], '\t')

    for idx, row in dfmerge.iterrows():
        lras_comb0 = [getattr(row, comb[0] + '_Rep1'), getattr(row, comb[0] + '_Rep2'), getattr(row, comb[0] + '_Rep3')]
        lras_comb1 = [getattr(row, comb[1] + '_Rep1'), getattr(row, comb[1] + '_Rep2'), getattr(row, comb[1] + '_Rep3')]
        testT_l = ttest_ind(lras_comb0, lras_comb1, alternative = 'less')
        testT_g = ttest_ind(lras_comb0, lras_comb1, alternative = 'greater')
        gL.writeLineByLineToFile(outFile, [idx, testT_l.statistic, testT_l.pvalue, testT_g.statistic, testT_g.pvalue, np.mean(lras_comb0), np.std(lras_comb0), np.mean(lras_comb1), np.std(lras_comb1)], '\t')
    outFile.close()


## Analyse output and compute log2 of fold change
allCombs = list(itt.combinations(lcellLines, 2))

for comb in allCombs:
    inputFile = 'ras_' + comb[0] +'_vs_' + comb[1] + '.csv'
    df = pd.read_csv(os.path.join(OUTDIR, inputFile), sep = '\t')
    df['result'] = np.nan

    filter = (df['pvalue_less'] <= 0.05) & (df['pvalue_greater'] > 0.05)
    df['result'].loc[filter] = -1

    filter = (df['pvalue_less'] > 0.05) & (df['pvalue_greater'] <= 0.05)
    df['result'].loc[filter] = 1

    filter = (df['pvalue_less'] <= 0.05) & (df['pvalue_greater'] <= 0.05)
    df['result'].loc[filter] = 2

    filter = (df['pvalue_less'] > 0.05) & (df['pvalue_greater'] > 0.05)
    df['result'].loc[filter] = 0

    filter = (np.isnan(df['pvalue_less'])) & (np.isnan(df['pvalue_greater']))
    df['result'].loc[filter] = 0

    df.to_csv(os.path.join(OUTDIR, inputFile), sep = '\t', index = False)

allCombs = list(itt.combinations(lcellLines, 2))

for comb in allCombs:
    inputFile = 'ras_' + comb[0] +'_vs_' + comb[1] + '.csv'
    df = pd.read_csv(os.path.join(OUTDIR, inputFile), sep = '\t')
    df['meanRatio'] = df['mean_' + comb[0]] / df['mean_' + comb[1]]
    df['meanRatio'] = df['meanRatio'].abs()
    df['meanRatio_log'] = np.log2(df['meanRatio'])
    filter = (df['meanRatio'] == 0)
    df['meanRatio_log'].loc[filter] = np.nan
    df.to_csv(os.path.join(OUTDIR, inputFile), sep = '\t', index = False)
