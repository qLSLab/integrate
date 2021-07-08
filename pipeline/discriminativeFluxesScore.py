#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import genericLib as gL
import dataVizLib as dvl
import integrateLib as iL

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
tStampRS = sys.argv[1]
modelId = 'ENGRO2'
nSamples = 1000000
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']

lAllDf = []
for cellLine in lcellLines:
    inputName = 'randomSampling_' + modelId + '_nSol_' + str(nSamples) + '_' + cellLine + '_' + tStampRS + '.csv'
    dfRS = pd.read_csv(os.path.join(OUTDIR, inputName), sep='\t')
    dfRS = iL.splt2Net(dfRS)
    lIndexes = dfRS.index.tolist()
    lNewIndexes = []
    for idx in lIndexes:
        lNewIndexes.append(cellLine + '_' + str(idx))
    dfRS['sampleID'] = lNewIndexes
    dfRS['cellLine'] = dfRS.sampleID.str.split('_', expand=True)[0]
    lAllDf.append(dfRS)

dfAll = pd.concat(lAllDf)

lFluxes = dfAll.columns.to_list()
lFluxes.remove('cellLine')
lFluxes.remove('sampleID')

dfFlux = dfAll[['cellLine', lFluxes[0]]]
(dfFluxReduced, nBins, bins) = dvl.restrict_q5_q95(dfFlux, lFluxes[0])
dfFLuxesScores = dvl.separationScoreDF(
    dfFluxReduced, lFluxes[0], lcellLines, bins)

lFluxes.remove(lFluxes[0])

for flux in lFluxes:
    dfFlux = dfAll[['cellLine', flux]]
    (dfFluxReduced, nBins, bins) = dvl.restrict_q5_q95(dfFlux, flux)
    dfFLuxesScores = pd.merge(dfFLuxesScores, dvl.separationScoreDF(
        dfFluxReduced, flux, lcellLines, bins),
        left_index=True, right_index=True)

dfFLuxesScores = dfFLuxesScores.T
dfFLuxesScores.sort_values(by='totalScore', inplace=True, ascending=False)

basenameRS = 'rs_' + tStampRS
outFileName = os.path.join(OUTDIR, basenameRS + '_fluxScores.csv')
dfFLuxesScores.to_csv(outFileName, sep='\t')
