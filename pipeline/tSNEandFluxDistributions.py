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
FIGUREDIR = workingDirs[4]

# setting input data
tStampRS = sys.argv[1]
nExtraction = 0
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']

## load net fluxes from random sampling
lOFdf = []
for c in lcellLines:
    inputName = 'randomSampling_' + modelId + '_nSol_' + str(nSamples) + '_' + cellLine + '_' + tStampRS + '.csv'
    df = pd.read_csv(os.path.join(OUTDIR, inputName), sep = '\t')
    dfRS = iL.splt2Net(df)
    lIndexes = df.index.tolist()
    lNewIndexes = []
    for idx in lIndexes:
        lNewIndexes.append(c + '_' + str(idx))
    df.index = lNewIndexes
    lOFdf.append(df)

allDf_OFdf = pd.concat(lOFdf)
# allDf_OFdf = allDf_OFdf[['ACONT', 'OCOAT1m']] # uncomment if you want to just show the 'ACONT' and 'OCOAT1m' reactions

inRSFile = 'sampledFluxes_rs_' + tStampRS + '_extr_' + str(nExtraction) + '.csv.bz2'
basenameRS = inRSFile[:-8]
dfAll = pd.read_csv(os.path.join(OUTDIR, inRSFile), compression = 'bz2', sep='\t', index_col= 'Unnamed: 0')
# dfAll = dfAll[['ACONT', 'OCOAT1m']]  # uncomment if you want to just show the 'ACONT' and 'OCOAT1m' reactions

inTsneFile = 'tsnewPca_rs_' + tStampRS + '_extr_' + nExtraction + '.csv'
basenameTsne = inTsneFile[:-4]
tagFile = os.path.join(OUTDIR, basenameTsne + '.tags')

fluxes = dfAll.columns
# fluxes = ['ACONT', 'OCOAT1m'] # uncomment if you want to just show the 'ACONT' and 'OCOAT1m' reactions

dfTsne = pd.read_csv(os.path.join(OUTDIR, inTsneFile), sep='\t', index_col=0)
dfTsne.index = dfAll.index.tolist()

dfFlux = pd.merge(dfTsne, dfAll, left_index=True, right_index=True)

allDf_OFdf['cellLine'] = [i.split('_')[0] for i in allDf_OFdf.index.tolist()]

for flux in fluxes:
    (dfToPlot, nBins, bins) = dvl.restrict_q5_q95(allDf_OFdf, flux)
    (dfToPlot_dfFlux, nBins_dfFlux, bins_dfFlux) = dvl.restrict_q5_q95(dfFlux, flux)
    outFileName = os.path.join(
        FIGUREDIR, basenameRS + '_ScatAndHist' + flux + '.pdf')
    dvl.scatAndHistOnTsne(dfToPlot_dfFlux, dfToPlot, flux, sTitle=flux,
                          tagFileName=tagFile, tagColumn='cellLine',
                          binDistrib=bins, file2save=outFileName)
