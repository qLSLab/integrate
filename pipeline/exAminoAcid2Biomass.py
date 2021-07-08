#!/usr/bin/env python
# coding: utf-8

import cobra as cb
import pandas as pd
import genericLib as gL
import os
import sys
import dataVizLib as dvl

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]

# setting input data
tStampRS = sys.argv[1]
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']
modelName = 'ENGRO2_MCF102A_wIrrRxns.xml'
dizMet2Rxn = {
    'leu__L_c': 'EX_leu__L_e_r',
    'ile__L_c': 'EX_ile__L_e_r',
    'tyr__L_c': 'EX_tyr__L_e_r',
    'val__L_c': 'EX_val__L_e_r'
}
biomassRxn = 'Biomass'
modelId = 'ENGRO2'
nSamples = 1000000

## Load random sampling input files
lOFdf = []
for c in lcellLines:
    inputName = 'randomSampling_' + modelId + '_nSol_' + str(nSamples) + '_' + cellLine + '_' + tStampRS
    df = pd.read_csv(os.path.join(OUTDIR, inputName + '.csv'), sep = '\t')
    df['EX_leu__L_e_r'] = -1 * df['EX_leu__L_e_r']
    df['EX_ile__L_e_r'] = -1 * df['EX_ile__L_e_r']
    df['EX_tyr__L_e_r'] = -1 * df['EX_tyr__L_e_r']
    df['EX_val__L_e_r'] = -1 * df['EX_val__L_e_r']
    lOFdf.append(df)

## Load input model
model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName))

bioMassRx = str(model.reactions.get_by_id(biomassRxn).reaction)
BMReags = bioMassRx.split('-->')[0].split('+')

## Compute contribution of amino acids to biomass synthesis reaction
dBioMassMetsExRxs = {}
for term in BMReags:
    metAndCoef = term.strip()
    coeff, met = metAndCoef.split()
    print(met, float(coeff))
    if met in dizMet2Rxn:
        dBioMassMetsExRxs[met[:-2]] = {
            'coef': float(coeff),
            'ExRxn': dizMet2Rxn[met],
            'BioM': met[:-2] + '_BioM'
        }

dfAll = pd.concat(lOFdf)
dfAll['sampleID'] = ['MCF102A_' + str(i) for i in dfMCF102ANet.index.tolist()] + ['SKBR3_' + str(i) for i in dfMCF102ANet.index.tolist()] + ['MCF7_' + str(i) for i in dfMCF102ANet.index.tolist()] + ['MDAMB231_' + str(i) for i in dfMCF102ANet.index.tolist()] + ['MDAMB361_' + str(i) for i in dfMCF102ANet.index.tolist()]
dfAll['cellLine'] = dfAll.sampleID.str.split('_', expand=True)[0]
fluxes = dfAll.columns[:-2]

columns = ['cellLine', 'sampleID']
columns += fluxes.to_list()
dfAll = dfAll[columns]

lBioMExRxns = list(dizMet2Rxn.values())
dfExFlux = dfAll[['cellLine', 'Biomass'] + lBioMExRxns]

for flux in dBioMassMetsExRxs:
    fluxInBioM = flux + '_BioM'
    dfExFlux[fluxInBioM] = dfExFlux.Biomass * dBioMassMetsExRxs[flux]['coef']

outFileName = os.path.join(OUTDIR, 'rs_' + tStampRS + '_BioMvsEX.csv')
dfExFlux.to_csv(outFileName)

for flux in dBioMassMetsExRxs:
    outFileName = os.path.join(FIGUREDIR,
                               'rs_' + tStampRS + '_BioMvsEX' + flux + '.png')
    dvl.vExcngInBiomass(dfExFlux,
                        dBioMassMetsExRxs[flux]['BioM'],
                        dBioMassMetsExRxs[flux]['ExRxn'],
                        dBioMassMetsExRxs[flux]['coef'],
                        lcellLines,
                        sTitle=flux,
                        file2save=outFileName)
