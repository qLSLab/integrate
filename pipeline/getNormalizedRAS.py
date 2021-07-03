#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import genericLib as gL

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
inputFileName = 'ENGRO2_RAS'
outputFileName = 'ENGRO2_wNormalizedRAS'

# Load RAS dataset
RAS = pd.read_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index_col = 'Rxn')
## Remaining NaN in the matrix will be substituted with 1 to indicate that a missing RAS value is originally present for these reactions due to emtpy GPR rules
RAS = RAS.fillna(1)

# Normalization for each reaction on the sample having the highest RAS. Here samples from the transcriptomics dataset in rawData directory are used.
# Change sample names and number of replicates (nReplicas variable) of each sample accordingly if the input transcriptomic dataset is different.
nReplicas = 3.0
RAS['media_MDAMB361']=(RAS['Cellule361_Rep1']+RAS['Cellule361_Rep2']+RAS['Cellule361_Rep3'])/nReplicas
RAS['media_MDAMB231']=(RAS['MDA-MB-231_Rep1']+RAS['MDA-MB-231_Rep2']+RAS['MDA-MB-231_Rep3'])/nReplicas
RAS['media_MCF7']=(RAS['MCF-7_Rep1']+RAS['MCF-7_Rep2']+RAS['MCF-7_Rep3'])/nReplicas
RAS['media_SKBR3']=(RAS['SK-BR-3_Rep1']+RAS['SK-BR-3_Rep2']+RAS['SK-BR-3_Rep3'])/nReplicas
RAS['media_MCF102A']=(RAS['MCF102A_Rep1']+RAS['MCF102A_Rep2']+RAS['MCF102A_Rep3'])/nReplicas

RASNormalized = RAS[['media_MDAMB361','media_MDAMB231','media_MCF7','media_SKBR3','media_MCF102A']]
maxValues = RASNormalized.max(axis=1)

RASNormalized['norm_MDAMB361']=RASNormalized['media_MDAMB361']/maxValues
RASNormalized['norm_MDAMB231']=RASNormalized['media_MDAMB231']/maxValues
RASNormalized['norm_MCF7']=RASNormalized['media_MCF7']/maxValues
RASNormalized['norm_SKBR3']=RASNormalized['media_SKBR3']/maxValues
RASNormalized['norm_MCF102A']=RASNormalized['media_MCF102A']/maxValues

RASNormalized = RASNormalized.fillna(1)

## Mask rows where all values were originally equal to 0, by putting 0: if all values in rows are equal to 0, then the maximum value in this row is 0 and 0 / 0 is
# undetermineted which means NaN that before was substituted with 1. However, it's important to mantain the 0 value to indicate that the RAS is not missing but it's 0 and
# consequently the reactions is off.
RASNormalized.loc[(RAS.eq(0).all(1))] = 0

# save output matrix
RASNormalized.to_csv(os.path.join(OUTDIR,outputFileName + '.csv'), sep = '\t', index = True)
