#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import genericLib as gL
import integrateLib as iL
import cobra as cb
import numpy as np
import re

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
gprRule = 'ENGRO2_GPR'
rnaSeqFileName = 'FPKM_Breast_forMarea.tsv'
modelId = 'ENGRO2'
regexOrgSpecific = r"([A-Z0-9.]+)"

# load model and extract genes list
model = cb.io.sbml.read_sbml_model(os.path.join(MODELDIR, modelId + ".xml"))
lModelGenes = []
for g in model.genes:
    lModelGenes.append(g.id)

# load RNA-seq data
rnaSeq = pd.read_csv(os.path.join(RAWDIR, rnaSeqFileName), sep = '\t', dtype = {'EnsID': str})
lCols = rnaSeq.columns.tolist()
lSamples = [col for col in lCols if col != 'EnsID']

# load input file with reactions and their GPR rules
dfRxnRule = pd.read_csv(os.path.join(OUTDIR, gprRule +'.csv'), sep = '\t', dtype = {'rule': str})

# Compute RAS score and reconstruct output file
outFile = gL.openFile(modelId + '_RAS', 'csv', dir=OUTDIR)
gL.writeLineByLineToFile(outFile, ['Rxn'] + lSamples, '\t')

for row in dfRxnRule.itertuples():
    if row.rule == '' or pd.isna(row.rule) is True:
        gL.writeLineByLineToFile(outFile, [row.id] + [float("NaN")] * len(lSamples), '\t')
    elif (row.rule != '' or pd.isna(row.rule) is False) and ' or ' not in row.rule and ' and ' not in row.rule: # if rule includes one gene
        dfSearchTranscript = rnaSeq.loc[rnaSeq['EnsID'] == row.rule]
        if dfSearchTranscript.empty is True:
            gL.writeLineByLineToFile(outFile, [row.id] + [float("NaN")] * len(lSamples), '\t')
        else:
            lSample_transcripts = []
            el2Write = [row.id]
            for s in lSamples:
                el2Write.append(dfSearchTranscript.iloc[0][s])
            gL.writeLineByLineToFile(outFile, el2Write, '\t')
    elif (row.rule != '' or pd.isna(row.rule) is False) and ' or ' in row.rule and ' and ' not in row.rule: # if GPR rule includes only OR operators
        dfgenes = gL.extractRegexFromItem(row.rule, regexOrgSpecific)
        lGenes = list(dfgenes[0])
        intersezione = gL.intersect(lGenes, list(rnaSeq['EnsID']))
        if len(intersezione) < len(lGenes):
            dfSearchTranscript = rnaSeq[rnaSeq['EnsID'].isin(intersezione)]
            dfSearchTranscript = dfSearchTranscript.reset_index(drop=True)
        else:
            dfSearchTranscript = rnaSeq[rnaSeq['EnsID'].isin(lGenes)]
            dfSearchTranscript = dfSearchTranscript.reset_index(drop=True)
        el2Write = [row.id]
        for s in lSamples:
            el2Write.append(dfSearchTranscript[s].sum())
        gL.writeLineByLineToFile(outFile, el2Write, '\t')
    elif (row.rule != '' or pd.isna(row.rule) is False) and ' or ' not in row.rule and ' and ' in row.rule: # if GPR rule includes only AND operators
        dfgenes = gL.extractRegexFromItem(row.rule, regexOrgSpecific)
        lGenes = list(dfgenes[0])
        intersezione = gL.intersect(lGenes, list(rnaSeq['EnsID']))
        if len(intersezione) < len(lGenes):
            dfSearchTranscript = rnaSeq[rnaSeq['EnsID'].isin(intersezione)]
            dfSearchTranscript = dfSearchTranscript.reset_index(drop=True)
        else:
            dfSearchTranscript = rnaSeq[rnaSeq['EnsID'].isin(lGenes)]
            dfSearchTranscript = dfSearchTranscript.reset_index(drop=True)
        el2Write = [row.id]
        for s in lSamples:
            el2Write.append(dfSearchTranscript[s].min())
        gL.writeLineByLineToFile(outFile, el2Write, '\t')
    elif (row.rule != '' or pd.isna(row.rule) is False) and ' or ' in row.rule and ' and ' in row.rule: # if GPR rule includes both AND and OR operators
        el2Write = [row.id]
        for s in lSamples:
            regola = row.rule
            lopenPar = [i.start() for i in re.finditer('\(', regola)]
            lclosePar = [i.start() for i in re.finditer('\)', regola)]
            while len(lopenPar) != 0 and len(lclosePar) != 0:
                closePar = lclosePar[0]
                j = 0
                while j < len(lopenPar) and lopenPar[j] < closePar:
                    j += 1
                openPar = lopenPar[j-1]
                if closePar == len(regola) - 1:
                    fine = True
                else:
                    fine = False
                tmp = regola[openPar+1:closePar].strip()
                if '(' in tmp:
                    otherOpenPar = regola.find('(')
                    realOpen = otherOpenPar + openPar + 1
                    if fine is True:
                        subrule = regola[realOpen:]
                    else:
                        subrule = regola[realOpen:closePar + 1]
                else:
                    if fine is True:
                        subrule = regola[openPar:]
                    else:
                        subrule = regola[openPar:closePar + 1]
                dfgenes = gL.extractRegexFromItem(subrule, regexOrgSpecific)
                lGenes = list(dfgenes[0])
                inside = gL.intersect(lGenes, list(rnaSeq['EnsID']))
                diff = iL.differenceKeepingDuplicates(lGenes, list(rnaSeq['EnsID']))
                diffNoModelGenes = iL.differenceKeepingDuplicates(diff, lModelGenes)
                scoreAlreadyComputed = []
                for el in diffNoModelGenes:
                    try:
                        if np.isnan(float(el)) == False:
                            scoreAlreadyComputed.append(float(el))
                    except:
                        print(el, '\tis not a number')
                out = rnaSeq[rnaSeq['EnsID'].isin(inside)]
                out = out.reset_index(drop=True)
                score = iL.computeScore(subrule, out, s, scoreAlreadyComputed)
                regola = regola.replace(subrule, str(score))
                lopenPar = [i.start() for i in re.finditer('\(', regola)]
                lclosePar = [i.start() for i in re.finditer('\)', regola)]
            dfgenes = gL.extractRegexFromItem(regola, regexOrgSpecific)
            lGenes = list(dfgenes[0])
            inside = gL.intersect(lGenes, list(rnaSeq['EnsID']))
            diff = iL.differenceKeepingDuplicates(lGenes, list(rnaSeq['EnsID']))
            diffNoModelGenes = iL.differenceKeepingDuplicates(diff, lModelGenes)
            scoreAlreadyComputed = []
            for el in diffNoModelGenes:
                try:
                    if np.isnan(float(el)) == False:
                        scoreAlreadyComputed.append(float(el))
                except:
                    print(el, '\tis not a number')
            out = rnaSeq[rnaSeq['EnsID'].isin(inside)]
            out = out.reset_index(drop=True)
            score = iL.computeScore(regola, out, s, scoreAlreadyComputed)
            el2Write.append(score)
        gL.writeLineByLineToFile(outFile, el2Write, '\t')

outFile.close()
