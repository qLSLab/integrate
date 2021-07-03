#!/usr/bin/env python
# coding: utf-8
import cobra as cb
import genericLib as gL
import os
import sys

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
modelId = 'ENGRO2'
modelName = 'ENGRO2_reversible_20210305'

# load model
model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName + '.xml'))

# generate output file
timeStamp = gL.getTimeStamp()
outFile = gL.openFile(modelId + '_gpr_' + timeStamp, 'csv', dir=OUTDIR)
gL.writeLineByLineToFile(outFile, ['id', 'rule'], '\t')

for r in model.reactions:
    gL.writeLineByLineToFile(outFile, [r.id, r.gene_reaction_rule], '\t')

outFile.close()
