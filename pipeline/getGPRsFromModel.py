#!/usr/bin/env python
# coding: utf-8
import cobra as cb
import genericLib as gL
import os

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
modelId = 'ENGRO2'

# load model
model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelId + '.xml'))

# generate output file
outFile = gL.openFile(modelId + '_GPR', 'csv', dir=OUTDIR)
gL.writeLineByLineToFile(outFile, ['id', 'rule'], '\t')

for r in model.reactions:
    gL.writeLineByLineToFile(outFile, [r.id, r.gene_reaction_rule], '\t')

outFile.close()
