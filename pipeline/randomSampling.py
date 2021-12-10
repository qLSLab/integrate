#!/usr/bin/env python
# coding: utf-8

import os
import sys
import cobra as cb
import genericLib as gL
import pandas as pd
import numpy as np
from cobra.sampling import OptGPSampler

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
nSamples = int(sys.argv[1])
nBatches = int(sys.argv[2])


biomassRxn = 'Biomass'
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']
modelId = 'ENGRO2'
ysiFileName = 'ysi_ratio.csv'

lacRxn = 'EX_lac__L_e'
glcRxn = 'EX_glc__D_e_r'
glnRxn = 'EX_gln__L_e_r'
gluRxn = 'EX_gluOUT__L_e'
lReplicas = ['_A','_B']


#Yield Information
valore=1e-3*180.15
maxYield=0.000167998*valore
minYield=3.90762E-05*valore


tStampRS = gL.getTimeStamp()

if __name__ == '__main__':
    dfYSI = pd.read_csv(os.path.join(RAWDIR, ysiFileName), sep = '\t', index_col = 'Ratio')
    dYSI = dfYSI.to_dict()

    for cellLine in lcellLines:
        print(cellLine)
        modelName = modelId + '_' + cellLine + '_wIrrRxns'
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName + '.xml'))
            
        model.reactions.get_by_id(biomassRxn).objective_coefficient = 0

        #Yield
        constrainUp = model.problem.Constraint(0.131972*model.reactions.Biomass.flux_expression-
                                              maxYield*model.reactions.EX_glc__D_e_r.flux_expression, 
                                                     lb=-1000, 
                                                     ub=0)
        constrainUp.name="yield_up" 

        constrainDown = model.problem.Constraint(0.131972*model.reactions.Biomass.flux_expression-
                                              minYield*model.reactions.EX_glc__D_e_r.flux_expression, 
                                                     lb=0, 
                                                     ub=1000)    
        constrainDown.name="yield_down" 

        model.add_cons_vars(constrainUp)
        model.add_cons_vars(constrainDown)      
        
        
        lacL = model.reactions.get_by_id(lacRxn)
        glc = model.reactions.get_by_id(glcRxn)
        gln = model.reactions.get_by_id(glnRxn)
        glu = model.reactions.get_by_id(gluRxn)

        lgluGln = []
        for replica in lReplicas:
            lgluGln.append(dYSI[cellLine + replica]['glu/gln'])
        gluGln_mean = np.mean(lgluGln)
        gluGln_std = np.std(lgluGln)
        gluGln = model.problem.Constraint(glu.flux_expression - gluGln_mean * gln.flux_expression, lb=-gluGln_std, ub=gluGln_std)
        model.add_cons_vars(gluGln)

        llacGlc = []
        for replica in lReplicas:
            llacGlc.append(dYSI[cellLine + replica]['lac/glc'])
        lacGlc_mean = np.mean(llacGlc)
        lacGlc_std = np.std(llacGlc)
        lacGlc = model.problem.Constraint(lacL.flux_expression - lacGlc_mean * glc.flux_expression, lb=-lacGlc_std, ub=lacGlc_std)
        model.add_cons_vars(lacGlc)

        llacGln = []
        for replica in lReplicas:
            llacGln.append(dYSI[cellLine + replica]['lac/gln'])
        lacGln_mean = np.mean(llacGln)
        lacGln_std = np.std(llacGln)
        lacGln = model.problem.Constraint(lacL.flux_expression - lacGln_mean * gln.flux_expression, lb=-lacGln_std, ub=lacGln_std)
        model.add_cons_vars(lacGln)

        #samples with batches
        optgp = OptGPSampler(model, processes=2,thinning=10,seed=20)
        k=0
        for s in optgp.batch(nSamples, nBatches):
            outputName = 'randomSampling_' + modelId + '_nSol_' + str(nSamples) +"_nBatch_"+ str(k) + '_' + cellLine
            s.to_csv(os.path.join(OUTDIR, outputName+ '.csv'), sep = '\t', index=False)
            k=k+1


