#!/usr/bin/env python
# coding: utf-8

import os
import cobra as cb
import pandas as pd
import genericLib as gL
import numpy as np
import sys

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
imposeYSI = 'Y' # (Y/N)
imposeMedium = 'Y' # (Y/N)
imposeRasConstraints = 'Y' # (Y/N)
rasNormFileName = 'ENGRO2_wNormalizedRAS'
ysiFileName = 'ysi_ratio.csv'
mediumFileName = 'medium.csv'
modelId = 'ENGRO2'
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']
biomassRxn = 'Biomass'
lacRxn = 'EX_lac__L_e'
glcRxn = 'EX_glc__D_e'
glnRxn = 'EX_gln__L_e'
gluRxn = 'EX_gluOUT__L_e'
lReplicas = ['_A','_B']

#Yield Information
valore=1e-3*180.15
maxYield=0.000167998*valore
minYield=3.90762E-05*valore

if __name__ == '__main__':
    ## load YSI dataset
    dfYSI = pd.read_csv(os.path.join(RAWDIR, ysiFileName), sep = '\t', index_col = 'Ratio')
    dYSI = dfYSI.to_dict()

    ## load RAS dataset
    RAS_norm = pd.read_csv(os.path.join(OUTDIR, rasNormFileName + '.csv'), sep = '\t', index_col = 'Rxn')
    lIdx_ras = RAS_norm.index.tolist()

    ## Load medium dataset
    medium = pd.read_csv(os.path.join(RAWDIR, mediumFileName), sep = '\t')

    ## Apply constraints
    for cellLine in lcellLines:
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelId + '.xml'))
        model.reactions.get_by_id(biomassRxn).objective_coefficient = 0
        #model.solver = 'gurobi'

        #Yield
        constrainUp = model.problem.Constraint(0.131972*model.reactions.Biomass.flux_expression+
                                              maxYield*model.reactions.EX_glc__D_e.flux_expression, 
                                                     lb=-1000, 
                                                     ub=0)
        constrainUp.name="yield_up" 

        constrainDown = model.problem.Constraint(0.131972*model.reactions.Biomass.flux_expression+
                                              minYield*model.reactions.EX_glc__D_e.flux_expression, 
                                                     lb=0, 
                                                     ub=1000)    
        constrainDown.name="yield_down" 

        model.add_cons_vars(constrainUp)
        model.add_cons_vars(constrainDown)        
            
        if imposeYSI == 'Y':
            lacL = model.reactions.get_by_id(lacRxn)
            glc = model.reactions.get_by_id(glcRxn)
            gln = model.reactions.get_by_id(glnRxn)
            glu = model.reactions.get_by_id(gluRxn)

            lgluGln = []
            for replica in lReplicas:
                lgluGln.append(dYSI[cellLine + replica]['glu/gln'])
            gluGln_mean = np.mean(lgluGln)
            gluGln_std = np.std(lgluGln)
            gluGln = model.problem.Constraint(glu.flux_expression - gluGln_mean * (gln.flux_expression * -1), lb=-gluGln_std, ub=gluGln_std)
            model.add_cons_vars(gluGln)

            llacGlc = []
            for replica in lReplicas:
                llacGlc.append(dYSI[cellLine + replica]['lac/glc'])
            lacGlc_mean = np.mean(llacGlc)
            lacGlc_std = np.std(llacGlc)
            lacGlc = model.problem.Constraint(lacL.flux_expression - lacGlc_mean * (-1 * glc.flux_expression), lb=-lacGlc_std, ub=lacGlc_std)
            model.add_cons_vars(lacGlc)

            llacGln = []
            for replica in lReplicas:
                llacGln.append(dYSI[cellLine + replica]['lac/gln'])
            lacGln_mean = np.mean(llacGln)
            lacGln_std = np.std(llacGln)
            lacGln = model.problem.Constraint(lacL.flux_expression - lacGln_mean * (-1 * gln.flux_expression), lb=-lacGln_std, ub=lacGln_std)
            model.add_cons_vars(lacGln)

        if imposeMedium == 'Y':
            for row in medium.itertuples():
                id = row.Rxn
                lb = getattr(row, cellLine + '_LB')
                ub = getattr(row, cellLine + '_UB')
                model.reactions.get_by_id(id).lower_bound = lb
                model.reactions.get_by_id(id).upper_bound = ub

        if imposeRasConstraints == 'Y':
            FVA_no = cb.flux_analysis.flux_variability_analysis(model).round(10)
            FVA_no.to_csv("FVA_"+cellLine+".csv")
            dfRASFVA = pd.merge(RAS_norm, FVA_no, how='inner', left_on='Rxn', right_index=True)
            dfRASFVA['minimum_abs'] = abs(dfRASFVA.minimum)
            dfRASFVA['maximum_abs'] = abs(dfRASFVA.maximum)
            modelRxns = []
            modelLb = []
            modelUb = []
            for r in model.reactions:
                modelRxns.append(r.id)
                modelLb.append(r.lower_bound)
                modelUb.append(r.upper_bound)
            dfModelConstraints = pd.DataFrame({'Rxn': modelRxns, 'Lb': modelLb, 'Ub': modelUb})
            dfRASFVAModelConstraints = pd.merge(dfRASFVA, dfModelConstraints, how='inner', on='Rxn')
            dfRASFVAModelConstraints['lb_sign'] = np.sign(dfRASFVAModelConstraints.Lb)
            dfRASFVAModelConstraints['ub_sign'] = np.sign(dfRASFVAModelConstraints.Ub)


            epsilon = 0
            dfRASFVAModelConstraints[cellLine + 'Low'] = (epsilon + ((dfRASFVAModelConstraints.minimum_abs - epsilon) * getattr(dfRASFVAModelConstraints, 'norm_' + cellLine))) * getattr(dfRASFVAModelConstraints, 'lb_sign')
            dfRASFVAModelConstraints[cellLine + 'Up'] = (epsilon + ((dfRASFVAModelConstraints.maximum_abs - epsilon) * getattr(dfRASFVAModelConstraints, 'norm_' + cellLine))) * getattr(dfRASFVAModelConstraints, 'ub_sign')

            maskMin = (dfRASFVAModelConstraints['minimum'] == 0) & (dfRASFVAModelConstraints[cellLine + 'Low'] < 0)
            dfRASFVAModelConstraints.loc[maskMin, cellLine + 'Low'] = 0
            maskMax = (dfRASFVAModelConstraints['maximum'] == 0) & (dfRASFVAModelConstraints[cellLine + 'Up'] > 0)
            dfRASFVAModelConstraints.loc[maskMax, cellLine + 'Up'] = 0

            for row in dfRASFVAModelConstraints.itertuples():
                ID = row.Rxn
                if ID == 'CARPEPT1tc' or ID == 'HIStiDF':
                    model.reactions.get_by_id(ID).lower_bound = float(getattr(row, 'minimum'))
                    model.reactions.get_by_id(ID).upper_bound = float(getattr(row, 'maximum'))
                else:
                    model.reactions.get_by_id(ID).lower_bound = float(getattr(row, cellLine + 'Low'))
                    model.reactions.get_by_id(ID).upper_bound = float(getattr(row, cellLine + 'Up'))

        model.reactions.get_by_id(biomassRxn).objective_coefficient = 1.0

        cb.io.write_sbml_model(model, os.path.join(MODELDIR, modelId + '_' + cellLine + '.xml'))
