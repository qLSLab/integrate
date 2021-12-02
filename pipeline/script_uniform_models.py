#!/usr/bin/env python
# coding: utf-8
import time
import cobra as cb
import pandas as pd
import genericLib as gL
import numpy as np
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import itertools as itt
import numpy as np
import itertools as itt
from cobra import Model, Reaction, Metabolite

def getRxns(model):
    lRxns = []
    for rxn in model.reactions:
        lRxns.append(rxn.id)
    return lRxns

def addbackwardForwardRxn(rId, lRxns, model):
    if rId not in lRxns:
        print('rId\t', rId)
        if rId.endswith('_f'):
            rGeneral = rId.split('_f')[0]
        elif rId.endswith('_b'):
            rGeneral = rId.split('_b')[0]
        print('rGeneral\t', rGeneral)
        if rGeneral in lRxns:
            model.reactions.get_by_id(rGeneral).id = model.reactions.get_by_id(rGeneral).id + '_f'
            model.repair()

            dReaction = model.reactions.get_by_id(model.reactions.get_by_id(rGeneral + '_f').id)._metabolites
            newdReaction = {}
            for k in dReaction:
                newdReaction[k] = -1 * dReaction[k]

            reaction = Reaction(rGeneral + '_b')
            reaction.name = model.reactions.get_by_id(rGeneral + '_f').name
            reaction.subsystem = model.reactions.get_by_id(rGeneral + '_f').subsystem
            reaction.lower_bound = 0
            reaction.upper_bound = 0
            reaction.add_metabolites(newdReaction)
            reaction.reaction
            model.add_reactions([reaction])
            # lRxns += [rGeneral + '_f', rGeneral + '_b']
            lRxns.remove(rGeneral)
            return model, [rGeneral, rGeneral + '_f', rGeneral + '_b'], lRxns
        else:
            return model, [], lRxns
    else:
        return model, [], lRxns

def addReverseRxn(rId, lRxns, model):
    if rId not in lRxns:
        rGeneral = rId.split('_r')[0]
        print('rGeneral\t', rGeneral)
        if rGeneral in lRxns:
            dReaction = model.reactions.get_by_id(rGeneral)._metabolites

            oldLb = model.reactions.get_by_id(rGeneral).lower_bound
            oldUb = model.reactions.get_by_id(rGeneral).upper_bound

            newdReaction = {}
            for k in dReaction:
                newdReaction[k] = -1 * dReaction[k]

            model.reactions.get_by_id(rGeneral).lower_bound = oldUb
            model.reactions.get_by_id(rGeneral).upper_bound = oldLb

            model.reactions.get_by_id(rGeneral).id = model.reactions.get_by_id(rGeneral).id + '_r'
            model.repair()

            model.reactions.get_by_id(rGeneral + '_r').subtract_metabolites(dReaction)
            model.reactions.get_by_id(rGeneral + '_r').add_metabolites(newdReaction)
            # lRxns += [rGeneral + '_r']
            lRxns.remove(rGeneral)
            return model, [rGeneral, rGeneral + '_r'], lRxns
        else:
            return model, [], lRxns
    else:
        return model, [], lRxns

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
LOGDIR = workingDirs[4]
FIGUREDIR = workingDirs[5]

start = time.time()
timeStamp = gL.getTimeStamp()


## uniformo reazioni dei 5 modelli splittati
name='ENGRO2_'
model_231 = cb.io.read_sbml_model(os.path.join(MODELDIR, name+'MDAMB231_wIrrRxns_original.xml'))
model_102A = cb.io.read_sbml_model(os.path.join(MODELDIR,name+ 'MCF102A_wIrrRxns_original.xml'))
model_361 = cb.io.read_sbml_model(os.path.join(MODELDIR, name+'MDAMB361_wIrrRxns_original.xml'))
model_SKBR3 = cb.io.read_sbml_model(os.path.join(MODELDIR,name+ 'SKBR3_wIrrRxns_original.xml'))
model_MCF7 = cb.io.read_sbml_model(os.path.join(MODELDIR, name+'MCF7_wIrrRxns_original.xml'))

lRxns_231 = getRxns(model_231)
lRxns_102A = getRxns(model_102A)
lRxns_361 = getRxns(model_361)
lRxns_SKBR3 = getRxns(model_SKBR3)
lRxns_MCF7 = getRxns(model_MCF7)

lComb = [lRxns_231, lRxns_102A, lRxns_361, lRxns_SKBR3, lRxns_MCF7]
lCombAll = itt.permutations(lComb, 2)
i = 0
lAllRxns = []
for comb in lCombAll:
    lDiffAB = gL.difference(comb[0],comb[1])
    lDiffBA = gL.difference(comb[1], comb[0])
    print('lDiffAB\t', lDiffAB)
    print('lDiffBA\t', lDiffBA, '\n')
    lAllRxns += lDiffAB
    lAllRxns += lDiffBA
    i+= 1

lAllRxns = gL.unique(lAllRxns)

lAllRxns = [r for r in lAllRxns if r.endswith('_f') == True or r.endswith('_b') == True or r.endswith('_r') == True]

idx = 0
while len(lAllRxns) != 0 and idx < len(lAllRxns):
    print('current Rxn:\t', lAllRxns[idx])
    if lAllRxns[idx].endswith('_f') or lAllRxns[idx].endswith('_b'):
        model_231, l2Remove1, lRxns_231 = addbackwardForwardRxn(lAllRxns[idx], lRxns_231, model_231)
        model_102A, l2Remove2, lRxns_102A = addbackwardForwardRxn(lAllRxns[idx], lRxns_102A, model_102A)
        model_361, l2Remove3, lRxns_361 = addbackwardForwardRxn(lAllRxns[idx], lRxns_361, model_361)
        model_SKBR3, l2Remove4, lRxns_SKBR3 = addbackwardForwardRxn(lAllRxns[idx], lRxns_SKBR3, model_SKBR3)
        model_MCF7, l2Remove5, lRxns_MCF7 = addbackwardForwardRxn(lAllRxns[idx], lRxns_MCF7, model_MCF7)
        l2Remove = gL.unique(l2Remove1 + l2Remove2 + l2Remove3 + l2Remove4 + l2Remove5)
        lAllRxns = gL.difference(lAllRxns, l2Remove)
        idx = 0
    elif lAllRxns[idx].endswith('_r'):
        model_231, l2Remove1, lRxns_231 = addReverseRxn(lAllRxns[idx], lRxns_231, model_231)
        model_102A, l2Remove2, lRxns_102A= addReverseRxn(lAllRxns[idx], lRxns_102A, model_102A)
        model_361, l2Remove3, lRxns_361 = addReverseRxn(lAllRxns[idx], lRxns_361, model_361)
        model_SKBR3, l2Remove4, lRxns_SKBR3 = addReverseRxn(lAllRxns[idx], lRxns_SKBR3, model_SKBR3)
        model_MCF7, l2Remove5, lRxns_MCF7 = addReverseRxn(lAllRxns[idx], lRxns_MCF7, model_MCF7)
        l2Remove = gL.unique(l2Remove1 + l2Remove2 + l2Remove3 + l2Remove4 + l2Remove5)
        lAllRxns = gL.difference(lAllRxns, l2Remove)
        idx = 0
    else:
        idx += 1

cb.io.write_sbml_model(model_231, os.path.join(MODELDIR, name+'MDAMB231_wIrrRxns.xml'))
cb.io.write_sbml_model(model_102A, os.path.join(MODELDIR, name+'MCF102A_wIrrRxns.xml'))
cb.io.write_sbml_model(model_361, os.path.join(MODELDIR, name+'MDAMB361_wIrrRxns.xml'))
cb.io.write_sbml_model(model_SKBR3, os.path.join(MODELDIR,name+ 'SKBR3_wIrrRxns.xml'))
cb.io.write_sbml_model(model_MCF7, os.path.join(MODELDIR, name+'MCF7_wIrrRxns.xml'))