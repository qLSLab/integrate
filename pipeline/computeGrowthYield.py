#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import genericLib as gL
import os
import sys
import numpy as np
import seaborn as sns

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]

tStampRS = sys.argv[1]
whichAnalysis = sys.argv[2]

if whichAnalysis == 'save':
    #### Step 1: SAVE DATA
    ##### Compute in silico Growth Yield
    dyield_silico_glc = {}
    dfMCF102ANet = pd.read_csv(os.path.join(OUTDIR, 'randomSampling_ENGRO2_1000000_MCF102A_' + tStampRS + '.csv'), sep = '\t')
    dfMCF102ANet['yield_glc'] = dfMCF102ANet['Biomass'] / dfMCF102ANet['EX_glc__D_e_r']
    yield_glc_mean = dfMCF102ANet['yield_glc'].mean()
    dyield_silico_glc['MCF102A_A'] = yield_glc_mean
    dyield_silico_glc['MCF102A_B'] = yield_glc_mean

    lcellLines = ['SKBR3', 'MCF7','MDAMB231', 'MDAMB361']
    for c in lcellLines:
        dfRS = pd.read_csv(os.path.join(OUTDIR, 'randomSampling_ENGRO2_1000000_' + c +'_' + tStampRS + '.csv'), sep = '\t')
        dfRS['yield_glc'] = dfRS['Biomass'] / dfRS['EX_glc__D_e_r']
        yield_glc_mean = dfRS['yield_glc'].mean()
        dyield_silico_glc[c + '_A'] = yield_glc_mean
        dyield_silico_glc[c + '_B'] = yield_glc_mean

    dfyield_silico_glc = pd.DataFrame(dyield_silico_glc.items(), columns = ['cellLine', 'yield_glc_silico'])

    ##### Compute wet Growth Yield
    lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']
    dfCellCount = pd.read_csv(os.path.join(RAWDIR, 'cellCounts.csv'), sep = '\t')
    dfCellCount = dfCellCount.drop(dfCellCount[dfCellCount.Time > 48].index)

    dfProteinContent_w = pd.read_csv(os.path.join(RAWDIR, 'proteinContent.csv'), sep = '\t')
    dfProteinContent_w = dfProteinContent_w.drop(dfProteinContent_w[dfProteinContent_w.Time > 48].index)

    # dati ricavati stessi tempi del YSI
    dProtContent = {'MCF102A': [175.7910906298, 110.460829493088], 'SKBR3': [74.7619047619048, 111.305683563748], 'MCF7': [87.926267281106, 87.8648233486943], 'MDAMB231': [95.6374807987711, 89.815668202765], 'MDAMB361': [38.6635944700461, 40.0768049155146]}
    dGlcConsumption = {'MCF102A': [4.458855544, 4.654227382], 'SKBR3': [2.917178061, 6.764445483], 'MCF7': [2.927506502, 2.958492487], 'MDAMB231': [6.759204908, 6.637699054], 'MDAMB361': [0.414875385, 0.407950783]}
    dGlnConsumption = {'MCF102A': [1.020263396, 0.969936553], 'SKBR3': [0.418782880, 1.032064126], 'MCF7': [1.219061914, 1.364396060], 'MDAMB231': [1.727965579, 1.560298921], 'MDAMB361': [0.946320027, 0.967281240]}

    dyield_wet_glc = {}
    for c in lcellLines:
        dyield_wet_glc[c + '_A'] = ((dProtContent[c][0]*1e-6) / ((dGlcConsumption[c][0]*1e-3) * 180.1559))
        dyield_wet_glc[c + '_B'] = ((dProtContent[c][1]*1e-6) / ((dGlcConsumption[c][1]*1e-3) * 180.1559))

    dfyield_wGlc = pd.DataFrame(dyield_wet_glc.items(), columns = ['cellLine', 'yield_glc_wet'])

    ## Join the computed wet and computational growth yield
    dfY_rl = pd.merge(dfyield_silico_glc, dfyield_wGlc, on = 'cellLine')
    dfY_rl = dfY_rl.set_index('cellLine')
    dfY_rl.to_csv(os.path.join(OUTDIR, 'yieldOnGlc_rs_' + tStampRS + '.csv'), sep = '\t', index = True)


elif whichAnalysis == 'plot':
    #### Step 1: PLOT DATA
    dColors = {
        'MCF102A': '#FF7F00',
        'SKBR3': '#0080FF',
        'MCF7': '#FF0000',
        'MDAMB231': '#4DAF4A',
        'MDAMB361': '#CE00FF'
    }

    ## corr plot metodo linear regression calcolo yield
    dfY_rl = pd.read_csv(os.path.join(OUTDIR, 'yieldOnGlc_rs_' + tStampRS + '.csv'), sep = '\t', index_col = 0)
    dfY_rl_corr_np = np.corrcoef(dfY_rl['yield_glc_silico'], dfY_rl['yield_glc_wet'])

    fig = plt.figure(figsize=(15,12))
    ax = fig.add_subplot(111)
    ax.set_axisbelow(True)
    plt.grid(color='w', linestyle='solid')
    ax.set_facecolor('#eaeaf2')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(-4,-4))
    ax.yaxis.get_offset_text().set_fontsize(40)
    # hide axis spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    color_dict = dict(dColors)
    ax = sns.regplot(x='yield_glc_silico', y ='yield_glc_wet', data = dfY_rl, ci = None)

    for row in dfY_rl.itertuples():
        try:
            ax.scatter(row.yield_glc_silico, row.yield_glc_wet, marker='o', color = dColors[row.Index[:-2]], label = row.Index[:-2], s = 400)
        except:
            ax.scatter(row.yield_glc_silico, row.yield_glc_wet, marker='o', color = dColors[row.Index], label = row.Index, s = 400)

    ax.set_title('Correlation = ' + "{:.2f}".format(dfY_rl_corr_np[0,1]), fontsize=55)
    plt.xlabel('Computational yield',fontsize= 45)
    plt.ylabel('Experimental yield',fontsize= 45)
    plt.xlim([min(list(dfY_rl['yield_glc_silico'])) - ((min(list(dfY_rl['yield_glc_silico']))*5)/100), max(list(dfY_rl['yield_glc_silico'])) + ((max(list(dfY_rl['yield_glc_silico']))*2)/100)])
    plt.xticks(fontsize= 30)
    plt.yticks(fontsize= 30)
    plt.legend(fontsize= 30)
    fig.subplots_adjust(left=0.2)
    plt.savefig(os.path.join(FIGUREDIR, 'yieldOnGlc_rs_' + tStampRS + '.png'))
