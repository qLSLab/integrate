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

## set dictionary for coloring each cell lines within plot
dColors = {
    'MCF102A': '#029e73',
    'SKBR3': '#d55e00',
    'MCF7': '#cc78bc',
    'MDAMB231': '#ece133',
    'MDAMB361': '#0173b2'
}

## compute computational glucose yield
lcellLines = ['MCF102A', 'SKBR3', 'MCF7','MDAMB231', 'MDAMB361']

biomassRxnName = 'Biomass'
glc = 'EX_glc__D_e_r'
gln = 'EX_gln__L_e_r'

dyield_silico_glc_MEDIAN = {}
for c in lcellLines:
    print(c)
    dfRS = pd.read_csv(os.path.join(OUTDIR, 'randomSampling_ENGRO2_nSol_10000_nBatch_0_' + c + '_MediumANDYsiANDRas.csv'), sep = '\t')
    dfRS['yield_glc'] = (0.131972 * dfRS[biomassRxnName]) / ((dfRS[glc]*1e-3) * 180.1559)
    yield_glc_median = dfRS['yield_glc'].median()
    dyield_silico_glc_MEDIAN[c + '_A'] = yield_glc_median
    dyield_silico_glc_MEDIAN[c + '_B'] = yield_glc_median


dfyield_silico_glc = pd.DataFrame(dyield_silico_glc_MEDIAN.items(), columns = ['cellLine', 'yield_glc_silico'])
dfyield_silico_glc.to_csv(os.path.join(OUTDIR, 'engro_yield_glc_MEDIAN_' + tStampRS + '.csv'), sep = '\t', index = True)

## compute experimental glucose yield
dProtContent = {'MCF102A': [175.7910906 - 40.84485407, 110.4608295 - 36.09831029], 'SKBR3': [74.76190476 - 45.42242704, 111.3056836 - 27.52688172], 'MCF7': [87.92626728 - 27.88018433, 87.86482335 - 31.1827957], 'MDAMB231': [95.6374808 - 39.46236559, 89.8156682 - 43.0890937], 'MDAMB361': [38.66359447- 27.75729647, 40.07680492 -28.57142857]} # t48-t0 del YSI
dGlcConsumption = {'MCF102A': [4.458855544, 4.654227382], 'SKBR3': [2.917178061, 6.764445483], 'MCF7': [2.927506502, 2.958492487], 'MDAMB231': [6.759204908, 6.637699054], 'MDAMB361': [0.414875385, 0.407950783]} # la terza misurazione di ogni linea cellulare è la media dei precedenti due
dGlnConsumption = {'MCF102A': [1.020, 0.970], 'SKBR3': [0.419, 1.032], 'MCF7': [1.219, 1.364], 'MDAMB231': [1.728, 1.560], 'MDAMB361': [0.946, 0.967]} # la terza misurazione di ogni linea cellulare è la media dei precedenti due

dyield_wet_glc = {}
for c in lcellLines:
    dyield_wet_glc[c + '_A'] = (dProtContent[c][0]*1e-6) / ((dGlcConsumption[c][0]*1e-3) * 180.1559)
    dyield_wet_glc[c + '_B'] = (dProtContent[c][1]*1e-6) / ((dGlcConsumption[c][1]*1e-3) * 180.1559)

dfyield_wGlc = pd.DataFrame(dyield_wet_glc.items(), columns = ['cellLine', 'yield_glc_wet'])

## merge experimental and computational yields and plot results
dfyield_silico_glc = pd.read_csv(os.path.join(OUTDIR, 'engro_yield_glc_MEDIAN_' + tStampRS + '.csv'), sep = '\t', index_col=0, dtype = {'yield_glc_silico': float, 'cellLine': str})
dfY_rl = pd.merge(dfyield_silico_glc, dfyield_wGlc, on = 'cellLine')
dfY_rl = dfY_rl.set_index('cellLine')
dfY_rl_corr = dfY_rl.corr(method='pearson')
dfY_rl_corr_np = np.corrcoef(dfY_rl['yield_glc_silico'], dfY_rl['yield_glc_wet'])

fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(111)
ax.set_axisbelow(True)
plt.grid(color='w', linestyle='solid')
ax.set_facecolor('#eaeaf2')
plt.ticklabel_format(axis='y', style='scientific', scilimits=(-4,-4))
plt.ticklabel_format(axis='x', style='scientific', scilimits=(-4,-4))
ax.yaxis.get_offset_text().set_fontsize(35)
ax.xaxis.get_offset_text().set_fontsize(35)
for spine in ax.spines.values():
    spine.set_visible(False)
color_dict = dict(dColors)
ax = sns.regplot(x='yield_glc_silico', y ='yield_glc_wet', data = dfY_rl, ci = None)
for row in dfY_rl.itertuples():
    try:
        ax.scatter(row.yield_glc_silico, row.yield_glc_wet, marker='o', color = dColors[row.Index[:-2]], label = row.Index[:-2], s = 400)
    except:
        ax.scatter(row.yield_glc_silico, row.yield_glc_wet, marker='o', color = dColors[row.Index], label = row.Index, s = 400)
#regression part
slope, intercept, r_value, p_value, std_err = stats.linregress(dfY_rl['yield_glc_silico'],dfY_rl['yield_glc_wet'])
r_value_s = spearmanr(dfY_rl['yield_glc_silico'],dfY_rl['yield_glc_wet'])
ax.set_title('Correlation = ' + "{:.2f}".format(r_value_s.correlation), fontsize=55)
plt.xlabel('Computational yield',fontsize= 45)
plt.ylabel('Experimental yield',fontsize= 45)
plt.xlim([min(list(dfY_rl['yield_glc_silico'])) - ((min(list(dfY_rl['yield_glc_silico']))*5)/100), max(list(dfY_rl['yield_glc_silico'])) + ((max(list(dfY_rl['yield_glc_silico']))*2)/100)])
plt.xticks(fontsize= 35)
plt.yticks(fontsize= 35)
plt.legend(fontsize= 30)
fig.subplots_adjust(left=0.2)
plt.savefig(os.path.join(FIGUREDIR, 'yieldENGRO_wetVsSilico_Glc_MEDIANSpearman_' + tStampRS + '.pdf'))
plt.close()
