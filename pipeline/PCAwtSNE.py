#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import genericLib as gL
import sys
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import preprocessing
import random
import integrateLib as iL

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
tStampRS = sys.argv[1]
modelId = 'ENGRO2'
nExtractions = 1
nExtracted = 10000
nSamples = 1000000
nComponents = 10
niterations = 2000
lcellLines = ['MCF102A', 'SKBR3', 'MCF7','MDAMB231', 'MDAMB361']

## load net fluxes from random sampling output files
dRS = {}
dIdx = {}

for c in lcellLines:
    inputName = 'randomSampling_' + modelId + '_nSol_' + str(nSamples) + '_' + cellLine + '_' + tStampRS + '.csv'
    df = pd.read_csv(os.path.join(OUTDIR, inputName), sep = '\t')
    dfNet = iL.splt2Net(df)
    dRS[c] = dfNet
    dIdx[c] = dfNet.index.tolist()

for extr in range(0, nExtractions):
    lOFdf = []
    for cellLine in lcellLines:
        lIdx2Extract = random.sample(dIdx[cellLine], nExtracted)
        lIdx2Extract.sort()
        dIdx[cellLine] = gL.difference(dIdx[cellLine], lIdx2Extract)
        df_subset = dRS[cellLine].iloc[lIdx2Extract]
        lIndexes = df_subset.index.tolist()
        lNewIndexes = []
        for idx in lIndexes:
            lNewIndexes.append(cellLine + '_' + str(idx))
        df_subset.index = lNewIndexes
        lOFdf.append(df_subset)

    allDf_OFdf = pd.concat(lOFdf)
    npallDf_OFdf = allDf_OFdf.to_numpy()

    ## Execute PCA
    scaler = preprocessing.StandardScaler().fit(npallDf_OFdf)
    X_scaled = scaler.transform(npallDf_OFdf)

    pca = PCA(n_components = nComponents)
    Xpca = pca.fit(X_scaled).transform(X_scaled)

    lComponents = ['PC' + str(el) for el in range(0,nComponents)]
    dfXpca_completo = pd.DataFrame(data = Xpca, columns = lComponents)
    lClustersName = []
    for cellLine in lcellLines:
        lClustersName += [cellLine] * nExtracted
    dfXpca_completo['Cluster'] = lClustersName

    dflabel = dfXpca_completo['Cluster']
    lComponents = ['PC' + str(el) for el in range(0,nComponents)]
    dfXpca = dfXpca_completo[lComponents]

    ## Define perplexity
    perp = (dfXpca.shape[0])**0.5

    ## Execute t-SNE
    tsne = TSNE(n_components=2, verbose=2, perplexity = perp, n_iter=niterations, n_jobs = -1, init = 'pca', random_state = 0)
    tsne_results = tsne.fit_transform(dfXpca)
    df_subset = pd.DataFrame({'cellLine': dflabel})
    df_subset['tsne-2d-one'] = tsne_results[:,0]
    df_subset['tsne-2d-two'] = tsne_results[:,1]

    df_subset.to_csv(os.path.join(OUTDIR, 'tsnewPca_rs_' + tStampRS + '_extr_' + str(extr) + '.csv'), sep = '\t', index = True)

    ## Save the nExtracted fluxes from the initial random sampling set
    allDf_OFdf.to_csv(os.path.join(OUTDIR, 'sampledFluxes_rs_' + tStampRS + '_extr_' + str(extr) + '.csv'), sep = '\t', index = True)
