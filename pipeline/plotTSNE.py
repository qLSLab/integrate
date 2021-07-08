import pandas as pd
import os
import genericLib as gL
import dataVizLib as dvl
import sys

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]
FIGUREDIR = workingDirs[4]

# setting input data
tStampRS = sys.argv[1]
nExtraction2Plot = 0
lcellLines = ['MCF102A', 'SKBR3', 'MCF7', 'MDAMB231', 'MDAMB361']

inTsneFile = 'tsnewPca_rs_' + tStampRS + '_extr_' + str(nExtraction2Plot) + '.csv'
basenameTsne = inTsneFile[:-4]
dfTsne = pd.read_csv(os.path.join(OUTDIR, inTsneFile), sep='\t', index_col=0)
dfTsne = dfTsne.rename(columns = {'cellLine':'sampleID'})
dfTsne['cellLine'] = dfTsne.sampleID.str.split('_', expand=True)[0]
dfTsne = dfTsne.set_index('sampleID')
tagFile = os.path.join(OUTDIR, basenameTsne + '.tags')
outFileName = os.path.join(FIGUREDIR, basenameTsne + '_tags.png')
dvl.tagsOnTsne(dfTsne, 'cellLine', dvl.dColors, tagFileName=tagFile,
               file2save=outFileName)
