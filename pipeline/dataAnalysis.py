#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import genericLib as gL
import integrateLib as iL
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sn
import matplotlib.patches as mpatches
import scanpy as sc
import ast
from sklearn import preprocessing
from sklearn.manifold import TSNE
from scanpy import AnnData

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]

# setting input data
cellProteinFile = 'CELLS-PROTEINS.csv'
proteinCellsOutputFigure = 'ProteinsvsCells.png'
cellsTimeOutputFigure = 'CellsVsTime.png'
proteinsTimeOutputFigure = 'ProteinsVsTime.png'
metsEngroVsMetabolomicsFile = 'metsEngroVsMetabolomics.csv'
metabolomicsLMFile = 'metabolomics_LM.csv'
tsneFigure = 'tsne.png'
dotplotFigure = 'dotplot.png'

# Set the color to draw the cell line
dColors = {
    'MCF102A': '#029E73',
    'SKBR3': '#D55E00',
    'MCF7': '#CC78BC',
    'MDAMB231': '#BAB029',
    'MDAMB361': '#0173B2'
}

sn.set(font_scale = 5)

# Load the data
Data=pd.read_csv(os.path.join(RAWDIR, cellProteinFile), sep=";")
Data=Data[Data["time (h)"]!=144]
Data['Line'] = Data['Line'].replace(['MDA-MB-361'],'MDAMB361')
Data['Line'] = Data['Line'].replace(['MDA-MB-231'],'MDAMB231')

DataProt=Data[Data["Feature"]!="CELLS"]
DataCells=Data[Data["Feature"]=="CELLS"]

DataProt2=DataProt[DataProt["time (h)"]!=72]
DataCells2=DataCells[DataCells["time (h)"]!=72]

namesLine=np.unique(DataProt["Line"])
j=0
for name in namesLine:
    A=list()
    B=list()
    df=pd.DataFrame()
    for i in range(1,4):
        A.extend(DataProt2[str(i)][DataProt['Line']==name].values)
        B.extend(DataCells2[str(i)][DataCells['Line']==name].values)
    df["Proteins"]=A
    df["Cells"]=B
    df["Line"]=name
    if j==0:
        df_total=df.copy()
    else:
        df_total=df_total.append(df,ignore_index=True)
    j=j+1


# Plot protein vs cells
i=0
sn.set(font_scale=5)
f, axes = plt.subplots(1, 5,figsize=(90,15))#,sharey=True)
for name in namesLine:
   ax=sn.regplot(x="Proteins",y="Cells",data=df_total[df_total["Line"]==name],ax=axes[i])
   ax=sn.scatterplot(x="Proteins",y="Cells",hue="Line",
                  data=df_total[df_total["Line"]==name],ax=axes[i],
                  palette={name:dColors[name]},legend=False,s=1000)
   ax.set_title(name,fontsize = 50)
   ax.ticklabel_format(axis='y',style='sci',scilimits=(1,5))
   ax.set(xlabel="Protein [$\mu$g]",ylabel="N° of cells")
   i=i+1

f.savefig(os.path.join(FIGUREDIR, proteinCellsOutputFigure))

# N° of cells vs Time
i=0
sn.set(font_scale=2)
f=plt.figure(figsize=(12,8))
for name in namesLine:
   A=list()
   for j in range(1,4):
        ax=sn.scatterplot(x="time (h)",y=str(j),hue="Line",
                      data=DataCells[DataCells["Line"]==name],
                      palette={name:dColors[name]},legend=True,s=100)
        A.append(DataCells[DataCells["Line"]==name][str(j)])
   A=np.array(A)
   sn.lineplot(x=[0,24,48,72],y=[np.mean(x) for x in A.T],color=dColors[name],legend=True)
   i=i+1

ax.set(xlabel="Time (h)",ylabel="N° of cells")
ax.ticklabel_format(axis='y',style='sci',scilimits=(1,5))

handles=list()
for name in dColors.keys():
    handles.append(mpatches.Patch(color=dColors[name], label=name))
ax.legend(handles=handles)
f.savefig(os.path.join(FIGUREDIR, cellsTimeOutputFigure))

# Protein vs Time
i=0
sn.set(font_scale=2)
f=plt.figure(figsize=(12,8))
for name in namesLine:
   A=list()
   for j in range(1,4):
        ax=sn.scatterplot(x="time (h)",y=str(j),hue="Line",
                      data=DataProt[DataProt["Line"]==name],
                      palette={name:dColors[name]},legend=True,s=100)
        A.append(DataProt[DataProt["Line"]==name][str(j)])
   A=np.array(A)
   sn.lineplot(x=[0,24,48,72],y=[np.mean(x) for x in A.T],color=dColors[name],legend=True)
   i=i+1

ax.set(xlabel="Time (h)",ylabel="Protein [$\mu$g]")
ax.ticklabel_format(axis='y',style='sci',scilimits=(1,5))

handles=list()
for name in dColors.keys():
    handles.append(mpatches.Patch(color=dColors[name], label=name))
ax.legend(handles=handles)

f.savefig(os.path.join(FIGUREDIR, proteinsTimeOutputFigure))

# Cluster analysis
# load conversion table between id metabolites in metamolomics data and id metabolites in ENGRO2 model and create dictionary
metIDConversion=pd.read_csv(os.path.join(RAWDIR, metsEngroVsMetabolomicsFile), sep=';')
metIDConversion['metId_engro'] = metIDConversion['metId_engro'].apply(ast.literal_eval)

#create dictionary Mass Spectometry to ENGRO
M_engro_dict=metIDConversion.set_index('metId_M')['metId_engro'].to_dict()

#create dictionary ENGRO to Mass Spectometry
engro_M_dict=iL.reverse_dict(M_engro_dict)

# Cluster analysis on metabolomic data
datasetClustering=pd.read_csv(os.path.join(RAWDIR, metabolomicsLMFile), sep=';',index_col=0)

lista=[el for el in datasetClustering.columns if len(M_engro_dict[el])>0]
datasetClustering=datasetClustering.loc[:,lista]
print(datasetClustering.shape)
# Normalize data
scaler = preprocessing.MaxAbsScaler().fit(datasetClustering)
X_scaled = scaler.transform(datasetClustering)

# Tsne computation
X_embedded = TSNE(n_components=2,random_state=0).fit_transform(X_scaled)

# Display the results of tsne
sn.set_style("darkgrid")
j=0
fig=plt.figure(figsize=(15,10))
for key in dColors.keys():
    elements=[i+18*j for i in range(18)]
    plt.scatter(X_embedded[elements,0],X_embedded[elements,1],s=300,label=key,color=dColors[key])
    j=j+1

plt.legend()
plt.tick_params(
axis='y', # changes apply to the x-axis
which='both', # both major and minor ticks are affected, # ticks along the bottom edge are off
left=False, # ticks along the top edge are off
labelleft=False,) # labels along the bottom edge are off

plt.tick_params(
axis='x', # changes apply to the x-axis
which='both', # both major and minor ticks are affected, # ticks along the bottom edge are off
bottom=False, # ticks along the top edge are off
labelbottom=False,) # labels along the bottom edge are off

plt.xlabel("TSNE 1")
plt.ylabel("TSNE 2")
fig.savefig(os.path.join(FIGUREDIR,tsneFigure))

# Dotplot
adata=AnnData(X_scaled)
adata.var.index=list(datasetClustering.columns)
adata.obs["clusters"]=[key for key in dColors.keys() for i in range(18) ]
adata.raw=adata

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack',n_comps=50)
sc.tl.rank_genes_groups(adata, "clusters", method='t-test')

df_markers=pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
df_marker_list=df_markers.T.values.flatten()

mpl.rcParams.update(mpl.rcParamsDefault)

dp=sc.pl.dotplot(adata, df_marker_list, groupby="clusters",expression_cutoff =0,
               swap_axes =True, dendrogram=False,return_fig=True,
                size_title ="% of samples above mean",
                 #cmap ="binary",
                 use_raw=False,
                colorbar_title ="Normalized Mean")

dp.savefig(os.path.join(FIGUREDIR, dotplotFigure))
