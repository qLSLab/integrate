# coding: utf-8

import os
import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import cobra as cb
import sklearn.metrics
from matplotlib.colors import ListedColormap
#from adjustText import adjust_text
from matplotlib.patches import Rectangle
import genericLib as gL

name="cohen"

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]

# setting input data
val = np.log2(1.2) # threshold
weight = "linear" # weight on the Kappa Cohen
resultsMetabolomicFile = 'resultsMetabolomic'
# name_file="outputs/resultsMetabolomic" # Directory of the data
# dir_output="outputs/"
rxn2Visualize = 'ACONT'
model_name = 'ENGRO2_irrev.xml'
namesLine=['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']
biomassRxn = 'Biomass'
mwuTestFile = 'mwuTest_'
rasFile = 'ras_'
dfconcordanceFile = 'df_concordance.csv'
concordanceHeatmapFigure = 'concordanceHeatmap.png'
concordanceHeatmaptotalFigure = 'concordanceHeatmap_total.png'
scatterplotFigure = 'scatterplot.png'
heatmapMeansFigure = 'heatmap_means.png'
meansFile = 'medie_Met.csv' # means
eps=0.00001

# Set all the combinations for the tests
tests=itertools.combinations_with_replacement(namesLine,2)
tests=[(test[0],test[1]) for test in tests if test[0]!=test[1]]

# Load and create the concordance datasets
# Load the metabolitc model
model = cb.io.read_sbml_model(os.path.join(MODELDIR, model_name))
keys=[reaction.id for reaction in model.reactions]
formule=[reaction.reaction for reaction in model.reactions]
dict_formule={key:value for key,value in zip(keys,formule)}
GPR_reactions=[reaction.id for reaction in model.reactions if reaction.gene_reaction_rule!='']

# Load RAS and FBA data
datasetsFBA=list()
datasetsRAS=list()
for test in tests:
    datasetsFBA.append(pd.read_csv(os.path.join(OUTDIR, mwuTestFile + test[0] + '_vs_' + test[1] + '.csv'),
                                sep="\t",index_col=0).sort_index())
    datasetsRAS.append(pd.read_csv(os.path.join(OUTDIR, rasFile + test[0] + '_vs_' + test[1] + '.csv'),
                                sep="\t",index_col=0).loc[GPR_reactions])

# Create the dataset of concordance and median/mean for FBA and RAS
df_concFBA=pd.DataFrame(index=datasetsFBA[0].index)
df_medianFBA=pd.DataFrame(index=datasetsFBA[0].index)

df_concRAS=pd.DataFrame(index=datasetsRAS[0].index)
df_meanRAS=pd.DataFrame(index=datasetsRAS[0].index)

# Perform statistical tests
for datasetFBA,datasetRAS,test in zip(datasetsFBA,datasetsRAS,tests):
    #RAS
    df_meanRAS[str(test)]=np.log2((datasetRAS["mean_"+test[0]]+eps)/(datasetRAS["mean_"+test[1]]+eps))#
    df_concRAS[str(test)]=datasetRAS["result"]
    df_medianFBA[str(test)]=np.log2((datasetFBA["median_"+test[0]]+eps)/(datasetFBA["median_"+test[1]]+eps))##
    pvalue_less=datasetFBA["pvalue_less"]
    pvalue_great=datasetFBA["pvalue_greater"]

    lista=list()
    for el,el2 in zip(pvalue_less,pvalue_great):
        if el<0.05:
            lista.append(-1)
        elif el2<0.05:
            lista.append(1)
        else:
            lista.append(0)
    df_concFBA[str(test)]=lista

# Load the dataset of mean and concordance for MET
df_met=pd.read_csv(os.path.join(OUTDIR, resultsMetabolomicFile + '_ttest.csv'), index_col = 0)
df_met_mean=pd.read_csv(os.path.join(OUTDIR, resultsMetabolomicFile + '_mean.csv'), index_col = 0).T

# Set to 0 the mean below specific threshold
df_met[np.abs(df_met_mean)<val]=0
df_concFBA[np.abs(df_medianFBA)<val]=0
df_concRAS[np.abs(df_concRAS)<val]=0

# Create index of Common reactions among the concordance datasets
index_commonMETvsFBA=df_met.index.intersection(df_concFBA.index)
index_commonMETvsRAS=df_met.index.intersection(df_concRAS.index)
index_commonFBAvsRAS=df_concFBA.index.intersection(df_concRAS.index)
indexCommon=index_commonMETvsFBA.intersection(index_commonMETvsRAS)

# Create the indexes of all the reactions in the common datasets
indexUnion=df_met.index.union(df_concRAS.index).union(df_concFBA.index)

# Sort the columns for comparison
df_concRAS2=df_concRAS[np.sort(df_concRAS.columns)]
df_meanRAS2=df_meanRAS[np.sort(df_meanRAS.columns)]

df_concFBA2=df_concFBA[np.sort(df_concFBA.columns)]
df_medianFBA2=df_medianFBA[np.sort(df_medianFBA.columns)]

df_met2=df_met[np.sort(df_met.columns[:-1])]
df_met_mean2=df_met_mean[np.sort(df_met_mean.columns[:-1])]

############################################################################
#######################statistical analysis
from statsmodels.stats.multitest import fdrcorrection,multipletests
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import kstest 
from statsmodels.graphics.gofplots import qqplot_2samples

ntest = 100

dfConcFRD_FBAvsMET=pd.DataFrame(index=index_commonMETvsFBA,columns=[str(i) for i in range(ntest)])
dfConcMETvsFBA=pd.DataFrame(index=index_commonMETvsFBA,columns=["METvsFBA"])

dfConcFRD_RASvsMET=pd.DataFrame(index=index_commonMETvsRAS,columns=[str(i) for i in range(ntest)])
dfConcRASvsMET=pd.DataFrame(index=index_commonMETvsRAS,columns=["RASvsMET"])

dfConcFRD_RASvsFBA=pd.DataFrame(index=index_commonFBAvsRAS,columns=[str(i) for i in range(ntest)])
dfConcRASvsFBA=pd.DataFrame(index=index_commonFBAvsRAS,columns=["RASvsFBA"])

#############################FDR for FBA vs MET and MET vs RAS and FBA vs RAS
df_met3=df_met2.copy()
valuesMet=df_met3.values.ravel()
#df_concRAS3=df_concRAS2.copy()
#valuesRas=df_concRAS3.values.ravel()

#kappa random
for ind in range(ntest):
    np.random.seed(ind)
    valori=np.random.choice(valuesMet,df_met3.shape,replace=False)
    df_met3.iloc[:,:]=valori
    #valori=np.random.choice(valuesRas,df_concRAS3.shape,replace=False)
    #df_concRAS3.iloc[:,:]=valori

    for reaction in index_commonMETvsFBA:
        valori1=list(df_met3.loc[reaction,:].values)
        valori2=list(df_concFBA2.loc[reaction,:].values)
        valori1.extend([-el for el in valori1])
        valori2.extend([-el for el in valori2])
        dfConcFRD_FBAvsMET.loc[reaction,str(ind)]=sklearn.metrics.cohen_kappa_score(valori1,valori2,labels=[-1,0,1],weights="linear",sample_weight=None)
    
    for reaction in index_commonMETvsRAS:
        valori1=list(df_met3.loc[reaction,:].values)
        valori2=list(df_concRAS2.loc[reaction,:].values)
        valori1.extend([-el for el in valori1])
        valori2.extend([-el for el in valori2])

        dfConcFRD_RASvsMET.loc[reaction,str(ind)]=sklearn.metrics.cohen_kappa_score(valori1,valori2,labels=[-1,0,1],weights="linear",sample_weight=None)
     
        
#######################################àà 
#kappa true
for reaction in index_commonMETvsFBA:
    valori1=list(df_met2.loc[reaction,:].values)
    valori2=list(df_concFBA2.loc[reaction,:].values)
    valori1.extend([-el for el in valori1])
    valori2.extend([-el for el in valori2])

    dfConcMETvsFBA.loc[reaction,:]=sklearn.metrics.cohen_kappa_score(valori1,valori2,labels=[-1,0,1],weights="linear",sample_weight=None)
    
for reaction in index_commonMETvsRAS:
    valori1=list(df_met2.loc[reaction,:].values)
    valori2=list(df_concRAS2.loc[reaction,:].values)
    valori1.extend([-el for el in valori1])
    valori2.extend([-el for el in valori2])

    dfConcRASvsMET.loc[reaction,:]=sklearn.metrics.cohen_kappa_score(valori1,valori2,labels=[-1,0,1],weights="linear",sample_weight=None)

###############################Plot of Cohen distribution (random vs experiment)
#Cumulative distribution
ecdfRPSvsFBA=ECDF(dfConcFRD_FBAvsMET.values.ravel())
ecdfRPSvsFBA_result=ECDF(dfConcMETvsFBA["METvsFBA"].values)
ecdfRASvsMET=ECDF(dfConcFRD_RASvsMET.values.ravel())
ecdfRASvsMET_result=ECDF(dfConcRASvsMET["RASvsMET"].values)

#####figure of cumulative distributions
#engro2
fig, ax=plt.subplots(1,1,figsize=(20, 10))
qqplot_2samples(dfConcFRD_FBAvsMET.values.ravel(),
                     dfConcMETvsFBA["METvsFBA"].values,
                    ylabel="Agreement by chance",
                    xlabel="INTEGRATE agreement",
                    line="45",
                    ax=ax)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.grid()
ax.set_title("RPSvsFFD ")
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)


fig.savefig(os.path.join(FIGUREDIR, "qqplotEngro2.png"),transparent=False,format="png")


##########################################fdr correction
############pvalues MET vs FBA
dfConcMETvsFBA=dfConcMETvsFBA.sort_values(by=dfConcMETvsFBA.columns[0],ascending=True)
dict_pvaluesMETvsFBA=dict()

valori=dfConcMETvsFBA.values
pvaluesRPSvsFBA=[1-ecdfRPSvsFBA(el[0]) for el in valori]
resultsRPSvsFBA_correction=fdrcorrection(pvaluesRPSvsFBA,0.05)
res=resultsRPSvsFBA_correction[1]

i=0
for el in dfConcMETvsFBA.index:
    dict_pvaluesMETvsFBA[el]=res[i]
    i=i+1

#print(dict_pvaluesMETvsFBA)
#for el,reaction in zip(res,dfConcMETvsFBA.index):
#    if el==True:    
#        print(reaction,dfConcMETvsFBA.loc[reaction,"METvsFBA"])

#############pvalues MET vs RAS
dfConcRASvsMET=dfConcRASvsMET.sort_values(by=dfConcRASvsMET.columns[0],ascending=True)
dict_pvaluesRASvsMET=dict()

valori=dfConcRASvsMET.values
pvaluesRASvsMET=[1-ecdfRASvsMET(el[0]) for el in valori]
resultsRASvsMet_correction=fdrcorrection(pvaluesRASvsMET,0.05)
res=resultsRASvsMet_correction[1]

i=0
for el in dfConcRASvsMET.index:
   dict_pvaluesRASvsMET[el]=res[i]
   i=i+1

###################################


#####Figure of a specific reaction
fig=plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size':20})
plt.scatter(df_medianFBA2.loc[rxn2Visualize],df_met_mean2.loc[rxn2Visualize])
plt.grid()
plt.xlabel("FFD fold change")
plt.ylabel("RPS fold change")
plt.xlim([-3.5,3.5])
plt.ylim([-3.5,3.5])
plt.title(rxn2Visualize)
plt.savefig(os.path.join(FIGUREDIR, rxn2Visualize + '_reaction.png'))

# Statistics tests
kappa_values=pd.DataFrame(index=indexUnion,columns=['RPSvsFBA cohen', 'RPSvsRAS cohen','FBAvsRAS cohen',
                                                    'RPSvsFBA pvalue','RPSvsRAS pvalue',#'FBAvsRAS pvalue'
                                                   ])



for r in kappa_values.index:
    compare=pd.DataFrame()
    if r in df_met2.index and r in df_concFBA2.index:
        list1=list(df_met2.loc[r])
        list2=list(df_concFBA2.loc[r])
        compare['MET']=list1+[el*-1 for el in list1]
        compare['FBA']=list2+[el*-1 for el in list2]
        
        kappa_values.loc[r,'RPSvsFBA cohen']=sklearn.metrics.cohen_kappa_score(compare['MET'], 
                                                                               compare['FBA'],
                                                                               labels=[-1,0,1],
                                                                               weights=weight, 
                                                                               sample_weight=None)
        if kappa_values.loc[r,'RPSvsFBA cohen']>=-1:
            kappa_values.loc[r,'RPSvsFBA pvalue']=1-ecdfRPSvsFBA(kappa_values.loc[r,'RPSvsFBA cohen'])
            kappa_values.loc[r,'RPSvsFBA adj-pvalue']=dict_pvaluesMETvsFBA[r]
        
    if r in df_met2.index and r in df_concRAS2.index:
        list1=list(df_met2.loc[r])
        list2=list(df_concRAS2.loc[r])
        compare['MET']=list1+[el*-1 for el in list1]
        compare['RAS']=list2+[el*-1 for el in list2]
             
        kappa_values.loc[r,'RPSvsRAS cohen']=sklearn.metrics.cohen_kappa_score(compare['MET'],
                                                                               compare['RAS'],
                                                                               labels=[-1,0,1],
                                                                               weights=weight,
                                                                               sample_weight=None)
        
        if kappa_values.loc[r,'RPSvsRAS cohen']>=-1:
            kappa_values.loc[r,'RPSvsRAS pvalue']=1-ecdfRASvsMET(kappa_values.loc[r,'RPSvsRAS cohen'])
            kappa_values.loc[r,'RPSvsRAS adj-pvalue']=dict_pvaluesRASvsMET[r]
  
    if r in df_concFBA2.index and r in df_concRAS2.index:
        list1=list(df_concFBA2.loc[r])
        list2=list(df_concRAS2.loc[r])
        compare['FBA']=list1+[el*-1 for el in list1]
        compare['RAS']=list2+[el*-1 for el in list2] 
        
        kappa_values.loc[r,'FBAvsRAS cohen']=sklearn.metrics.cohen_kappa_score(compare['FBA'],compare['RAS'],
                                               labels=[-1,0,1],weights=weight, sample_weight=None)
        #if kappa_values.loc[r,'FBAvsRAS cohen']>=-1:
        #    kappa_values.loc[r,'RPSvsRAS pvalue']=1-ecdfRASvsFBA(kappa_values.loc[r,'RPSvsFBA cohen'])
        
kappa_values=kappa_values.sort_values(by='RPSvsRAS cohen',ascending=False)
kappa_values["Formula"]=[str(dict_formule[el]) for el in kappa_values.index]

#save the concordance dataset
kappa_values.to_csv(os.path.join(OUTDIR, dfconcordanceFile))

####
kappa_values["Formula2"]=[str(el)+':'+str(dict_formule[el]) for el in kappa_values.index]
kappa_values.fillna(-3,inplace=True)  ## just the create a mask map for the next figure

##############################################################################
# Recreate figures of the paper
# Heatmap Cohen
val1=0.2
val2=0.2
fig, ax=plt.subplots(1, 3, figsize=(16, 40), gridspec_kw={'width_ratios': [0.6,2.4, 0.7],'wspace': 0})

df=kappa_values.loc[((kappa_values['RPSvsFBA '+name]>val1)  ),:].sort_values('RPSvsFBA '+name,ascending=False)
df_cohen=df.loc[:,['RPSvsFBA '+name, 'RPSvsRAS '+name]]
mask = df_cohen.values==-3
a=sns.heatmap(df_cohen,annot=True,mask=mask, xticklabels=["RPSvsFFD","RPSvsRAS"],
            cbar=True,
            cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.8),
            yticklabels=[],cmap="viridis",
            ax=ax[2],
            edgecolor="k",
            linewidths=0.1,
            annot_kws={"fontsize":18})

a.set_xticklabels(a.get_xmajorticklabels(), fontsize = 15)

df_new=df
df_new["diff"]=[1 if el2>val2 else 0 for el2 in df_cohen['RPSvsRAS '+name]]
ax[2].xaxis.tick_top()
cmap=ListedColormap(["salmon","grey"])

# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))

# Set alpha
my_cmap[0,-1] = 0.5
my_cmap[1,-1] = 0.3

# Create new colormap
my_cmap = ListedColormap(my_cmap)

###reactions
g=sns.heatmap(np.array(df_new["diff"].values).reshape(df_new["diff"].shape[0],1),ax=ax[0],
             cbar=True,
             cbar_kws = dict(use_gridspec=False,location="bottom", shrink=0.8),
              cmap=my_cmap,
             vmax=1,vmin=0,annot = np.array(df_new["diff"].index).reshape(df_new["diff"].shape[0],1),fmt = '',
             xticklabels=["Reaction"],yticklabels=[],edgecolor="k",linewidths=0.1,annot_kws={"color":"black","fontsize":20})

g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15)

ax[0].xaxis.tick_top()

###formula
g=sns.heatmap(np.array(df_new["diff"].values).reshape(df_new["diff"].shape[0],1),ax=ax[1],
             cbar=True,
             cbar_kws = dict(use_gridspec=False,location="bottom", shrink=0.8),
              cmap=my_cmap,
             vmax=1,vmin=0,annot = np.array(df_new["Formula"].values).reshape(df_new["diff"].shape[0],1),fmt = '',
             xticklabels=["Formula"],yticklabels=[],edgecolor="k",linewidths=0.1,annot_kws={"color":"black","fontsize":15})

g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15)

ax[1].xaxis.tick_top()

fig.savefig(os.path.join(FIGUREDIR, concordanceHeatmapFigure))

############################################################################
# Heatmap of all reactions
fig, ax=plt.subplots(1, 2, figsize=(40, 50), gridspec_kw={'width_ratios': [1, 0.7],'wspace': 0})

df=kappa_values.loc[((kappa_values['RPSvsFBA '+name]!=-3) | (kappa_values['RPSvsRAS '+name]!=-3) ),:].sort_values('RPSvsFBA '+name,ascending=False)

#list of significant reaction
lista=list()
lista_index=list()
for index,el1,el2 in zip(df['RPSvsFBA '+name].index,df['RPSvsFBA '+name],df['RPSvsRAS '+name]):
    if el1>val1 and el2>val2:
        lista.append(2)
        lista_index.append(index)
    elif el1>val1 and el2<=val2:
        lista.append(1)
        lista_index.append(index)
    elif el1<-val1 and el2<-val2:
        lista.append(0)
        lista_index.append(index)
df=df.loc[lista_index,:]

df_cohen=df.loc[:,['RPSvsFBA '+name, 'RPSvsRAS '+name,"RPSvsFBA pvalue","RPSvsRAS pvalue"]] 
mask = df_cohen.values==-3

a=sns.heatmap(df_cohen,annot=True,mask=mask, xticklabels=["RPSvsFFD Cohen","RPSvsRAS Cohen","RPSvsFFD pvalue","RPSvsRAS pvalue"],
            cbar=True,
            cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.8),
            yticklabels=[],cmap="viridis",
            ax=ax[1],edgecolor="k",linewidths=0.1,annot_kws={"fontsize":18})

a.set_xticklabels(a.get_xmajorticklabels(), fontsize = 15)

df_new=df
df_new["diff"]=lista


ax[1].xaxis.tick_top()

cmap=ListedColormap(["white","salmon","grey"])

# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))

# Set alpha
my_cmap[0,-1] = 1
my_cmap[1,-1] = 0.3
my_cmap[2,-1] = 0.3

# Create new colormap
my_cmap = ListedColormap(my_cmap)

g=sns.heatmap(np.array(df_new["diff"].values).reshape(df_new["diff"].shape[0],1),ax=ax[0],
             cbar=True,
             cbar_kws = dict(use_gridspec=False,location="bottom", shrink=0.8),
              cmap=my_cmap,
             vmax=2,vmin=0,annot = np.array(df_new["diff"].index).reshape(df_new["diff"].shape[0],1),fmt = '',
               linecolor='black',
              edgecolor="k",linewidths=0.1,annot_kws={"color":"black","fontsize":20})

g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15)

ax[0].xaxis.tick_top()
fig.savefig(os.path.join(FIGUREDIR, concordanceHeatmaptotalFigure))

################################################################################
# Scatter plot reactions
df=kappa_values.sort_values('RPSvsFBA '+name,ascending=False)
df_cohen=df.loc[:,['RPSvsFBA '+name, 'RPSvsRAS '+name,'FBAvsRAS '+name]]

f = plt.figure(figsize=(13,23))
plt.rcParams.update({'font.size':20})

ax = f.add_subplot(111)
ax.axes.get_yaxis().set_visible(True)
ax.axes.get_xaxis().set_visible(True)

ax.add_patch(Rectangle((val1, val2), 1-val1, -1-val2, facecolor= 'salmon',fill=True,alpha=0.5,zorder=0))
ax.add_patch(Rectangle((val1, val2), 1-val1, 1-val2, facecolor= 'grey',fill=True,alpha=0.3))
ax.add_patch(Rectangle((-val1, -val1), -1+val1, -1+val2, facecolor= 'white',fill=True,alpha=0.3))

ax.add_patch(Rectangle((val1, -val2), -1-val1, 1+val2, hatch='.',facecolor= 'white',fill=True,alpha=0.2,zorder=0))
ax.add_patch(Rectangle((val1, -val2), -val1, -1+val2, hatch='.',facecolor= 'white',fill=True,alpha=0.2,zorder=0))

ax.text(0.15, 0.9, 'Transcriptional and metabolic control',fontsize=15,backgroundcolor='white',c='black')
ax.text(0.35, -0.9, 'Metabolic control',fontsize=18,backgroundcolor='white',c='black')
ax.text(-0.9, -0.9, 'Transcriptional control',fontsize=18,backgroundcolor='white',c='black')
ax.text(-0.9, 0.9, 'Other',fontsize=18,backgroundcolor='white',c='black')

im=ax.scatter(df_cohen["RPSvsFBA "+name],
           df_cohen["RPSvsRAS "+name],100,df_cohen['FBAvsRAS '+name],zorder=1,cmap="viridis")

Texts=[]
for x,y,z in zip(df_cohen["RPSvsFBA "+name],df_cohen["RPSvsRAS "+name],df_cohen.index):
    if (x>=val1 and y>=-1) or (x<=-val1 and y<=-val2 and y>=-1):
        Texts.append(plt.text(x,y,z.split(':')[0],fontsize=15))
#adjust_text(Texts)#,arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
f.colorbar(im, ax=ax,orientation="horizontal")
im.set_clim(-1,1)

ax.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax.hlines(0,-1,1,colors="k")
ax.vlines(0,-1,1,colors="k")
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_xlabel("RPSvsFFD",fontsize=25)
ax.set_ylabel("RPSvsRAS",fontsize=25)
ax.grid()#linestyle = '--', linewidth = 0.5)
f.savefig(os.path.join(FIGUREDIR, scatterplotFigure))

# Heatmap mean values
df2_cohen=df_cohen[(df_cohen["RPSvsFBA "+name]>=0.2) | (df_cohen["RPSvsRAS "+name]>=0.2)]
index=df2_cohen.index
index=[ind.split(":")[0] for ind in index]
df_mean_RPS=df_met_mean2.loc[index]
df_FBA=pd.DataFrame()

    df_FBA[str(test[0])]=datasetFBA["median_"+test[0]]
    df_FBA[str(test[1])]=datasetFBA["median_"+test[1]]

df_RPS=pd.read_csv(os.path.join(OUTDIR, meansFile),index_col=0)
df_RPS=df_RPS.T

dfRPS=df_RPS.loc[index].T.div(df_RPS.loc[index].max(1)).T
dfFBA=df_FBA.loc[index].T.div(df_FBA.loc[index].max(1)).T

fig,ax=plt.subplots(1,2, figsize=(25, 15), gridspec_kw={'width_ratios': [1, 1.2],'wspace': 0.1})
a=sns.heatmap(dfRPS,cmap='viridis',ax=ax[0],xticklabels=list(dfFBA.columns.values),cbar=False)
b=sns.heatmap(dfFBA,cmap='viridis' ,ax=ax[1],yticklabels=[],xticklabels=list(dfFBA.columns.values),cbar=True)

a.set_xticklabels(a.get_xmajorticklabels(), fontsize = 18)
b.set_xticklabels(b.get_xmajorticklabels(), fontsize = 18)
a.set_yticklabels(a.get_ymajorticklabels(), fontsize = 18)
b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 18)

cbar = b.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)

ax[0].xaxis.tick_top()
ax[1].xaxis.tick_top()

fig.savefig(os.path.join(FIGUREDIR, heatmapMeansFigure),transparent=False,format="png")
