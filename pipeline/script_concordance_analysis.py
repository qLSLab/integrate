#!/usr/bin/env python
# coding: utf-8

# # Concordance analysis

# Load the libraries

# In[ ]:


import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import cobra as cb


# In[ ]:


val=np.log2(1.2)     #threshold
weight="linear"      #weight on the Kappa cohen
name_file="outputs/resultsMetabolomic"  #directory of the data
dir_output="outputs/"
model_name="models/ENGRO2_irrev.xml"
# Set the name of cell lines

# In[ ]:


namesLine=['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']


# Set all the combinations for the tests

# In[ ]:


tests=itertools.combinations_with_replacement(namesLine,2)
tests=[(test[0],test[1]) for test in tests if test[0]!=test[1]]
tests


# ## Load and create the concordance datasets

# Load the metabolitc model

# In[ ]:


model = cb.io.read_sbml_model(model_name)
keys=[reaction.id for reaction in model.reactions]
formule=[reaction.reaction for reaction in model.reactions]
dict_formule={key:value for key,value in zip(keys,formule)}


# In[ ]:


model.reactions.get_by_id("Biomass").objective_coefficient = 0


# In[ ]:


GPR_reactions=[reaction.id for reaction in model.reactions if reaction.gene_reaction_rule!='']


# Load RAS and FBA data

# In[ ]:


datasetsFBA=list()
datasetsRAS=list()
for test in tests:
    datasetsFBA.append(pd.read_csv('rawData/concordance_data/mwuTest_'+test[0]+'_vs_'+test[1]+'.csv',
                                sep="\t",index_col=0).sort_index())
    datasetsRAS.append(pd.read_csv('rawData/concordance_data/ras_'+test[0]+'_vs_'+test[1]+'.csv',
                                sep="\t",index_col=0).loc[GPR_reactions])


# Create the dataset of concordance and median/mean for FBA and RAS

# In[ ]:


df_concFBA=pd.DataFrame(index=datasetsFBA[0].index)
df_medianFBA=pd.DataFrame(index=datasetsFBA[0].index)

df_concRAS=pd.DataFrame(index=datasetsRAS[0].index)
df_meanRAS=pd.DataFrame(index=datasetsRAS[0].index)


# Perform statistical tests

# In[ ]:


eps=0.00001
for datasetFBA,datasetRAS,test in zip(datasetsFBA,datasetsRAS,tests):
    #RAS
    df_meanRAS[str(test)]=np.log2((datasetRAS["mean_"+test[0]]+eps)/(datasetRAS["mean_"+test[1]]+eps))#
    df_concRAS[str(test)]=datasetRAS["result"]
    #FBA
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

# In[ ]:


df_met=pd.read_csv(name_file+"_ttest.csv",index_col=0)
df_met_mean=pd.read_csv(name_file+"_mean.csv",index_col=0).T


# Set to 0 the mean below specific threshold

# In[ ]:


df_met[np.abs(df_met_mean)<val]=0
df_concFBA[np.abs(df_medianFBA)<val]=0
df_concRAS[np.abs(df_concRAS)<val]=0


# Create index of Common reactions among the concordance datasets

# In[ ]:


index_commonMETvsFBA=df_met.index.intersection(df_concFBA.index)
index_commonMETvsRAS=df_met.index.intersection(df_concRAS.index)
index_commonFBAvsRAS=df_concFBA.index.intersection(df_concRAS.index)


# In[ ]:


indexCommon=index_commonMETvsFBA.intersection(index_commonMETvsRAS)


# Create the indexes of all the reactions in the common datasets

# In[ ]:


indexUnion=df_met.index.union(df_concRAS.index).union(df_concFBA.index)


# Sort the columns for comparison

# In[ ]:


df_concRAS2=df_concRAS[np.sort(df_concRAS.columns)]
df_meanRAS2=df_meanRAS[np.sort(df_meanRAS.columns)]


# In[ ]:


df_concFBA2=df_concFBA[np.sort(df_concFBA.columns)]
df_medianFBA2=df_medianFBA[np.sort(df_medianFBA.columns)]


# In[ ]:


df_met2=df_met[np.sort(df_met.columns[:-1])]
df_met_mean2=df_met_mean[np.sort(df_met_mean.columns[:-1])]


# ## Anaysis of a single reaction

# Analysis of a specific reaction

# In[ ]:


df_concRAS2.loc["ACONT"]


# In[ ]:


reaction="ACONT"#
print(df_medianFBA2.loc[reaction])
print(df_met_mean2.loc[reaction])
print(stats.pearsonr(df_medianFBA2.loc[reaction],df_met_mean2.loc[reaction]))
fig=plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size':20})
plt.scatter(df_medianFBA2.loc[reaction],df_met_mean2.loc[reaction])
plt.grid()
plt.xlabel("FFD fold change")
plt.ylabel("RPS fold change")
plt.xlim([-3.5,3.5])
plt.ylim([-3.5,3.5])
plt.title("ACONT")

fig.savefig(dir_output+"ACONT_reaction.png")


# ## Statistics tests

# In[ ]:


import sklearn.metrics
kappa_values=pd.DataFrame(index=indexUnion,columns=['LMvsFBA cohen', 'LMvsRAS cohen','FBAvsRAS cohen',
                                                       'LMvsFBA pearson','LMvsRAS pearson','FBAvsRAS pearson'])
for r in kappa_values.index:
    compare=pd.DataFrame()
    if r in df_met2.index and r in df_concFBA2.index:
        compare['MET']=df_met2.loc[r]
        compare['FBA']=df_concFBA2.loc[r]
        kappa_values.loc[r,'LMvsFBA cohen']=sklearn.metrics.cohen_kappa_score(compare['MET'], compare['FBA'],
                                                                           labels=[-1,0,1],weights=weight, sample_weight=None)
        kappa_values.loc[r,'LMvsFBA pearson'],p=stats.pearsonr(df_met_mean2.loc[r],df_medianFBA2.loc[r])
    if r in df_met2.index and r in df_concRAS2.index:
        compare['MET']=df_met2.loc[r]
        compare['RAS']=df_concRAS2.loc[r]                                              
        kappa_values.loc[r,'LMvsRAS cohen']=sklearn.metrics.cohen_kappa_score(compare['MET'], compare['RAS'],
                                                                           labels=[-1,0,1],weights=weight, sample_weight=None)
        kappa_values.loc[r,'LMvsRAS pearson'],p=stats.pearsonr(df_met_mean2.loc[r],df_meanRAS2.loc[r])
    if r in df_concFBA2.index and r in df_concRAS2.index:
        compare['FBA']=df_concFBA2.loc[r]    
        compare['RAS']=df_concRAS2.loc[r]                                           
        kappa_values.loc[r,'FBAvsRAS cohen']=sklearn.metrics.cohen_kappa_score(compare['FBA'],compare['RAS'],
                                               labels=[-1,0,1],weights=weight, sample_weight=None)
        kappa_values.loc[r,'FBAvsRAS pearson'],p=stats.pearsonr(df_medianFBA2.loc[r],df_meanRAS2.loc[r])


# In[ ]:


kappa_values=kappa_values.sort_values(by='LMvsRAS cohen',ascending=False)
kappa_values.index=[str(el)+':'+str(dict_formule[el]) for el in kappa_values.index]
kappa_values


# In[ ]:


kappa_values.to_csv(dir_output+"df_concordance.csv")


# In[ ]:


kappa_values.fillna(-3,inplace=True)  ## just the create a mask map for the next figure


# ## Figures of the paper

# ### Heatmap Cohen

# In[ ]:


from matplotlib.colors import ListedColormap

fig, ax=plt.subplots(1, 2, figsize=(16, 30), gridspec_kw={'width_ratios': [3, 0.7],'wspace': 0})

df=kappa_values.loc[((kappa_values['LMvsFBA cohen']>0.2) | (kappa_values['LMvsRAS cohen']>0.2) ),:].sort_values('LMvsFBA cohen',ascending=False)
df_cohen=df.loc[:,['LMvsFBA cohen', 'LMvsRAS cohen']]

mask = df_cohen.values==-3

a=sns.heatmap(df_cohen,annot=True,mask=mask, xticklabels=["RPSvsFFD","RPSvsRAS"],
            cbar=True,
            cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.8),
            yticklabels=[],cmap="viridis",
            ax=ax[1],edgecolor="k",linewidths=0.1,annot_kws={"fontsize":18})

a.set_xticklabels(a.get_xmajorticklabels(), fontsize = 15)

df_new=df
df_new["diff"]=[1 if np.sign(el)==np.sign(el2)  else 0 for el,el2 in zip(df_cohen['LMvsFBA cohen'],df_cohen['LMvsRAS cohen'])]


ax[1].xaxis.tick_top()

cmap=ListedColormap(["salmon","grey"])

# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))

# Set alpha
my_cmap[0,-1] = 0.5
my_cmap[1,-1] = 0.3

# Create new colormap
my_cmap = ListedColormap(my_cmap)


g=sns.heatmap(np.array(df_new["diff"].values).reshape(29,1),ax=ax[0],
             cbar=True,
             cbar_kws = dict(use_gridspec=False,location="bottom", shrink=0.8),              
              cmap=my_cmap,
             vmax=1,vmin=0,annot = np.array(df_new["diff"].index).reshape(29,1),fmt = '',
             xticklabels=["Reactions"],yticklabels=[],edgecolor="k",linewidths=0.1,annot_kws={"color":"black","fontsize":20})

g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15)


ax[0].xaxis.tick_top()


# In[ ]:


fig.savefig(dir_output+"concordanceHeatmap.png")


# ## Heatmap of all reactions

# In[ ]:


from matplotlib.colors import ListedColormap

fig, ax=plt.subplots(1, 2, figsize=(20, 50), gridspec_kw={'width_ratios': [3, 0.7],'wspace': 0})

df=kappa_values.loc[((kappa_values['LMvsFBA cohen']!=-3) | (kappa_values['LMvsRAS cohen']!=-3) ),:].sort_values('LMvsFBA cohen',ascending=False)
df_cohen=df.loc[:,['LMvsFBA cohen', 'LMvsRAS cohen']]

mask = df_cohen.values==-3

a=sns.heatmap(df_cohen,annot=True,mask=mask, xticklabels=["RPSvsFFD","RPSvsRAS"],
            cbar=True,
            cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.8),
            yticklabels=[],cmap="viridis",
            ax=ax[1],edgecolor="k",linewidths=0.1,annot_kws={"fontsize":18})

a.set_xticklabels(a.get_xmajorticklabels(), fontsize = 15)

df_new=df
df_new["diff"]=[1 if np.sign(el)==np.sign(el2)  else 0 for el,el2 in zip(df_cohen['LMvsFBA cohen'],df_cohen['LMvsRAS cohen'])]


ax[1].xaxis.tick_top()

cmap=ListedColormap(["salmon","grey"])

# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))

# Set alpha
my_cmap[0,-1] = 0.5
my_cmap[1,-1] = 0.3

# Create new colormap
my_cmap = ListedColormap(my_cmap)


g=sns.heatmap(np.array(df_new["diff"].values).reshape(81,1),ax=ax[0],
             cbar=True,
             cbar_kws = dict(use_gridspec=False,location="bottom", shrink=0.8),              
              cmap=my_cmap,
             vmax=1,vmin=0,annot = np.array(df_new["diff"].index).reshape(81,1),fmt = '',
             xticklabels=["Reactions"],yticklabels=[],edgecolor="k",linewidths=0.1,annot_kws={"color":"black","fontsize":20})

g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15)


ax[0].xaxis.tick_top()


# In[ ]:


fig.savefig(dir_output+"concordanceHeatmap_total.png")


# ### Scatter plot reactions

# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib.patches import Rectangle

df=kappa_values.sort_values('LMvsFBA cohen',ascending=False)
df_cohen=df.loc[:,['LMvsFBA cohen', 'LMvsRAS cohen','FBAvsRAS cohen']]


f = plt.figure(figsize=(13,23))
plt.rcParams.update({'font.size':20})
ax = f.add_subplot(111)



ax.add_patch(Rectangle((0, 0), 1, -1, facecolor= 'salmon',fill=True,alpha=0.5,zorder=0))
ax.add_patch(Rectangle((0, 0), 1, 1, facecolor= 'grey',fill=True,alpha=0.3))

#ax.add_patch(Rectangle((-1, -1), 1, 1, facecolor= 'lightpink',fill=True,alpha=0.5))

ax.text(0.05, 0.9, 'Transcriptional and metabolic control',fontsize=18,backgroundcolor='white',c='black')
ax.text(0.05, -0.9, 'Metabolic control',fontsize=18,backgroundcolor='white',c='black')
ax.text(-0.9, -0.9, 'Transcriptional control',fontsize=18,backgroundcolor='white',c='black')

im=ax.scatter(df_cohen["LMvsFBA cohen"],
           df_cohen["LMvsRAS cohen"],100,df_cohen['FBAvsRAS cohen'],zorder=1,cmap="viridis")

Texts=[]
for x,y,z in zip(df_cohen["LMvsFBA cohen"],df_cohen["LMvsRAS cohen"],df_cohen.index):
    if (x>=-1 and x<=-0.25 and y>=-1) or (x>=0.25 and y>=-1):
        Texts.append(plt.text(x,y,z.split(':')[0],fontsize=15))
adjust_text(Texts)#,arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
f.colorbar(im, ax=ax,orientation="horizontal")
im.set_clim(-1,1)

ax.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax.hlines(0,-1,1,colors="k")
ax.vlines(0,-1,1,colors="k")
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_xlabel("RPS vs FFD cohen",fontsize=25)
ax.set_ylabel("RPS vs RAS cohen",fontsize=25)
ax.grid(linestyle = '--', linewidth = 0.5)


# In[ ]:


f.savefig(dir_output+"scatterplot.png")


# ## Heatmap mean values

# In[ ]:


df2_cohen=df_cohen[(df_cohen["LMvsFBA cohen"]>=0.2) | (df_cohen["LMvsRAS cohen"]>=0.2)]


# In[ ]:


index=df2_cohen.index
index=[ind.split(":")[0] for ind in index]
df_mean_RPS=df_met_mean2.loc[index]


# In[ ]:


df_FBA=pd.DataFrame()


# In[ ]:


eps=0.00001
for datasetFBA,datasetRAS,test in zip(datasetsFBA,datasetsRAS,tests):
    df_FBA[str(test[0])]=datasetFBA["median_"+test[0]]
    df_FBA[str(test[1])]=datasetFBA["median_"+test[1]]


# In[ ]:


df_RPS=pd.read_csv("data\medie_Met.csv",index_col=0)
df_RPS=df_RPS.T
df_RPS


# In[ ]:


dfRPS=df_RPS.loc[index].T.div(df_RPS.loc[index].max(1)).T
dfFBA=df_FBA.loc[index].T.div(df_FBA.loc[index].max(1)).T


# In[ ]:


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


# In[ ]:


fig.savefig(dir_output+"heatmap_means.png",transparent=False,format="png")


# In[ ]:




