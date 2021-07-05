#!/usr/bin/env python
# coding: utf-8

# Load the Python libraries

# In[1]:


import cobra as cb
import numpy as np
import os
import sys
import ast
import matplotlib.pyplot 
import scipy
import itertools
import pandas as pd


# In[2]:


data_quality_filter=1                 #quality filter
valLog=np.log2(1.2)                   #filter on the media
namefile="outputs/resultsMetabolomic"  #name of the file where the data are written
metabolic_model='data\ENGRO2_rev.xml'
metabolic_data='data/metabolomics_LM.csv'
dict_to_convert_metnames='data/metsEngroVsMetabolomics.csv'
output_means="outputs\medie_Met.csv"          #means


# Load the stechiometric matrix of the metabolic model (Engro2)

# In[3]:


model = cb.io.read_sbml_model(metabolic_model)
S=cb.util.create_stoichiometric_matrix(model, array_type='DataFrame', dtype=None)
S


# Load the metabolic data (liquid spectrometry)

# In[4]:


metabolomicsLM=pd.read_csv(metabolic_data,sep=';',index_col=0)


# Create the dictionary to pass from the id reactions of Engro2 to Id reaction of metabolomic and viceversa

# In[5]:


#function to reverse a dictionary
def reverse_dict(d):
    dinv = {}
    for k, v in d.items():
        for l in v:
            dinv[l]=k
    return dinv


# In[6]:


#load conversion table between id metabolites in metamolomics data and id metabolites in ENGRO2 model and create dictionary
metIDConversion=pd.read_csv(dict_to_convert_metnames,sep=';')
metIDConversion['metId_engro'] = metIDConversion['metId_engro'].apply(ast.literal_eval)

#create dictionary Mass Spectometry to ENGRO
M_engro_dict=metIDConversion.set_index('metId_M')['metId_engro'].to_dict()

#create dictionary ENGRO to Mass Spectometry
engro_M_dict=reverse_dict(M_engro_dict)


# Create the column about the replicate

# In[7]:


#average by group
metabolomicsLM['Replicate']=metabolomicsLM.index.str.split('_').str[0]
MetMediaLM=metabolomicsLM


# Select the part of stechiometric matrix for the only the available metabolites

# In[8]:


S2=S.loc[[m for m in S.index.values if m in engro_M_dict.keys() and engro_M_dict[m] in MetMediaLM.columns]]
S2


# Create two function to compute the concentration

# In[9]:


def getLineConc(MetMedia,S,S2,engro_met_dict,pairName):
    
    #output dataframe
    delta_fluxes=pd.Series(index=S.columns,dtype =float)   
    data_quality_substracts=pd.Series(index=S.columns,dtype =float)
     
    for r in S2.columns: #for every reaction 
        
        n_substrates=sum(S.loc[:,r]<0);     #n° of reactions (c<0) in the stechiometric matrix (S)     
        quantified_substrates=len(S2[r][S[r]<0])

        columnS=-S2[r][S[r]<0]              #negative coefficient in the S2 matrix
        if quantified_substrates>0:
            fw_value=1; 
            for m in columnS.index.values: 
                m_concTot=MetMedia.loc[pairName,engro_met_dict[m]]
                fw_value=fw_value*(m_concTot**columnS[m]) 
        else:
            fw_value=0
        
        delta_fluxes[r]=fw_value

        
        if n_substrates>0:    
            data_quality_substracts[r]=quantified_substrates/n_substrates
        else:
            data_quality_substracts[r]=np.nan

    return delta_fluxes,data_quality_substracts


# In[10]:


def get_conc(MetMedia,S,S2,engro_met_dict):

    delta_fluxes_df=pd.DataFrame()
    data_quality_df=pd.DataFrame()

    for name in MetMedia.index.values: #for each cell line     
        delta_fluxes_df[name],data_quality_df=getLineConc(MetMedia,S,S2,engro_met_dict,name)

    return delta_fluxes_df,data_quality_df


# Compute concentration and quality of the data

# In[11]:


delta_fluxes_df,data_quality_df=get_conc(MetMediaLM,S,S2,engro_M_dict)


# How many reactions have all the metabolites?

# In[12]:


perc100=np.sum(data_quality_df==1)
print("perc100:",perc100)


# In[13]:


delta_fluxes_df


# # Statistical tests

# In[14]:


delta_fluxes_df_2=delta_fluxes_df.T
delta_fluxes_df_2["Replicate"]=delta_fluxes_df_2.index.str.split('_').str[0]
delta_fluxes_df_2=delta_fluxes_df_2.groupby(by="Replicate").agg(lambda x:np.mean(x))
delta_fluxes_df_2["Line"]=delta_fluxes_df_2.index.str.split(' ').str[0]
delta_fluxes_df_2


# Group by the replicates of each cell line

# In[15]:


delta_fluxes_df_3=delta_fluxes_df_2.groupby("Line").agg(lambda x:list(x))


# In[16]:


delta_fluxes_df_mean=delta_fluxes_df_2.groupby("Line").agg(lambda x:np.mean(x))
delta_fluxes_df_mean


# Save the file with all the means 

# In[17]:


delta_fluxes_df_mean.to_csv(output_means)


# In[18]:


delta_fluxes_df_3


# Rename indexes for compatability reasons with other files

# In[19]:


namesLine=['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']


# In[20]:


delta_fluxes_df_3.index=namesLine


# Create all the combinations to perform the statistical tests

# In[21]:


names=delta_fluxes_df_3.index.values
tests=itertools.combinations_with_replacement(names,2)
tests=[(test[0],test[1]) for test in tests if test[0]!=test[1]]


# In[22]:


tests


# In[23]:


test_df=pd.DataFrame(columns=delta_fluxes_df_3.columns[0:-1],index=tests)
mean_df=pd.DataFrame(columns=delta_fluxes_df_3.columns[0:-1],index=tests)


# ## Perform statistical tests  

# In[24]:


from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu


# In[25]:


import scipy.stats as stats
pval=0.05
eps=0.0001
names=delta_fluxes_df_3.index.values
for el in delta_fluxes_df_3.columns[:-1]:
    column=delta_fluxes_df_3[el]
    for test in tests:
        list1=np.round(column.loc[test[0]],4)
        list2=np.round(column.loc[test[1]],4)
        
        mean1=np.mean(list1)+eps
        mean2=np.mean(list2)+eps
        
        try:
            statistic_less,pvalue_less=ttest_ind(list1,list2,alternative= 'less')            
        except:
            pvalue_less=1

        try:
            statistic_greater,pvalue_greater=ttest_ind(list1,list2,alternative= 'greater')
        except:
            pvalue_greater=1
        #result of statistical tests   
        if pvalue_less<=pval:
             test_df[el][test]=-1
            
        elif pvalue_greater<=pval:
            test_df[el][test]=1
        else:
            test_df[el][test]=0
        
        mean_df[el][test]= np.log2(mean1/mean2)


# Set to 0 all the reactions where the differences between the means is low (<20% fold change)

# In[26]:


test_df2=test_df.copy()
test_df2[np.abs(mean_df)<=valLog]=0
test_df2


# Filter all the reactions with a low data quality

# In[27]:


#test_df3=test_df2.append(data_quality_df.to_frame().T)
#test_df3


# In[28]:


test_df3=test_df2.copy()
test_df3.loc["qc_substracts"]=data_quality_df
test_df3


# In[29]:


test_df4=test_df3.T[(test_df3.T["qc_substracts"]>=data_quality_filter)]
test_df4


# In[30]:


test_df4.to_csv(namefile+"_ttest.csv")


# In[31]:


mean_df2=mean_df.copy()
mean_df2.loc["qc_substracts"]=data_quality_df


# In[32]:


mean_df2=mean_df2.T[mean_df2.T["qc_substracts"]>=data_quality_filter]
mean_df2


# In[33]:


mean_df2.T.to_csv(namefile+"_mean.csv")

