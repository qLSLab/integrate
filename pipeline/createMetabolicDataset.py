#!/usr/bin/env python
# coding: utf-8

import cobra as cb
import numpy as np
import os
import ast
import itertools
import pandas as pd
import genericLib as gL
import integrateLib as iL
from scipy.stats import ttest_ind

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]

# setting input data
data_quality_filter = 1 # quality filter
valLog = np.log2(1.2)   # filter on the media
namefile = 'resultsMetabolomic'  # name of the file where the data are written
metabolic_model = 'ENGRO2_irrev.xml'
metabolic_data = 'metabolomics_LM.csv'
dict_to_convert_metnames = 'metsEngroVsMetabolomics.csv'
output_means = 'medie_Met.csv' # means


# Load the stoichiometric matrix of the input metabolic model
model = cb.io.read_sbml_model(os.path.join(MODELDIR, metabolic_model))
S = cb.util.create_stoichiometric_matrix(model, array_type='DataFrame', dtype=None)

# Load the metabolic data (Liquid chromatography–mass spectrometry)
metabolomicsLM = pd.read_csv(os.path.join(RAWDIR, metabolic_data), sep=';', index_col=0)

# Create the dictionary to pass from the id reactions of Engro2 to Id reaction of metabolomic and viceversa
# Load conversion table between ID of metabolites in metabolomics data and ID of metabolites in the input model and create dictionary

metIDConversion=pd.read_csv(os.path.join(RAWDIR, dict_to_convert_metnames), sep=';')
metIDConversion['metId_engro'] = metIDConversion['metId_engro'].apply(ast.literal_eval)

#create dictionary Mass Spectometry to ENGRO
M_engro_dict=metIDConversion.set_index('metId_M')['metId_engro'].to_dict()

#create dictionary ENGRO to Mass Spectometry
engro_M_dict=iL.reverse_dict(M_engro_dict)

# Create the column about the replicate

#average by group
metabolomicsLM['Replicate']=metabolomicsLM.index.str.split('_').str[0]
MetMediaLM=metabolomicsLM


# Select the part of stechiometric matrix for the only the available metabolites
S2=S.loc[[m for m in S.index.values if m in engro_M_dict.keys() and engro_M_dict[m] in MetMediaLM.columns]]

# Compute concentration and quality of the data
delta_fluxes_df,data_quality_df=iL.get_conc(MetMediaLM,S,S2,engro_M_dict)

# How many reactions have all the metabolites?

perc100=np.sum(data_quality_df==1)
print("How many reactions have all the metabolites?\t",perc100)

# Statistical tests
delta_fluxes_df_2=delta_fluxes_df.T
delta_fluxes_df_2["Replicate"]=delta_fluxes_df_2.index.str.split('_').str[0]
delta_fluxes_df_2=delta_fluxes_df_2.groupby(by="Replicate").agg(lambda x:np.mean(x))
delta_fluxes_df_2["Line"]=delta_fluxes_df_2.index.str.split(' ').str[0]
delta_fluxes_df_2

# Group by the replicates of each cell line
delta_fluxes_df_3=delta_fluxes_df_2.groupby("Line").agg(lambda x:list(x))
delta_fluxes_df_mean=delta_fluxes_df_2.groupby("Line").agg(lambda x:np.mean(x))

# Save the file with all the means
delta_fluxes_df_mean.to_csv(os.path.join(OUTDIR, output_means))

# Rename indexes for compatability reasons with other files
namesLine=['MCF102A','MDAMB231','MDAMB361','MCF7','SKBR3']
delta_fluxes_df_3.index=namesLine

# Create all the combinations to perform the statistical tests
names=delta_fluxes_df_3.index.values
tests=itertools.combinations_with_replacement(names,2)
tests=[(test[0],test[1]) for test in tests if test[0]!=test[1]]

test_df=pd.DataFrame(columns=delta_fluxes_df_3.columns[0:-1],index=tests)
mean_df=pd.DataFrame(columns=delta_fluxes_df_3.columns[0:-1],index=tests)


# Perform statistical tests
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
test_df2=test_df.copy()
test_df2[np.abs(mean_df)<=valLog]=0

# Filter all the reactions with a low data quality
test_df3=test_df2.copy()
test_df3.loc["qc_substracts"]=data_quality_df
test_df4=test_df3.T[(test_df3.T["qc_substracts"]>=data_quality_filter)]
test_df4.to_csv(os.path.join(OUTDIR,  namefile + '_ttest.csv'))

mean_df2=mean_df.copy()
mean_df2.loc["qc_substracts"]=data_quality_df
mean_df2=mean_df2.T[mean_df2.T["qc_substracts"]>=data_quality_filter]
mean_df2.T.to_csv(os.path.join(OUTDIR,  namefile + '_mean.csv'))
