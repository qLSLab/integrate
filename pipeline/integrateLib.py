import pandas as pd
import numpy as np

def computeScore(subrule, out, s, scoreAlreadyComputed):
    if ' or ' in subrule and ' and ' not in subrule:
        if out.empty is False and len(scoreAlreadyComputed) == 0:
            score = out[s].sum()
        if out.empty is False and len(scoreAlreadyComputed) != 0:
            score = out[s].sum() + sum(scoreAlreadyComputed)
        elif out.empty is True and len(scoreAlreadyComputed) == 0:
            score = float("NaN")
        elif out.empty is True and len(scoreAlreadyComputed) != 0:
            score = sum(scoreAlreadyComputed)
    elif ' or ' not in subrule and ' and ' in subrule:
        if out.empty is False and len(scoreAlreadyComputed) != 0:
            score = min(list(out[s]) + scoreAlreadyComputed)
        elif out.empty is False and len(scoreAlreadyComputed) == 0:
            score = min(list(out[s]))
        elif out.empty is True and len(scoreAlreadyComputed) == 0:
            score = float("NaN")
        elif out.empty is True and len(scoreAlreadyComputed) != 0:
            score = min(scoreAlreadyComputed)
    return score

def differenceKeepingDuplicates(v1, v2):
    outputDifference = []
    while v1:
        # Get first item (and remove).
        item = v1.pop(0)
        if item in v2:
            v2.remove(item)
        else:
            outputDifference.append(item)
    return outputDifference

def reverse_dict(d):
    '''
    Aim: reverse a dictionary
    '''
    dinv = {}
    for k, v in d.items():
        for l in v:
            dinv[l]=k
    return dinv

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


def get_conc(MetMedia,S,S2,engro_met_dict):
    '''
    Aim: compute the concentration
    '''
    delta_fluxes_df=pd.DataFrame()
    data_quality_df=pd.DataFrame()

    for name in MetMedia.index.values: #for each cell line
        delta_fluxes_df[name],data_quality_df=getLineConc(MetMedia,S,S2,engro_met_dict,name)

    return delta_fluxes_df,data_quality_df

def splt2Net(df):
    '''
    Process random sampling output file by computing the net flux for reversible reactions
    '''
    lCols = df.columns.tolist()
    lCols_fb = []
    for col in lCols:
        if (col.endswith('_b') is True or col.endswith('_f') is True) and col[:-2] not in lCols_fb:
            lCols_fb.append(col[:-2])
    for col in lCols_fb:
        df[col] = df[col + '_f'] - df[col + '_b']
        df = df.drop(columns=[col + '_f', col + '_b'])
    print('dopo\n', df.shape, '\n')
    return df
