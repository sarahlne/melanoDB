import pandas as pd
import numpy as np

def f(x):
    return(x*x)

def set_ge_df_value_col(samp, df_ge, list_genes):
    samp_tab = pd.DataFrame(index = list_genes)
    #print(samp)
    for gene in list_genes:
        if (np.shape(df_ge[(df_ge['sample_id']==samp) & (df_ge['HGNC']==gene)]) !=0):
            if (np.any(df_ge[(df_ge['sample_id']==samp) & (df_ge['HGNC']==gene)])):
                val = df_ge[(df_ge['sample_id']==samp) & (df_ge['HGNC']==gene)].value.to_numpy()[0]
                #print(val)
                samp_tab.loc[gene, samp]=val
    return(samp_tab)