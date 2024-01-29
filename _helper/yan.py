import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json


############################################################################################
#
#                                        Yan class
#                                         
############################################################################################

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.nan):
            return None
        else:
            return super(NpEncoder, self).default(obj)

class Yan():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path):
        """
        Read the csv file extract from the clinical study (Fichier_melanome_RL.xlsx) and returns a list of dictionnaries
        where terms represent proteins filled with their informations (age, sex, stage, etc...). 

        :returns: list of all patients in the study
        :rtype: list
        """
        list_patients = []

        # import table
        table = pd.read_csv(file_path, header=0)
        
        # change values of some table columns 
        table_sex = table['Sex']
        table_sex = table_sex.replace(['Male','Female'],['male', 'female'])

        table_OS = table['OS Censor'].replace([0,1],[1,0])
        table_OS = table_OS.replace([0,1], ['alive', 'dead'])
        table_PFS = table['PFS (Months)']
        table_PFS_Censor = table['PFS Censor'].replace([0,1],[1,0])
        table_LDH = table['LDH'].replace(['Normal','Elevated'], ['normal', 'elevated'])

        table_AJCC_Stage = table['Stage'].replace(['M1A','M1B', 'M1C', 'Unresectable', 'IIIC', 'Unresectable Stage IIIC'], ['IV', 'IV', 'IV', 'IV', 'III', 'III'])

        # create dictionnaries
        for ind in table.index:
            patient_dict=dict(
                patient_ID = table['Patient ID'][ind],
                original_patientID = table['Patient ID'][ind],
                internal_patientID="YR_"+f"{table['Patient ID'][ind]}",
                sex = table_sex[ind],
                age = table['Age'][ind],
                stage = table_AJCC_Stage[ind],
                M_stage = table['Stage'][ind],
                LDH = table_LDH[ind],
                os_statut = table_OS[ind],
                os_months = table['OS (Months)'][ind],
                pfs_statut = str(int(table_PFS_Censor[ind])),
                pfs = table_PFS[ind],
                braf_mut = table['BRAF V600 Mut'][ind],
                disease_control_rate = table['BORR'][ind],
                drug = table['Rx'][ind],
                brain_metastasis = np.NaN,
                immunotherapy_treatment = np.NaN,
                CNA='no',
                SNV='no',
                GEX='yes',
                source = dict(
                    title = 'Genomic Features of Exceptional Response in Vemurafenib + Cobimetinib–treated Patients with BRAFV600-mutated Metastatic Melanoma',
                    author =  'Yibing Yan, Antoni Ribas',
                    journal =  'Clinical Cancer Research',
                    location = 'San Francisco, California (United States)',
                    date = 2019)  
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
            
            convert_dict = json.dumps(patient_dict, cls=NpEncoder)

            list_patients.append(json.loads(convert_dict))

        print('the list of "Yan" patients has been created')
        return(list_patients)
    
    def parse_gene_expression_from_file(self, file_gene_expr, file_ids_match, file_geneids_HGNC):
        list_yan_gene_expr = []
        table_gene_expr = pd.read_csv(file_gene_expr, sep=';', encoding='latin1')
        table_gene_expr.set_index('Unnamed: 0', inplace=True)
        df_ids = pd.read_csv(file_ids_match)
        df_ids = df_ids.astype({'Patient ID': 'str'})
        df_ids['Patient ID'] = [df_ids['Patient ID'][i].replace('.0', '') for i in range(0, len(df_ids))]
        dict_ids = {elem[0]:elem[1] for elem in df_ids.to_dict(orient='tight')['data']}
        dict_ids = {v: k for k, v in dict_ids.items()}
        df_hgnc=pd.read_csv(file_geneids_HGNC)
        df_hgnc = df_hgnc.astype({'query':'str'})
        df_hgnc['query'] = ['GeneID:'+df_hgnc['query'][i] for i in range(0,len(df_hgnc['query']))]
        dict_hgnc={elem[0]:elem[1] for elem in df_hgnc.to_dict(orient='tight')['data']}

        for j in range(0, len(table_gene_expr.columns)):
            if (table_gene_expr.columns[j] in dict_ids.keys()):
                pat_ID = "YR_"+dict_ids[table_gene_expr.columns[j]]
            else:
                pat_ID=np.NAN
            samp_ID = table_gene_expr.columns[j]

            for i in range(0, len(table_gene_expr)):
                val=table_gene_expr.iloc[i,j]
                if (table_gene_expr.index[i] in dict_hgnc.keys()):
                    hgnc=dict_hgnc[table_gene_expr.index[i]]
                else:
                    hgnc=np.NAN

                gene_expr_dict= dict(
                    patientID = pat_ID,
                    sample_id = samp_ID,
                    HGNC = hgnc,
                    GeneID= table_gene_expr.index[i],
                    description = np.NaN,
                    value = val, 
                    temporality = 'pre treatment',
                    source = dict(
                            title = 'Genomic Features of Exceptional Response in Vemurafenib + Cobimetinib–treated Patients with BRAFV600-mutated Metastatic Melanoma',
                            author =  'Yibing Yan, Antoni Ribas',
                            journal =  'Clinical Cancer Research',
                            location = 'San Francisco, California (United States)',
                            date = 2019)
                )
                for (key, value) in gene_expr_dict.items():
                        if (pd.isna(value)):
                            gene_expr_dict[key]=None
                        if (isinstance(value, str)):
                            if(value=='NAN'):
                                gene_expr_dict[key]=None
                convert_dict = json.dumps(gene_expr_dict, cls=NpEncoder)
                list_yan_gene_expr.append(json.loads(convert_dict))
        print('the list of "Yan" gene expressions has been created')
        return(list_yan_gene_expr)


