import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json
import math

from math import floor


############################################################################################
#
#                                        LOUVEAU class
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


class Louveau():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path):
        """
        Read the csv file extract from the clinical study (DATA_RESISTANCE_BITHERAPIE_extract_051022.csv) and returns a list of dictionnaries
        where terms represent patients filled with their informations (age, sex, stage, etc...). 

        :returns: list of all patients in the study
        :rtype: list
        """
        list_patients = []

        # import table
        table = pd.read_excel(file_path, header=0)
        # exclude patients with no inhibitors as first or second line therapy ('ligne_bitherap1'!=1 | 'ligne_bitherap1'!=2)
        table = table[(table['ligne_bitherap1']==1) | (table['ligne_bitherap1']==2)].reset_index(drop=True)
        # change values of some table columns 
        table_sex = table['Sexe']
        table_sex = table_sex.replace(['M','F'],['male', 'female'])

        table_age = [floor(x/365) for x in table['bitherap1_age_init_jours']]

        table_stage = [re.split(' ', x)[0] for x in table['classif_j1biotherap1']] 
        table_M_stage = [re.split(' ', x)[1] for x in table['classif_j1biotherap1']]

        table_LDH = []
        for elem in table['ldh_j1biotherap1']:
            if(elem > 225):
                table_LDH.append('elevated')
            else:
                if (pd.isna(elem)):
                    table_LDH.append(elem)
                else:
                    table_LDH.append('normal')
        
        table_OS_statut = table['deces'].replace([0,1],['alive','dead'])
        table_OS_month = [round((x-y)/30, 2) for x, y in zip(table['age_derniere_nouvelle_jours'], table['bitherap1_age_init_jours'])]
        # Fetch PFS month
        list_PFS_month = []
        list_PFS_statut = []
        for i in range(len(table)):
            if(np.isnan(table.at[i,'bitherap1_age_prog_jours'])==True):
                PFS_month = round((table['bitherap1_age_arret_jours'][i]-table['bitherap1_age_init_jours'][i])/30.5,1)
                PFS_statut = 0.0 # no prog; reason 1 toxicity, reason 2 unknow stop; reason 3 the study stopped
            else:
                PFS_month = round((table['bitherap1_age_prog_jours'][i]-table['bitherap1_age_init_jours'][i])/30.5,1)
                PFS_statut = table['bitherap1_prog'][i]
            list_PFS_month.append(PFS_month)
            list_PFS_statut.append(PFS_statut)
        
        table_BRAF_mut = [re.sub('BRAF', '', x) for x in table['Statut_BRAFJ1']]
        table_DCR = table['bitherap1_meilleure_reponse'].replace(['RP','RC'],['PR','CR'])
        table_drug = table['nom_bitherap1'].replace(['vemu+cobi', 'dabra+trame'], ['vemurafenib + cobimetinib', 'dabrafenib + trametinib'])
        table_brain_met = table['Metastases_cerebrales'].replace([0,1],['no', 'yes'])

        list_prior_treatment = ['yes' if x==2 else 'no' for x in table['ligne_bitherap1']]

        pat_prefix = 'LM_'
        sample_prefix= 'LMSAM_'

        #create dictionaries
        for ind in table.index:
            patient_dict=dict(
                patient_ID = pat_prefix+str(table['Numero'][ind]),
                original_patientID = str(table['Numero'][ind]),
                internal_patientID=pat_prefix+str(table['Numero'][ind]),
                sex = table_sex[ind],
                age = table_age[ind],
                stage = table_stage[ind],
                M_stage = table_M_stage[ind],
                LDH = table_LDH[ind],
                os_statut = table_OS_statut[ind],
                os_months = table_OS_month[ind],
                pfs_statut = str(int(list_PFS_statut[ind])),
                pfs = list_PFS_month[ind],
                braf_mut = table_BRAF_mut[ind],
                disease_control_rate = table_DCR[ind],
                prelevement_temporality = np.NaN,
                drug = table_drug[ind],
                brain_metastasis = table_brain_met[ind],
                immunotherapy_treatment = 'no',
                prior_mapk_treatment = list_prior_treatment[ind],
                CNA='no',
                SNV='yes',
                GEX='yes',
                source = dict(
                    title = 'Baseline Genomic Features in BRAFV600-Mutated Metastatic Melanoma Patients Treated with BRAF Inhibitor + MEK Inhibitor in Routine Care',
                    author =  'Baptiste Louveau, Samia Mourah',
                    journal =  'Cancers (Basel)',
                    location = 'Paris (France)',
                    date = 2019) 
            )
            
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)

            list_patients.append(json.loads(convert_dict))
        print('the list of "Louveau" patients has been created')
        return(list_patients)
    
    def parse_mutations_from_file(self, file_mutation_path):
        list_louveau_mutations = []
        table_mutations = pd.read_excel(file_mutation_path)
        table_mutations = table_mutations[(table_mutations['ligne_bitherap1']==1) | (table_mutations['ligne_bitherap1']==2)].reset_index(drop=True)

        pat_prefix = 'LM_'
        sample_prefix= 'LMSAM_'

        for ind in table_mutations.index:
            list_mutations_pat = table_mutations['Analyse NGS_cur'][ind]
            if isinstance(list_mutations_pat, str):
                list_mutations_pat = table_mutations['Analyse NGS_cur'][ind].split('\n')
                for elem in list_mutations_pat:
                    t = re.search('(?<=p.)(\w*)',elem)
                    if not isinstance(t, type(None)):
                        hgvsp_short = t.group(1)
                    else:
                        hgvsp_short = np.NaN
                    snp_dict = dict(
                        patientID = pat_prefix+str(table_mutations['Numero'][ind]),
                        sample_id = sample_prefix+str(table_mutations['Numero'][ind]),
                        HGNC = elem.split(':')[0].replace(' ', ''),
                        HGVSp_short = hgvsp_short,
                        mutated = 'yes',
                        temporality = 'pre treatment',
                        source = dict(
                            title = 'Baseline Genomic Features in BRAFV600-Mutated Metastatic Melanoma Patients Treated with BRAF Inhibitor + MEK Inhibitor in Routine Care',
                            author =  'Baptiste Louveau, Samia Mourah',
                            journal =  'Cancers (Basel)',
                            location = 'Paris (France)',
                            date = 2019)
                    )
                    for (key, value) in snp_dict.items():
                        if (pd.isna(value)):
                            snp_dict[key]=None
                        if (isinstance(value, str)):
                            if(value=='NAN'):
                                snp_dict[key]=None

                    convert_dict = json.dumps(snp_dict, cls=NpEncoder)
                    list_louveau_mutations.append(json.loads(convert_dict))
        print('the list of "Louveau" mutations has been created')
        return(list_louveau_mutations)
    
    def parse_gene_expression_from_file(self, file_gene_expr):
        list_louveau_gene_expr = []
        table_gene_expr = pd.read_csv(file_gene_expr, sep=';', encoding='latin1')
        table_gene_expr = table_gene_expr[table_gene_expr['ligne_bitherap1']==1].reset_index(drop=True)

        pat_prefix = 'LM_'
        sample_prefix= 'LMSAM_'

        for i in range(0,len(table_gene_expr)):
            pat_ID = pat_prefix+str(table_gene_expr['Numero'][i])
            samp_ID = sample_prefix+str(table_gene_expr['Numero'][i])

            for j in range(24, len(table_gene_expr.columns)):
                if (len(table_gene_expr.iloc[i,j])!=0):
                        val=float(table_gene_expr.iloc[i,j].replace(",","."))
                else:
                    val=np.NAN
                gene_expr_dict= dict(
                    patientID = pat_ID,
                    sample_id = samp_ID,
                    HGNC = table_gene_expr.columns[j],
                    GeneID=np.NaN,
                    description = np.NaN,
                    value = val, 
                    temporality = 'pre treatment',
                    source = dict(
                            title = 'Baseline Genomic Features in BRAFV600-Mutated Metastatic Melanoma Patients Treated with BRAF Inhibitor + MEK Inhibitor in Routine Care',
                            author =  'Baptiste Louveau, Samia Mourah',
                            journal =  'Cancers (Basel)',
                            location = 'Paris (France)',
                            date = 2019)
                )
                for (key, value) in gene_expr_dict.items():
                        if (pd.isna(value)):
                            gene_expr_dict[key]=None
                        if (isinstance(value, str)):
                            if(value=='NAN'):
                                gene_expr_dict[key]=None

                convert_dict = json.dumps(gene_expr_dict, cls=NpEncoder)
                list_louveau_gene_expr.append(json.loads(convert_dict))
        print('the list of "Louveau" gene expressions has been created')
        return(list_louveau_gene_expr)



