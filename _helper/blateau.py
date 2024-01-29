import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json


############################################################################################
#
#                                        Blateau class
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


class Blateau():
    
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path, cond):
        """
        Read the csv file extract from the clinical study (Fichier_melanome_RL.xlsx) and returns a list of dictionnaries
        where terms represent proteins filled with their informations (age, sex, stage, etc...). 

        :returns: list of all patients in the study
        :rtype: list
        """

        list_patients = []
        list_snps = []

        # import table
        table = pd.read_excel(file_path, header=0)
        table = table.drop(table.index[53:89],0, inplace=False)
        
        # change values of some table columns
        
        table_sexe = table['sexe']
        table_sexe = table_sexe.replace([0,1], ['male', 'female'])

        table_AJCC = table["AJCC à l'initiation du traitement\n0=I; 1=II; 2=III; 3=IV"]
        table_AJCC = table_AJCC.replace([0,1,2,3],['I','II','III','IV'])

        table_LDH = table['LDH (normales = 0, augmentées = 1)']
        table_LDH = table_LDH.replace([0,1], ['normal', 'elevated'])

        table_OS = table['Statut OS (0 = neg; 1 = death)']
        table_OS = table_OS.replace([0,1], ['alive', 'dead'])

        table_PFS_statut = table['Progression sous traitement (0=neg; 1 = event)'].convert_dtypes(convert_floating=False)

        table_disease_control_rate = table['Statut (0= stable disease, 1= progression, 2= partial response, 3= complete response)']
        table_disease_control_rate = table_disease_control_rate.replace([0,1,2,3],['SD', 'PD', 'PR', 'CR'])

        table_braf_mut = table['Mut BRAF (0=V600E, 1=V600K)']
        table_braf_mut = table_braf_mut.replace([0,1],['V600E', 'V600K'])

        table_brain_met = table['Méta cérébrales (0=neg, 1= event)']
        table_brain_met = table_brain_met.replace([0,1],['no', 'yes'])

        table_immuno_T = table['Immunothérapie (0=neg, 1=event)']
        table_immuno_T = table_immuno_T.replace([0,1],['no', 'yes'])

        table_prior = table['Anti BRAF 1ère ligne (0=oui; 1=non)'].replace({0.0: 'no', 1.0:'yes'})

        list_dabrafenib = table.index[table['DABRAFENIB']==True].to_list()
        list_vemurafenib = table.index[table['VEMURAFENIB']==True].to_list()
        list_dabrafenib_trametinib = table.index[table['DABRAFENIB + TRAMETINIB']==True].to_list()
        list_vemurafenib_cobimetinib = table.index[table['VEMURAFENIB + COBIMETINIB']==True].to_list()
        
        table_agent = pd.Series(['nan' for i in range(len(table))])
        table_agent[list_dabrafenib]='dabrafenib'
        table_agent[list_vemurafenib]='vemurafenib'
        table_agent[list_dabrafenib_trametinib]='dabrafenib + trametinib'
        table_agent[list_vemurafenib_cobimetinib] = 'vemurafenib + cobimetinib'
        "Clark level0=II; 1=III; 2=IV; 3=V"

        pat_prefix = 'BS'
        sample_prefix= 'BSAM'

        # create dictionnaries for patients and their snps

        # PATIENTS
        for ind in table.index:
            patient_dict=dict(
                patient_ID=pat_prefix+"_"+f"{ind:03}",
                original_patientID = f"{ind:01}",
                internal_patientID=pat_prefix+"_"+f"{ind:03}",
                sex = table_sexe[ind],
                age = table['age'][ind],
                stage = table_AJCC[ind],
                LDH = table_LDH[ind],
                os_statut = table_OS[ind],
                os_months = (table['OS'][ind])/30,
                pfs_statut = str(table_PFS_statut[ind]),
                pfs = table['PFS en mois'][ind],
                braf_mut = table_braf_mut[ind],
                disease_control_rate = table_disease_control_rate[ind],
                drug = table_agent[ind],
                brain_metastasis = table_brain_met[ind],
                immunotherapy_treatment = table_immuno_T[ind],
                prior_mapk_treatment = table_prior[ind],
                CNA='no',
                SNV='yes',
                GEX='no',
                source = dict(
                    title = 'TERT Promoter Mutation as an Independent Prognostic Marker for Poor Prognosis MAPK Inhibitors-Treated Melanoma',
                    author =  'Pauline Blateau, Jerome Solassol',
                    journal =  'Cancers',
                    location = 'Montpellier (France)',
                    date = 2020)
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)
            list_patients.append(json.loads(convert_dict))

        # SNVS
            braf_stat=table.iloc[ind,39]
            if braf_stat==0:
                braf_stat='V600E'
            elif braf_stat==1:
                braf_stat='V600K'
            mutations=table.iloc[ind,41:]
            mutations['BRAF']=1
            mutations=mutations.replace(to_replace=(0,1), value=('no', 'yes'))
            for i in range(0,len(mutations)):
                if (mutations.index[i])=='BRAF':
                    HGVSp = braf_stat
                else:
                    HGVSp = None
                snp_dict=dict(
                    sample_id=sample_prefix+f"{ind:03}",
                    patientID=pat_prefix+"_"+f"{ind:03}",
                    HGNC = mutations.index[i],
                    HGVSp_short = HGVSp,
                    mutated = mutations[i],
                    temporality = 'pre treatment',
                    source = dict(
                        title = 'TERT Promoter Mutation as an Independent Prognostic Marker for Poor Prognosis MAPK Inhibitors-Treated Melanoma',
                        author =  'Pauline Blateau, Jerome Solassol',
                        journal =  'Cancers',
                        location = 'Montpellier (France)',
                        date = 2020)
                )
                list_snps.append(snp_dict)

        if (cond=='patients'):
            print('the list of "Blateau" patients has been created')
            return(list_patients)
        elif(cond=='mutations'):
            print('the list of "Blateau" mutations has been created')
            return(list_snps)
