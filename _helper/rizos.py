import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json


############################################################################################
#
#                                        Rizos class
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

class Rizos():
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
        table = pd.read_excel(file_path)

        # change values of some table columns
        table_sex = [re.split('/', x)[1] for x in table['Age/ sex']]
        table_sex = ['male' if x=='M' else 'female' for x in table_sex]

        table_age = [int(re.split('/', x)[0]) for x in table['Age/ sex']]
        table_DCR = table['RECIST Response category'].replace(['No RES', 'NoRes', 'NoRes', 'RES'], ['PD', 'PD', 'PD','CR'])
        table_BRAF_mut = table['BRAF genotype']
        table_M_stage = 'M1'+table['M stage']
        table_drug = table['Drug'].replace(['Dab','Vem'],['dabrafenib','vemurafenib'])
        table_PFS_month = table['PFS (weeks)']/4
        table_OS_statut = ['dead' if 'e' in x else 'alive' for x in table['OS (weeks)']]
        table_OS_month = [float(re.sub('e', '', x))/4 for x in table['OS (weeks)']]
        table_brain_met = [1 if 'Brain' in x else 0 for x in table['Site']]
        for i in range(len(table_brain_met)):
            if table_brain_met[i] == 0:
                table_brain_met[i] = 'no'
            else:
                table_brain_met[i] = 'yes'

        pat_prefix = 'RIZ_'
        sample_prefix= 'RIZSAM_'
        for ind in table.index:
            patient_dict=dict(
                patient_ID = pat_prefix+str(table['Patient'][ind]),
                sex = table_sex[ind],
                age = table_age[ind],
                stage = 'IV',
                M_stage = table_M_stage[ind],
                LDH = np.NaN,
                os_statut = np.NaN,
                os_months = table_OS_month[ind],
                pfs_statut = np.NaN,
                pfs = table_PFS_month[ind],
                braf_mut = table_BRAF_mut[ind],
                disease_control_rate = table_DCR[ind],
                prelevement_temporality = np.NaN,
                drug = table_drug[ind],
                brain_metastasis = table_brain_met[ind],
                immunotherapy_treatment = np.NaN,
                seq_data = 'no',
                seq_type = np.NaN,
                source = dict(
                    title = 'BRAF Inhibitor Resistance Mechanisms in Metastatic Melanoma: Spectrum and Clinical Impact',
                    author =  'Helen Rizos, Georgina V.Long',
                    journal =  'Clinical Cancer Research',
                    location = 'Sydney (Australia)',
                    date = 2019) 
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
            convert_dict = json.dumps(patient_dict, cls=NpEncoder)

            list_patients.append(json.loads(convert_dict))

        print('the list of "Rizos" patients has been created')
        return(list_patients)
            