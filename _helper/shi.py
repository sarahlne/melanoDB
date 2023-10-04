import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json

from math import floor


############################################################################################
#
#                                        SHI class
#                                         
############################################################################################*

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
        

class Shi():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path):

        list_patients = []

        table = pd.read_excel(file_path)

        table_sex = table['Sex'].replace(['M','F'],['male','female'])
        table_age = table['Age']

        table_PFS_month = round(table['Progression\nfree survival\n(days)']/30,2)
        table_stage = ['III' if 'III' in x else 'IV' for x in table['Stage']]
        table_M_stage = [np.NAN if 'III' in x else x.upper() for x in table['Stage']]
        table_drug = ['vemurafenib' if 'Vemu' in x else x for x in table['Dose of\nBRAF\ninhibitor\n(mg)']]
        for i in range(len(table_drug)):
            if 'Dabra' in table_drug[i]:
                table_drug[i]='dabrafenib'

        pat_prefix = 'SHI_'
        sample_prefix = 'SHISAM_'

        #create dictionaries
        for ind in table.index:
            patient_dict=dict(
                patient_ID = pat_prefix+str(table['Patient ID'][ind]),
                sex = table_sex[ind],
                age = table_age[ind],
                stage = table_stage[ind],
                M_stage = table_M_stage[ind],
                LDH = np.NaN,
                os_statut = np.NaN,
                os_months = np.NaN,
                pfs_statut = np.NaN,
                pfs = table_PFS_month[ind],
                braf_mut = np.NAN,
                disease_control_rate = np.NAN,
                prelevement_temporality = np.NaN,
                drug = table_drug[ind],
                brain_metastasis = 'no',
                immunotherapy_treatment = np.NaN,
                seq_data = 'no',
                seq_type = np.NAN,
                source = dict(
                    title = 'Acquired resistance and clonal evolution in melanoma during BRAF inhibitor therapy',
                    author =  'Hubing Shi, Roger S. Lo',
                    journal =  'Cancers Discovery',
                    location = 'Los Angeles (United State)',
                    date = 2014) 
            )
            
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)

            list_patients.append(json.loads(convert_dict))

        print('the list of "Shi" patients has been created')
        return(list_patients)