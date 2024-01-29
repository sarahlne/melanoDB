import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json


############################################################################################
#
#                                        Van Allen class
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

class VanAllen():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path):
        """
        Read the csv file extract from the clinical study and returns a list of dictionnaries
        where terms represent proteins filled with their informations (age, sex, stage, etc...). 

        :returns: list of all patients in the study
        :rtype: list
        """
        list_patients = []

        # import table
        table = pd.read_csv(file_path, sep = '\t', header=4)
        table = table.drop([1])
        table=table.reset_index(drop=True)

        table['MEDICATION']=table['MEDICATION'].replace(np.nan, 'unknow')
        #table = table.drop(table.index[53:89],0, inplace=False)
        
        # change values of some table columns 
        table_sex = table['SEX']
        table_sex = table_sex.replace(['Male','Female'],['male', 'female'])

        table_pfs_statut = table['EARLY_RESISTANCE'].replace(['No', 'Yes'],['1','1'])

        # create dictionnaries
        for ind in table.index:
            patient_dict=dict(
                patient_ID = table['PATIENT_ID'][ind],
                original_patientID = table['PATIENT_ID'][ind],
                internal_patientID="VS_"+table['PATIENT_ID'][ind],
                sex = table_sex[ind],
                age = table['AGE'][ind],
                stage = np.NaN,
                LDH = np.NaN,
                os_statut = np.NaN,
                os_months = np.NaN,
                pfs_statut = str(table_pfs_statut[ind]),
                pfs = (table['DURATION_OF_THERAPY_WEEKS'][ind])/4,
                braf_mut = "V600",
                disease_control_rate = table['TREATMENT_BEST_RESPONSE'][ind],
                prelevement_temporality = np.NaN,
                drug = table['MEDICATION'][ind],
                brain_metastasis = np.NaN,
                immunotherapy_treatment = np.NaN,
                immunotherapy_mol = np.NaN,
                CNA='no',
                SNV='yes',
                GEX='no',
                source = dict(
                    title = 'The Genetic Landscape of Clinical Resistance to RAF Inhibition in Metastatic Melanoma',
                    author =  'Eliezer M. Van Allen, Dirk Schadendorf',
                    journal =  'Cancer Discovery',
                    location = 'United States, Germany',
                    date = 2019)  
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None

            list_patients.append(patient_dict)

        print('the list of "Van Allen" patients has been created')
        return(list_patients)

    def parse_mutation_from_file(self, file_mutation_path):
        list_van_allen_mutations = []

        table_mutations_extended = pd.read_csv(file_mutation_path, header=0, sep='\t')

        table_IDs = table_mutations_extended[['Tumor_Sample_Barcode']]
        list_patients = []
        list_treatment = []
        for ind in table_IDs.index:
            list_patients.append(re.match('Pat_\d{2}', table_IDs['Tumor_Sample_Barcode'][ind]).group(0))
            if(re.match('Pat_\d{2}_Pre', table_IDs['Tumor_Sample_Barcode'][ind])):
                list_treatment.append('pre treatment')
            elif(re.match('Pat_\d{2}_Post', table_IDs['Tumor_Sample_Barcode'][ind])):
                list_treatment.append('post treatment')


        table_IDs.insert(loc=0, column='Patient_ID', value=list_patients)
        table_IDs.insert(loc=0, column='Treatment', value=list_treatment)

        #create mutations dictionaries
        for ind in table_mutations_extended.index:
            snp_dict = dict(
                sample_id = table_mutations_extended['Tumor_Sample_Barcode'][ind],
                patientID = "VS_"+table_IDs['Patient_ID'][ind],
                HGNC = table_mutations_extended['Hugo_Symbol'][ind],
                consequence = table_mutations_extended['Consequence'][ind],
                HGVSp = table_mutations_extended['HGVSp'][ind],
                HGVSp_short = table_mutations_extended['HGVSp_Short'][ind],
                variant_classification = table_mutations_extended['Variant_Classification'][ind],
                variant_type = table_mutations_extended['Variant_Type'][ind],
                chromosome = table_mutations_extended['Chromosome'][ind],
                start_position = table_mutations_extended['Start_Position'][ind],
                end_position = table_mutations_extended['End_Position'][ind],
                strand = table_mutations_extended['Strand'][ind],
                ref_allele = table_mutations_extended['Reference_Allele'][ind],
                tumor_allele_1 = table_mutations_extended['Tumor_Seq_Allele1'][ind],
                tumor_allele_2 = table_mutations_extended['Tumor_Seq_Allele2'][ind],
                mutated = 'yes',
                temporality = table_IDs['Treatment'][ind],
                source = dict(
                    title = 'The Genetic Landscape of Clinical Resistance to RAF Inhibition in Metastatic Melanoma',
                    author =  'Eliezer M. Van Allen, Dirk Schadendorf',
                    journal =  'Cancer Discovery',
                    location = 'United States, Germany',
                    date = 2019)
            )

            for (key, value) in snp_dict.items():
                if (pd.isna(value)):
                    snp_dict[key]=None

            convert_dict = json.dumps(snp_dict, cls=NpEncoder)
            list_van_allen_mutations.append(json.loads(convert_dict))

        print('the list of "Van Allen" mutations has been created')
        return(list_van_allen_mutations)