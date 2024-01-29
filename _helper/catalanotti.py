import sys
import re
import os
import pathlib
import pandas as pd
import numpy as np
import json


############################################################################################
#
#                                        Catalanotti class
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


class Catalanotti():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """

    def parse_xlsx_from_file(self, file_path, file_clinical_sample_path, file_baits): #file_path2
        """
        Read the csv file extract from the clinical study and returns a list of dictionnaries
        where terms represent proteins filled with their informations (age, sex, stage, etc...). 

        :returns: list of all patients in the study
        :rtype: list
        """
        list_patients = []

        # import table
        table = pd.read_excel(file_path, header=1)
        idx = table[table.isnull().all(1)].index.tolist()
        table = table.drop(idx)
        table=table.reset_index(drop=True)
        #sup_table = pd.read_excel(file_path2, header=1)
        #table = table.drop(table.index[53:89],0, inplace=False)
        
        # change values of some table columns 
        table_sex = table['Sex']
        table_sex = table_sex.replace(['M','F'],['male', 'female'])

        table_LDH = table['LDH (1=high, 0=WNL)']
        table_LDH = table_LDH.replace([0,1], ['normal', 'elevated'])

        table_OS = table['OS NEW Alive (0) Dead (1)']
        table_OS = table_OS.replace([0,1], ['alive', 'dead'])

        clark_level_dict = {'II': 0, 'III': 1, 'IV': 2, 'V': 3}
        table_stage = table['Stage'].tolist()
        table_stage = [elem[0:3] for elem in table_stage]
        table_stage = [clark_level_dict[elem] for elem in table_stage]

        table['Response'] = table['Response'].replace({'N.E.': np.NaN})

        table_brain_met = table['Brain mets (Yes=1; no=0)']
        table_brain_met = table_brain_met.replace([0,1],['no', 'yes'])

        table_immuno_bool = table['Immunotherapy'].tolist()
        for i in range(len(table_immuno_bool)):
            if(table_immuno_bool[i]=='no'):
                table_immuno_bool[i] = 'no'
            else:
                table_immuno_bool[i] = 'yes'
        dict_prior_treatment = {'High-dose IL-2': 'IL-2 (High-dose IL-2)','ipilimumab(adjuvant trial)' : 'Ipi (ipilimumab(adjuvant trial))','ipilimumab (adjuvant trial)' : 'Ipi (ipilimumab(adjuvant trial))', ' interferon (adjuvant setting)': 'IFN (interferon (adjuvant setting))','MK-3475  ': "MEKi (MK-3475)", 'MK-3475 + Nivolumab': "MEKi (MK-3475)"}
        prior_chemo = ['Chemotheraypy' if x==1 else 'no' for x in table['Chemotherapy (Yes=1; no=0)']]
        prior_immuno = [dict_prior_treatment[x] if x in dict_prior_treatment.keys() else 'no' for x in table['Immunotherapy']]
        list_prior_treatment=[]

        for i in range(0, len(prior_immuno)):
            if(prior_immuno[i]=='no'):
                list_prior_treatment.append(prior_chemo[i])
            else:
                list_prior_treatment.append(prior_immuno[i])

        table_baits = pd.read_csv(file_baits, header=0, sep='\t')
        table_baits = table_baits.rename(columns={'SAMPLE_ID': 'Berger ID'})
        table = pd.merge(table, table_baits, on='Berger ID')

        table_clinical_sample = pd.read_csv(file_clinical_sample_path, header=4, sep='\t')
        

        dict_patient_sample = {}
        for i in range(0,len(table_clinical_sample)):
            dict_patient_sample[table_clinical_sample['SAMPLE_ID'][i]] = table_clinical_sample['PATIENT_ID'][i]

        table_drug = table['Drug']
        table_drug = table_drug.replace(['vem+cob (leukemia patient)'], ['vemurafenib + cobimetinib'])
        
        # create dictionnaries
        for ind in table.index:
            patient_dict=dict(
                patient_ID = dict_patient_sample[table['Berger ID'][ind]],
                original_patientID = dict_patient_sample[table['Berger ID'][ind]],
                internal_patientID="CS_"+dict_patient_sample[table['Berger ID'][ind]],
                sex = table_sex[ind],
                age = table['Age'][ind],
                stage = re.sub("[ABCc]","",table['Stage'][ind]), # remettre le C
                M_stage = str(table['M status'][ind]).upper(),
                sample_ID = table['Berger ID'][ind],
                LDH = table_LDH[ind],
                os_statut = table_OS[ind],
                os_months = table['OS (Months)'][ind],
                pfs_statut = str(int(table['POD on vem 1=YES; 0=NO'][ind])),
                pfs = table['PFS (months)'][ind],
                braf_mut = table['BRAF Mutation'][ind],
                disease_control_rate = table['Response'][ind],
                prelevement_temporality = table['timing'][ind],
                drug = table_drug[ind],
                brain_metastasis = table_brain_met[ind],
                immunotherapy_treatment = table_immuno_bool[ind],
                prior_mapk_treatment = list_prior_treatment[ind],
                immunotherapy_mol = table['Immunotherapy'][ind],
                CNA='yes',
                SNV='yes'+' ('+table['mutations'][ind]+')',
                GEX='no',
                source = dict(
                    title = 'PTEN Loss-of-Function Alterations Are Associated With Intrinsic Resistance to BRAF Inhibitors in Metastatic Melanoma',
                    author =  'Federica Catalanotti, David B. Solit',
                    journal =  'JCO Precision Oncology',
                    location = 'New York (United States)',
                    date = 2017) 
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)

            list_patients.append(json.loads(convert_dict))

        print('the list of "Catalanotti" patients has been created')
        return(list_patients)


    def parse_mutations_from_file(self, file_mutation_path, file_clinical_sample_path):
        list_catalanotti_mutations = []
        
        table_mutations_extended = pd.read_csv(file_mutation_path, header=0, sep='\t')

        table_clinical_sample = pd.read_csv(file_clinical_sample_path, header=4, sep='\t')
        dict_patient_sample = {}
        for i in range(0,len(table_clinical_sample)):
            dict_patient_sample[table_clinical_sample['SAMPLE_ID'][i]] = table_clinical_sample['PATIENT_ID'][i]

        #create mutations dictionaries
        for ind in table_mutations_extended.index:
            snp_dict = dict(
                sample_id = table_mutations_extended['Tumor_Sample_Barcode'][ind],
                patientID = "CS_"+dict_patient_sample[table_mutations_extended['Tumor_Sample_Barcode'][ind]],
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
                temporality = 'pre treatment',
                source = dict(
                    title = 'PTEN Loss-of-Function Alterations Are Associated With Intrinsic Resistance to BRAF Inhibitors in Metastatic Melanoma',
                    author =  'Federica Catalanotti, David B. Solit',
                    journal =  'JCO Precision Oncology',
                    location = 'New York (United States)',
                    date = 2017)
            )

            for (key, value) in snp_dict.items():
                if (pd.isna(value)):
                    snp_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        snp_dict[key]=None

            convert_dict = json.dumps(snp_dict, cls=NpEncoder)
            list_catalanotti_mutations.append(json.loads(convert_dict))

        print('the list of "Catalanotti" mutations has been created')
        return(list_catalanotti_mutations)
    

    def parse_CNAs_from_file(self, file_CNA_path, file_clinical_sample_path):
        list_catalanotti_CNAs = []

        # import table
        table_CNAs = pd.read_csv(file_CNA_path, header=0, sep='\t')

        table_clinical_sample = pd.read_csv(file_clinical_sample_path, header=4, sep='\t')
        dict_patient_sample = {}
        for i in range(0,len(table_clinical_sample)):
            dict_patient_sample[table_clinical_sample['SAMPLE_ID'][i]] = table_clinical_sample['PATIENT_ID'][i]

        list_samples = list(table_CNAs.columns)
        list_samples.remove('Hugo_Symbol')

        #create CNAs dictionaries
        for samp in list_samples:
            for ind in table_CNAs.index:
                cnas_dict= dict(
                    sample_id = samp,
                    patientID = "CS_"+dict_patient_sample[samp],
                    HGNC = table_CNAs['Hugo_Symbol'][ind],
                    value = table_CNAs[samp][ind],
                    temporality = 'pre treatment',
                    source = dict(
                        title = 'PTEN Loss-of-Function Alterations Are Associated With Intrinsic Resistance to BRAF Inhibitors in Metastatic Melanoma',
                        author =  'Federica Catalanotti, David B. Solit',
                        journal =  'JCO Precision Oncology',
                        location = 'New York (United States)',
                        date = 2017)
                )
                for (key, value) in cnas_dict.items():
                    if (pd.isna(value)):
                        cnas_dict[key]=None
                    if (isinstance(value, str)):
                        if(value=='NAN'):
                            cnas_dict[key]=None

                convert_dict = json.dumps(cnas_dict, cls=NpEncoder)
                list_catalanotti_CNAs.append(json.loads(convert_dict))

        print('the list of "Catalanotti" CNAs has been created')
        return(list_catalanotti_CNAs)

