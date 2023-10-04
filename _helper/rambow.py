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
#                                        RAMBOW class
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
        

class Rambow():
    """
    This module allows to get list of dictionnaries where terms represents patients
    """
    def parse_xlsx_from_file(self, file_path):
        list_patients = []

        #import table 

        table = pd.read_excel(file_path)
        table = table.drop(25, axis=0)
        table.dropna(axis = 0, how = 'all', inplace = True)
        table = table.reset_index(drop=True)
        table['Patient ID CCR']=table['Patient ID CCR'].astype('str')

        table_ids = pd.read_excel(file_path, 'Microarray Key')
        table_ids['CCR Research ID'] = [x.replace('.0', '') for x in table_ids['CCR Research ID'].astype('str')]
        table_ids['NatureComm ID'] = [x.replace('.0', '') for x in table_ids['NatureComm ID'].astype('str')]

        dict_match_ids_rizos = {}
        dict_match_samples_rizos = {}
        dict_match_ids_long = {}
        dict_match_samples_long = {}
        
        for i in range(0, len(table_ids)):
            if (table_ids['CCR Research ID'][i] != 'nan'):
                if (table_ids['CCR Research ID'][i] not in dict_match_ids_rizos.keys()):
                    dict_match_ids_rizos[table_ids['CCR Research ID'][i]] = table_ids['Patient'][i]
            else:
                if (table_ids['NatureComm ID'][i] not in dict_match_ids_long.keys()):
                    dict_match_ids_long[table_ids['NatureComm ID'][i]] = table_ids['Patient'][i]

        for i in range(0, len(table_ids)):
            if (table_ids['CCR Research ID'][i] != 'nan'):
                if (table_ids['CCR Research ID'][i] not in dict_match_samples_rizos.keys()):
                    dict_match_samples_rizos[table_ids['CCR Research ID'][i]]={table_ids['Bx TYPE'][i]: table_ids['HR array sample_deidentified'][i]}
                else:
                    dict_match_samples_rizos[table_ids['CCR Research ID'][i]][table_ids['Bx TYPE'][i]]= table_ids['HR array sample_deidentified'][i]
            else:
                if (table_ids['NatureComm ID'][i] not in dict_match_samples_long.keys()):
                    dict_match_samples_long[table_ids['NatureComm ID'][i]] = {table_ids['Patient'][i]: table_ids['HR array sample_deidentified'][i]}
                else:
                    dict_match_samples_long[table_ids['NatureComm ID'][i]][table_ids['Patient'][i]]= table_ids['HR array sample_deidentified'][i]


        # change value of some table columns
        table_sex = table['sex'].replace(['M','F'],['male', 'female'])
        #table_age = table['Age at stating trial'],
        table_age = [round(x,0) for x in table['Age at stating trial']]
        table_LDH = table['LDH @ drug start (0=normal, 1=elevated)'].replace([0,1],['normal', 'elevated'])
        table_drug = table['drug type'].replace(['dab','vem','dab/tram'],['dabrafenib','vemurafenib','dabrafenib + trametinib'])

        table_PFS_month = [round(x/30,5) for x in table['PFS (days)']]

        table_OS_statut = [(1 if '>' in x else 0) for x in table['OS (days, > indicates still alive)'].astype('str')]
        table_OS_statut = ['alive' if x==1 else 'dead' for x in table_OS_statut]
        table_OS_month = [re.sub('>', '', x) for x in table['OS (days, > indicates still alive)'].astype('str')]
        table_OS_month = [round(int(x)/30,2) for x in table_OS_month]

        table_DCR = table['RECIST category (CR/RES/No RES if no RECIST) '].replace(['No RES'], ['PD']) ### a modifié

        table_BRAF_mut = table['BRAF genotype']
        table_brain_met = table['Active brain mets @ drug start (0=no, 1=yes)'].replace([0,1], ['no','yes'])
        table_stage = ['IV'] * len(table)
        table_M_stage = 'M1'+table['M stage at drug start']

        #create dictionaries
        for ind in table.index:
            if (ind < 24):
                patient_id = dict_match_ids_rizos[table['Patient ID CCR'][ind]]
                if (table['Patient ID CCR'][ind] in dict_match_samples_rizos.keys()):
                    samples = dict_match_samples_rizos[table['Patient ID CCR'][ind]]
                    temp = list(dict_match_samples_rizos[table['Patient ID CCR'][ind]].keys())
                else:
                    samples=np.NAN
                    temp=np.NAN
                article = dict(
                    title = 'BRAF Inhibitor Resistance Mechanisms in Metastatic Melanoma: Spectrum and Clinical Impact',
                    author =  'Helen Rizos, Georgina V.Long',
                    journal =  'Clinical Cancer Research',
                    location = 'Sydney (Australia)',
                    date = 2014)
            else:
                patient_id = dict_match_ids_long[table['Patient ID CCR'][ind]]
                if (table['Patient ID CCR'][ind] in dict_match_samples_long.keys()):
                    samples = dict_match_samples_long[table['Patient ID CCR'][ind]]
                    temp = list(dict_match_samples_long[table['Patient ID CCR'][ind]].keys())
                else:
                    samples=np.NAN
                    temp=np.NAN
                article=dict(
                    title = 'Increased MAPK reactivation in early resistance to dabrafenib/trametinib combination therapy of BRAF-mutant metastatic melanoma',
                    author =  'Georgina V.Long, Helen Rizos',
                    journal =  'Nature Communications',
                    location = 'Sydney (Australia)',
                    date = 2014)
                
            patient_dict=dict(
                patient_ID = patient_id,
                samples=samples,
                sex = table_sex[ind],
                age = table_age[ind],
                stage = table_stage[ind],
                M_stage = table_M_stage[ind],
                LDH = table_LDH[ind],
                os_statut = table_OS_statut[ind],
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
                source = article
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)
            list_patients.append(json.loads(convert_dict))
        print('the list of "Rambow (Sydney, Rizos&Long)" patients has been created')
        return(list_patients)
    
    def parse_kwong_xlsx_from_file(self, file_path, file_supp_path):
        list_patients= []
        # Get supplement data information
        table_supp=pd.read_excel(file_supp_path, header=3)
        table_supp= table_supp.iloc[:26,:]
        table_supp=table_supp.dropna(subset=['A.1','B.1','C.1'], how='all')
        table_supp=table_supp[['PT','Age','RX','Site of disease', 'Response', 'Time To Progression (mos)']].rename(columns={'PT': 'Patient'})
        
        # Get main clinical data informations
        table_pat=pd.read_excel(file_path, 'key')
        dict_samples = {}
        for ind in table_pat.index:
            if (table_pat['Patient'][ind] in dict_samples.keys()):
                dict_samples[table_pat['Patient'][ind]] = dict_samples[table_pat['Patient'][ind]] + [table_pat['RPKM'][ind]]
            else:
                dict_samples[table_pat['Patient'][ind]] = [table_pat['RPKM'][ind]]
        table_pat=table_pat.dropna(subset=['sex'])

        table_clin_extended = table_pat.merge(table_supp, how = 'left', on = 'Patient')

        list_brain_met = ['yes' if 'br' in  table_clin_extended['Site of disease'][ind] else 'no' for ind in table_clin_extended.index]
        list_pfs = [np.round(table_clin_extended['Time To Progression (mos)'][ind],1) if type(table_clin_extended['Time To Progression (mos)'][ind])!=str else np.round(float(table_clin_extended['Time To Progression (mos)'][ind][0:2]),1) for ind in table_clin_extended.index]

        for ind in range(0,len(table_clin_extended.index)):
            patient_dict=dict(
                patient_ID = table_clin_extended['Patient'][ind],
                samples=np.NaN,
                sex = table_clin_extended['sex'].replace({'M': 'male', 'F': 'female'})[ind],
                age = np.round(table_clin_extended['age'][ind]),
                stage = 'IV',
                M_stage = np.NaN,
                LDH = np.NaN,
                os_statut = table_clin_extended['status'][ind],
                os_months = np.round(table_clin_extended['OS_days_1jan2018'][ind]/30.5,1),
                pfs_statut = '1',
                pfs = list_pfs[ind],
                braf_mut = np.NaN,
                disease_control_rate = table_clin_extended['Response'][ind][0:2],
                prelevement_temporality = np.NaN,
                drug = table_clin_extended['RX_short'].replace({'vem':'vemurafenib', 'dab+tra': 'dabrafenib + trametinib'})[ind],
                brain_metastasis = list_brain_met[ind],
                immunotherapy_treatment = 'no',
                seq_data = 'no',
                seq_type = np.NAN,
                source = dict(
                    title = 'Co-clinical assessment identifies patterns of BRAF inhibitor resistance in melanoma',
                    author =  'Lawrence N. Kwong, Lynda Chin',
                    journal =  'Journal of Clinical Investigations ',
                    location = 'Boston (United States)',
                    date = 2015)
            )
            for (key, value) in patient_dict.items():
                if (pd.isna(value)):
                    patient_dict[key]=None
                if (isinstance(value, str)):
                    if(value=='NAN'):
                        patient_dict[key]=None

            convert_dict = json.dumps(patient_dict, cls=NpEncoder)
            list_patients.append(json.loads(convert_dict))
        print('the list of "Rambow (Boston, Kwong)" patients has been created')
        return(list_patients)

    def parse_shi_xlsx_from_file(self, file_path, file_extend_path):
        list_patients = []
        # Get main clinical data
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

        table['Best\noverall\nresponse'] = table['Best\noverall\nresponse']*100
        table['Best\noverall\nresponse'] = np.round(table['Best\noverall\nresponse'])
        
        list_DCR = []

        for i in range(0,len(table['Best\noverall\nresponse'])):
            if(table['Best\noverall\nresponse'][i]==-100):
                list_DCR.append('CR')
            elif(table['Best\noverall\nresponse'][i]<=-30):
                list_DCR.append('PR')
            elif(table['Best\noverall\nresponse'][i]>20):
                list_DCR.append('PD')
            else:
                list_DCR.append('SD')

            

        # Get extended clinical data
        table_extend = pd.read_excel(file_extend_path, 'A-ClinicalTable', skiprows=3 ,skipfooter=13)
        cols_to_remove = [p for p in list(table_extend.columns) if 'Unnamed' in p]
        table_extend.drop(cols_to_remove, axis=1, inplace=True)
        table_extend = table_extend.set_index('Patient initial ID')
        table_extend.index = pd.Series(table_extend.index).fillna(method='ffill')
        table_extend['Patient '] = pd.Series(table_extend['Patient ']).fillna(method='ffill')
        table_extend['Patient '] = table_extend['Patient '].astype(int)
        table_extend['Patient '] = table_extend['Patient '].astype(str)
        table_extend['Patient ']= ['Pt' + a for a in table_extend['Patient ']]
        table_extend['Biopsy timing'] = table_extend['Biopsy timing'].replace({'B': 'baseline'})
        ind_to_keep = [p for p in list(table_extend.index) if 'Shi' in p and '&' not in p]
        ind_to_keep = list(dict.fromkeys(ind_to_keep))
        table_extend = table_extend.loc[ind_to_keep]
        list_samples = list(zip(table_extend.index, table_extend['Biopsy timing']))
        dict_samples = {}
        for i in list_samples:  
            dict_samples.setdefault(i[0],[]).append(i[1])

        list_pat = list(set(list(zip(table_extend.index, table_extend['Patient ']))))
        dict_pat = {}
        for i in list_pat:  
            dict_pat.setdefault(i[0],[]).append(i[1])

        list_pat_samples = list(zip(table_extend['Patient '], table_extend['Biopsy timing']))
        list_pat_samples = [k+'_'+v for k,v in list_pat_samples]


        list_all = list(zip(table_extend.index, table_extend['Patient '], table_extend['Biopsy timing']))
        list_all=[(k,v1+'_'+v2) for k, v1, v2 in list_all]
        dict_all={}
        for i in list_all:
            dict_all.setdefault(i[0],[]).append(i[1])

        pat_prefix = 'SHI_'
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
                pfs_statut = '1',
                pfs = table_PFS_month[ind],
                braf_mut = np.NAN,
                disease_control_rate = list_DCR[ind],
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

        print('the list of "Rambow (LA, Shi&Lo)" patients has been created')
        return(list_patients)
    

    def parse_gene_expression_from_file(self, file_path, file_gene_expr):
        list_ramb_gene_expr = []

        table_ids = pd.read_excel(file_path, 'Microarray Key')
        table_ids['CCR Research ID'] = [x.replace('.0', '') for x in table_ids['CCR Research ID'].astype('str')]
        table_ids['NatureComm ID'] = [x.replace('.0', '') for x in table_ids['NatureComm ID'].astype('str')]

        dict_match_ids_rizos = {}
        dict_match_samples_rizos = {}
        dict_match_ids_long = {}
        dict_match_samples_long = {}
        
        for i in range(0, len(table_ids)):
            if (table_ids['CCR Research ID'][i] != 'nan'):
                if (table_ids['CCR Research ID'][i] not in dict_match_ids_rizos.keys()):
                    dict_match_ids_rizos[table_ids['CCR Research ID'][i]] = table_ids['Patient'][i]
            else:
                if (table_ids['NatureComm ID'][i] not in dict_match_ids_long.keys()):
                    dict_match_ids_long[table_ids['NatureComm ID'][i]] = table_ids['Patient'][i]

        for i in range(0, len(table_ids)):
            if (table_ids['CCR Research ID'][i] != 'nan'):
                if (table_ids['CCR Research ID'][i] not in dict_match_samples_rizos.keys()):
                    #dict_match_samples_rizos[table_ids['CCR Research ID'][i]]={table_ids['Bx TYPE'][i]: table_ids['HR array sample_deidentified'][i]}
                    dict_match_samples_rizos[table_ids['CCR Research ID'][i]]=[table_ids['HR array sample_deidentified'][i]]
                else:
                    #dict_match_samples_rizos[table_ids['CCR Research ID'][i]][table_ids['Bx TYPE'][i]]= table_ids['HR array sample_deidentified'][i]
                    dict_match_samples_rizos[table_ids['CCR Research ID'][i]].append(table_ids['HR array sample_deidentified'][i])
            else:
                if (table_ids['NatureComm ID'][i] not in dict_match_samples_long.keys()):
                    #dict_match_samples_long[table_ids['NatureComm ID'][i]] = {table_ids['Bx TYPE'][i]: table_ids['HR array sample_deidentified'][i]}
                    dict_match_samples_long[table_ids['NatureComm ID'][i]] = [table_ids['HR array sample_deidentified'][i]]
                else:
                    #dict_match_samples_long[table_ids['NatureComm ID'][i]][table_ids['Bx TYPE'][i]]= table_ids['HR array sample_deidentified'][i]
                    dict_match_samples_long[table_ids['NatureComm ID'][i]].append(table_ids['HR array sample_deidentified'][i])

        for num in list(dict_match_ids_rizos.keys()):
            if num in dict_match_samples_rizos:
                p=dict_match_ids_rizos[num]
                dict_match_ids_rizos[p]=dict_match_samples_rizos[num]
                del dict_match_ids_rizos[num]
        for num in list(dict_match_samples_long.keys()):
            if num in dict_match_samples_long:
                p=dict_match_ids_long[num]
                dict_match_ids_long[p] = dict_match_samples_long[num]
                del dict_match_ids_long[num]

        table_gene_expr = pd.read_excel(file_gene_expr, 'Sheet1')

        for j in range(2, len(table_gene_expr.columns)):
                samp_ID = list(table_gene_expr.columns)[j]
                pat_ID=None
                src= None
                for key, val in dict_match_ids_rizos.items():
                    if(samp_ID in val):
                        pat_ID=key
                        src=dict(
                                title = 'BRAF Inhibitor Resistance Mechanisms in Metastatic Melanoma: Spectrum and Clinical Impact',
                                author =  'Helen Rizos, Georgina V.Long',
                                journal =  'Clinical Cancer Research',
                                location = 'Sydney (Australia)',
                                date = 2014)
                        break
                if not (pat_ID):
                    for key, val in dict_match_ids_long.items():
                        if(samp_ID in val):
                            pat_ID=key
                            src=dict(
                                    title = 'Increased MAPK reactivation in early resistance to dabrafenib/trametinib combination therapy of BRAF-mutant metastatic melanoma',
                                    author =  'Georgina V.Long, Helen Rizos',
                                    journal =  'Nature Communications',
                                    location = 'Sydney (Australia)',
                                    date = 2014)
                if('Pre' in samp_ID):
                    temp = 'pre treatment'
                elif('Prog' in samp_ID):
                    temp = 'progression'
                elif('Post' in samp_ID):
                    temp = 'post treatment'
                else:
                    temp=np.NaN
                print(samp_ID)
                print(pat_ID)
                print(temp)
                for i in range(0, len(table_gene_expr)):
                    val = table_gene_expr.iloc[i,j]
                    hgnc = table_gene_expr.Name[i]
                    desc = table_gene_expr['Description-deidentified'][i]
                    gene_expr_dict= dict(
                        patient_ID = pat_ID,
                        sample_ID = samp_ID,
                        HGNC = hgnc,
                        gene_ID= np.NaN,
                        description = desc,
                        value = val, 
                        temporality = temp,
                        source = src
                    )
                    for (key, value) in gene_expr_dict.items():
                            if (pd.isna(value)):
                                gene_expr_dict[key]=None
                            if (isinstance(value, str)):
                                if(value=='NAN'):
                                    gene_expr_dict[key]=None
                    convert_dict = json.dumps(gene_expr_dict, cls=NpEncoder)
                    list_ramb_gene_expr.append(json.loads(convert_dict))
        print('the list of "Rambow" gene expressions has been created')
        return(list_ramb_gene_expr)
    
    def parse_kwong_gene_expression_from_file(self, file_path, file_gene_expr):
        list_kwong_gene_expr = []
        # Get main clinical data informations
        table_pat=pd.read_excel(file_path, 'key')
        dict_samples = {}
        for ind in table_pat.index:
            if (table_pat['Patient'][ind] in dict_samples.keys()):
                dict_samples[table_pat['Patient'][ind]] = dict_samples[table_pat['Patient'][ind]] + [table_pat['RPKM'][ind]]
            else:
                dict_samples[table_pat['Patient'][ind]] = [table_pat['RPKM'][ind]]

        table_gene_expr = pd.read_csv(file_gene_expr, sep='\t', lineterminator='\n', encoding='latin1')
        table_gene_expr['Symbol']=table_gene_expr['Symbol'].replace({' ': np.nan})
        table_gene_expr = table_gene_expr.dropna(subset=('Symbol')).reset_index(drop=True)

        dict_temp_samples= {}
        for j in range(1, len(table_gene_expr.columns)):
            if('A' in table_gene_expr.columns[j]):
                dict_temp_samples[table_gene_expr.columns[j]]='pre treatment'
            elif('B' in table_gene_expr.columns[j]):
                dict_temp_samples[table_gene_expr.columns[j]]='on treatment'
            elif('C' in table_gene_expr.columns[j]):
                dict_temp_samples[table_gene_expr.columns[j]]='progression'
            else:
                dict_temp_samples[table_gene_expr.columns[j].replace('\r','')]=np.NAN

        for j in range(1, len(table_gene_expr.columns)):
            p=np.NaN   
            for key, val in dict_samples.items():
                if(table_gene_expr.columns[j] in val):
                    pat_ID = list(dict_samples.keys())[list(dict_samples.values()).index(val)]
                    p=pat_ID
                else:
                    pat_ID=np.NaN
            samp_ID=table_gene_expr.columns[j].replace('\r','')

            for i in range(0, len(table_gene_expr)):
                val = table_gene_expr.iloc[i,j]
                hgnc = table_gene_expr['Symbol'][i]
                if(p):
                    pat_ID=p
                else:
                    pat_ID=np.NaN
                gene_expr_dict= dict(
                    patient_ID = pat_ID,
                    sample_ID = samp_ID,
                    HGNC = hgnc,
                    gene_ID= np.NaN,
                    description = np.NaN,
                    value = val, 
                    temporality = dict_temp_samples[samp_ID],
                    source = dict(
                                title = 'Co-clinical assessment identifies patterns of BRAF inhibitor resistance in melanoma',
                                author =  'Lawrence N. Kwong, Lynda Chin',
                                journal =  'Journal of Clinical Investigations ',
                                location = 'Boston (United States)',
                                date = 2015)
                )
                for (key, value) in gene_expr_dict.items():
                        if (pd.isna(value)):
                            gene_expr_dict[key]=None
                        if (isinstance(value, str)):
                            if(value=='NAN'):
                                gene_expr_dict[key]=None
                convert_dict = json.dumps(gene_expr_dict, cls=NpEncoder)
                list_kwong_gene_expr.append(json.loads(convert_dict))
        print('the list of "Rambow (Boston)" gene expressions has been created')
        return(list_kwong_gene_expr)
    

    def parse_shi_gene_expression_from_file(self, file_extend_path, file_gene_expr):
        list_shi_gene_expr = []
        # Get extended clinical data
        table_extend = pd.read_excel(file_extend_path, 'A-ClinicalTable', skiprows=3 ,skipfooter=13)
        cols_to_remove = [p for p in list(table_extend.columns) if 'Unnamed' in p]
        table_extend.drop(cols_to_remove, axis=1, inplace=True)
        table_extend = table_extend.set_index('Patient initial ID')
        table_extend.index = pd.Series(table_extend.index).fillna(method='ffill')
        table_extend['Patient '] = pd.Series(table_extend['Patient ']).fillna(method='ffill')
        table_extend['Patient '] = table_extend['Patient '].astype(int)
        table_extend['Patient '] = table_extend['Patient '].astype(str)
        table_extend['Patient ']= ['Pt' + a for a in table_extend['Patient ']]
        table_extend['Biopsy timing'] = table_extend['Biopsy timing'].replace({'B': 'baseline'})
        ind_to_keep = [p for p in list(table_extend.index) if 'Shi' in p and '&' not in p]
        ind_to_keep = list(dict.fromkeys(ind_to_keep))
        table_extend = table_extend.loc[ind_to_keep]

        list_all = list(zip(table_extend.index, table_extend['Patient '], table_extend['Biopsy timing']))
        list_all=[(k,v1+'-'+v2) for k, v1, v2 in list_all]
        dict_all={}
        for i in list_all:
            dict_all.setdefault(i[0],[]).append(i[1])
        samples = [v for k,v in list_all]
        table_gene_expr = pd.read_csv(file_gene_expr, sep='\t', encoding='latin1').set_index('Gene')
        
        for j in range(0, len(table_gene_expr.columns)):
            if(table_gene_expr.columns[j] in samples):
                p=np.NaN   
                for key, val in dict_all.items():
                    if(table_gene_expr.columns[j] in val):
                        pat_ID = list(dict_all.keys())[list(dict_all.values()).index(val)]
                        p=pat_ID
                    else:
                        pat_ID=np.NaN
                samp_ID = table_gene_expr.columns[j]
                if('baseline' in samp_ID):
                    temp = 'pre treatment'
                elif ('DP' in samp_ID):
                    temp = 'progression'
                else:
                    temp= np.NAN
                for i in range(0, len(table_gene_expr)):
                    val = table_gene_expr.iloc[i,j]
                    hgnc = table_gene_expr.index[i]
                    if(p):
                        pat_ID=p
                    else:
                        pat_ID=np.NaN
                    gene_expr_dict= dict(
                        patient_ID = pat_ID,
                        sample_ID = samp_ID,
                        HGNC = hgnc,
                        gene_ID= np.NaN,
                        description = np.NaN,
                        value = val, 
                        temporality = temp,
                        source = dict(
                            title = 'Acquired resistance and clonal evolution in melanoma during BRAF inhibitor therapy',
                            author =  'Hubing Shi, Roger S. Lo',
                            journal =  'Cancers Discovery',
                            location = 'Los Angeles (United State)',
                            date = 2014)
                    )
                    for (key, value) in gene_expr_dict.items():
                            if (pd.isna(value)):
                                gene_expr_dict[key]=None
                            if (isinstance(value, str)):
                                if(value=='NAN'):
                                    gene_expr_dict[key]=None
                    convert_dict = json.dumps(gene_expr_dict, cls=NpEncoder)
                    list_shi_gene_expr.append(json.loads(convert_dict))
        print('the list of "Rambow (LA, Shi&Lo)" gene expressions has been created')
        return(list_shi_gene_expr)