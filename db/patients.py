import json

from datetime import datetime
from peewee import CharField, FloatField, ForeignKeyField
from peewee import Field, IntegerField, DateField, DateTimeField, CharField, ForeignKeyField, Model as PWModel
from peewee import SqliteDatabase
from playhouse.sqlite_ext import JSONField

#import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

#db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])
db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])

class Patients(PWModel):
    """
    This class represents Patients.

    :property pt_id: id of the patients
    :type pt_id: CharField 
    :property sex: sex of the pt term
    :type name: CharField 
    :property namespace: namespace of the pt term
    :type namespace: CharField 
    """

    id = IntegerField(primary_key=True)
    data = JSONField(null=True)
    creation_datetime = DateTimeField(default=datetime.now)
    #save_datetime = DateTimeField()
    type = CharField(null=True, index=True)

    patient_ID = CharField(null=True, index=True)
    sex = CharField(null=True, index=True)
    age = IntegerField(null=True, index=True)
    AJCC_stage = IntegerField(null=True, index=True)
    M_stage = CharField(null=True, index=True)
    LDH = CharField(null=True, index=True)
    OS_statut = CharField(null=True, index=True)
    OS_month = FloatField(null=True, index=True)
    PFS_statut=CharField(null=True, index=True)
    PFS_month = FloatField(null=True, index=True)
    drug = CharField(null=True, index=True)
    disease_control_rate = CharField(null=True, index=True)
    BRAF_mut=CharField(null=True, index=True)
    brain_metastasis = CharField(null=True, index=True)
    immunotherapy_treatment = IntegerField(null=True, index=True)
    sequencing_data = CharField(null=True, index=True)
    sequencing_type = CharField(null=True, index=True)
    source = JSONField(null=True)

    #_table_name = 'Patients'

    def fetch_patients_and_create(cls, settings_data):
        from _helper.blateau import Blateau
        from _helper.catalanotti import Catalanotti
        from _helper.vanallen import VanAllen
        from _helper.yan import Yan
        from _helper.rizos import Rizos
        from _helper.louveau import Louveau
        from _helper.rambow import Rambow
        from _helper.shi import Shi
    

        b = Blateau()
        list_pat_blateau = b.parse_xlsx_from_file(settings_data['blateau_file'], 'patients')

        c = Catalanotti()
        list_pat_catalanotti = c.parse_xlsx_from_file(settings_data['catalanotti_file'], settings_data['catalanotti_clinical_sample'])

        v = VanAllen()
        list_pat_vanallen = v.parse_xlsx_from_file(settings_data['van_allen_file'])

        y = Yan()
        list_pat_yan = y.parse_xlsx_from_file(settings_data['yan_ribas_file'])

        m = Louveau()
        list_pat_louveau = m.parse_xlsx_from_file(settings_data['louveau_file'])
        #print(list_pat_louveau)

        r = Rambow()
        list_pat_rambow = r.parse_xlsx_from_file(settings_data['rambow_clinical'])
        #print(list_pat_rambow)

        riz = Rizos()
        #list_pat_rizos = riz.parse_xlsx_from_file(settings_data['rizos_clinical'])
        #print(list_pat_rizos)
        
        shi = Shi()
        list_pat_shi = shi.parse_xlsx_from_file(settings_data['shi_clinical'])
        #print(list_pat_shi)

        cls.create_patients_table_from_list(list_pat_blateau)
        cls.create_patients_table_from_list(list_pat_catalanotti)
        cls.create_patients_table_from_list(list_pat_vanallen)
        cls.create_patients_table_from_list(list_pat_yan)
        cls.create_patients_table_from_list(list_pat_louveau)
        #cls.create_patients_table_from_list(list_pat_rambow)
        #cls.create_patients_table_from_list(list_pat_rizos)
        cls.create_patients_table_from_list(list_pat_shi)

    @classmethod
    def create_patients_table_from_list(cls, list_pt):
        for dico in list_pt:
            pat = cls(data=dico)
            pat.type = type(pat)
            pat.set_pat_id(dico["patient_ID"])
            pat.set_sex(dico["sex"])
            pat.set_age(dico["age"])
            pat.set_AJCC_stage(dico["stage"])
            if ('M_stage' in pat.data.keys()):
                pat.set_M_stage(pat.data['M_stage'])
            pat.set_LDH(dico['LDH'])
            pat.set_os_statut(dico['os_statut'])
            pat.set_os_month(dico['os_months'])
            pat.set_pfs_statut(dico['pfs_statut'])
            pat.set_pfs(dico['pfs'])
            pat.set_drug(dico['drug'])
            pat.set_DCR(dico['disease_control_rate'])
            pat.set_BRAF_mut(dico['braf_mut'])
            pat.set_brain_met(dico['brain_metastasis'])
            pat.set_immuno_treatment(dico['immunotherapy_treatment'])
            pat.set_sequencing_data(dico['seq_data'])
            pat.set_sequencing_type(dico['seq_type'])
            pat.set_source(dico['source'])
            pat.save()
        #patients = [cls(data = dict_) for dict_ in list_pat]


        # for pat in patients:
        #     pat.set_sex(pat.data["sex"])
        #     pat.set_age(pat.data["age"])
        #     pat.set_stage(pat.data["stage"])
        #     #pat.save()
            
        #cls.save_all(patients)

    def set_pat_id(self, pat_id):
        if(pat_id):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.patient_ID = pat_id

    def set_sex(self, sex):
        if(sex):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.sex = sex

    def set_age(self, age):
        if(age):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.age = age

    def set_AJCC_stage(self, stage):
        if(stage is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.AJCC_stage = stage

    def set_M_stage(self, m_stage):
        if(m_stage is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.M_stage = m_stage

    def set_LDH(self, LDH):
        if(LDH is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            
            self.LDH = LDH

    def set_os_statut(self, os_statut):
        if(os_statut is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.OS_statut = os_statut

    def set_os_month(self, os_month):
        if(os_month):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.OS_month = round(os_month,1)
    
    def set_pfs_statut(self,pfs_statut):
        if(pfs_statut):
            self.PFS_statut = pfs_statut
    
    def set_pfs(self, pfs):
        if(pfs):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.PFS_month = round(pfs,1)

    def set_drug(self, drug):
        if(drug):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            if ('leukemia' in drug):
                drug = 'vemurafenib + cobimetinib'
            self.drug = drug.lower()

    def set_DCR(self, DCR):
        if(DCR):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            if('/' in DCR):
                DCR = DCR.split('/')[0]

            self.disease_control_rate = DCR.replace(" ", "")

    def set_BRAF_mut(self, BRAF_mut):
        if(BRAF_mut is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.BRAF_mut = BRAF_mut
    
    def set_brain_met(self, brain_met):
        if(brain_met is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.brain_metastasis = brain_met

    def set_immuno_treatment(self, immuno_treatment):
        if(immuno_treatment is not None):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.immunotherapy_treatment = immuno_treatment
    
    def set_sequencing_data(self, seq_data):
        if(seq_data):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.sequencing_data = seq_data

    def set_sequencing_type(self, seq_type):
        if(seq_type):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.sequencing_type = seq_type

    def set_source(self, source):
        if(source):
            """
            Sets the name of the go term

            :param name: The name
            :type name: str
            """
            self.source = source
    

    class Meta():
        database = db
        table_name = 'patients'