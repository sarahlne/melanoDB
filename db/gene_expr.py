import json
import re 
import time 
import zlib 

from datetime import datetime
from peewee import CharField, FloatField, ForeignKeyField
from peewee import Field, IntegerField, DateField, DateTimeField, CharField, ForeignKeyField, Model as PWModel
from peewee import SqliteDatabase
from peewee import chunked
from playhouse.sqlite_ext import JSONField

from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

#import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

db = SqliteDatabase(settings_data['db_dir']+settings_data['db_name'])
#db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])

class GENE_EXPR(PWModel):
    id = IntegerField(primary_key=True)
    creation_datetime = DateTimeField(default=datetime.now)
    patientID = CharField(null=True, index=True)
    sample_id = CharField(null=True, index=True)
    HGNC = CharField(null=True, index=True)
    GeneID = CharField(null=True, index=True)
    description = CharField(null=True, index=True)
    value = FloatField(null=True, index=True)
    temporality = CharField(null=True, index=True)
    source = CharField(null=True, index=True)

    def fetch_gene_expr_and_create(cls, settings_data):
        from _helper.louveau import Louveau
        from _helper.yan import Yan
        from _helper.rambow import Rambow

        l = Louveau()
        y = Yan()
        r = Rambow()

        list_gene_expr_louveau = l.parse_gene_expression_from_file(settings_data['louveau_gene_expr_file'])
        list_gene_expr_yan = y.parse_gene_expression_from_file(settings_data['yan_ribas_gene_expr_file'],settings_data['yan_ribas_IDs_match'],settings_data['yan_ribas_hgnc_match'])
        list_gene_expr_rambow_rizos_long = r.parse_gene_expression_from_file(settings_data['rambow_clinical_rizos&long'],settings_data['rambow_gene_expr_rizos&long'])
        list_gene_expr_rambow_kwong = r.parse_kwong_gene_expression_from_file(settings_data['rambow_clinical_kwong'], settings_data['rambow_gene_expr_file_kwong'])
        list_gene_expr_rambow_shi = r.parse_hugo_gene_expression_from_file(settings_data['rambow_clinical_extd_shi'], settings_data['rambow_gene_expr_shi'])

        cls.create_gene_expr_from_list(list_gene_expr_louveau)
        cls.create_gene_expr_from_list(list_gene_expr_yan)
        cls.create_gene_expr_from_list(list_gene_expr_rambow_rizos_long)
        cls.create_gene_expr_from_list(list_gene_expr_rambow_kwong)
        cls.create_gene_expr_from_list(list_gene_expr_rambow_shi)
        
    
    ############################################################################################
    @classmethod
    def create_gene_expr_from_list(cls, list_expr):
        with db.atomic():
            for batch in chunked(list_expr,100):
                cls.insert_many(batch).execute()
    
    class Meta():
        database = db
        table_name = 'gene_expressions'
    