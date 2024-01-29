import json
import os
from datetime import datetime
from peewee import CharField, FloatField, ForeignKeyField
from peewee import Field, IntegerField, DateField, DateTimeField, CharField, ForeignKeyField, Model as PWModel
from peewee import SqliteDatabase
from peewee import chunked
from playhouse.sqlite_ext import JSONField

#import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

db = SqliteDatabase(settings_data['db_dir']+settings_data['db_name'])
#db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])


class CNA(PWModel):

    """
    This class represents CNAs (Single Nucleotide Polymorphism).

    :property patienr_id: id of the patients
    :type pt_id: CharField 
    :property sample_id: id of the biological sample
    """

    id = IntegerField(primary_key=True)
    creation_datetime = DateTimeField(default=datetime.now)
    patientID = CharField(null=True, index=True)
    sample_id = CharField(null=True, index=True)
    HGNC = CharField(null=True, index=True)
    value = IntegerField(null=True, index=True)
    temporality = CharField(null=True, index=True)
    source = CharField(null=True, index=True)

    def fetch_cnas_and_create(cls, settings_data):
        from _helper.catalanotti import Catalanotti

        c = Catalanotti()
        list_CNAs_catalanotti = c.parse_CNAs_from_file(settings_data['catalanotti_CNA'], file_clinical_sample_path = settings_data['catalanotti_clinical_sample'])

        cls.create_cnas_table_from_list(list_CNAs_catalanotti)

    @classmethod
    def create_cnas_table_from_list(cls, list_cnas):
        with db.atomic():
            for batch in chunked(list_cnas,100):
                cls.insert_many(batch).execute()

    class Meta():
        database = db
        table_name = 'cnas'

        