import json

from datetime import datetime
from peewee import CharField, FloatField, ForeignKeyField
from peewee import Field, IntegerField, DateField, DateTimeField, CharField, ForeignKeyField, Model as PWModel
from peewee import SqliteDatabase
from playhouse.sqlite_ext import JSONField

#import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])

class Mutations(PWModel):


     class Meta():
        database = db
        table_name = 'mutations'