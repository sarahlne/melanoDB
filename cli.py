############# imports
import sys
import re
import os
import pathlib
import pandas as pd
import json
import logging

from datetime import date, datetime
from peewee import Field, IntegerField, DateField, DateTimeField, CharField, ForeignKeyField, Model as PWModel
from peewee import SqliteDatabase

# from db.patients import Patients
# from db.mutations import Mutations
# from db.snp import SNP
# from db.gene_expr import GENE_EXPR
# from db.patients import db

############# import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

############# import Timer
from timeit import default_timer
import time

############# log initiation
log_file = settings_data['log_file']
user = "DBCreator"
# if not os.path.exists(log_file):
#     os.mkdir(log_file)

fh = logging.FileHandler(log_file, mode='w')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(" %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info(f"Hello {user}")
logger.info(f"Create tables ...")

# create data folder if not exist
if not (os.path.exists(settings_data['db_dir'])):
    os.makedirs(settings_data['db_dir'])