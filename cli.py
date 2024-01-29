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

from db.patients import Patients
from db.mutations import Mutations
from db.snv import SNV
from db.cnas import CNA
from db.gene_expr import GENE_EXPR
from db.patients import db

############# import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

############# import Timer
from timeit import default_timer
import time

############# log initiation
log_file = settings_data['log_file']
user = "DBCreator"

fh = logging.FileHandler(log_file, mode='w')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(" %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info(f"Hello {user}")
logger.info(f"Create tables ...")

############# create data folder if not exist
if not (os.path.exists(settings_data['db_dir'])):
    os.makedirs(settings_data['db_dir'])

############# connect db
db.connect()

############# drop tables
db.drop_tables([Patients])
db.drop_tables([CNA])
db.drop_tables([SNV])
db.drop_tables([GENE_EXPR])

############# create tables
db.create_tables([Patients])
db.create_tables([CNA])
db.create_tables([SNV])
db.create_tables([GENE_EXPR])

p1 = Patients()
c1 = CNA()
snv1=SNV()
g1 = GENE_EXPR()

# ----------------- Create Patients ----------------- #
logger.info("Step 1 | Loading Patients ...")
start_time = time.time()
p1.fetch_patients_and_create(settings_data['clinical_studies'])
len_patients = Patients.select().count()
elapsed_time = time.time() - start_time
logger.info("... done in {:10.2f} min for #patients = {}".format(elapsed_time/60, len_patients))

# # ----------------- Create CNAs ----------------- #
logger.info("Step 2 | Loading CNAs ...")
start_time = time.time()
c1.fetch_cnas_and_create(settings_data['clinical_studies'])
len_cnas = CNA.select().count()
elapsed_time = time.time() - start_time
logger.info("... done in {:10.2f} min for #cnas = {}".format(elapsed_time/60, len_cnas))

# ----------------- Create SNPs ----------------- #
logger.info("Step 3 | Loading SNPs ...")
start_time = time.time()
snv1.fetch_snvs_and_create(settings_data['clinical_studies'])
len_snps = SNV.select().count()
print(len_snps)
elapsed_time = time.time() - start_time
logger.info("... done in {:10.2f} min for #snps = {}".format(elapsed_time/60, len_snps))

# # ----------------- Fetch uniprot IDs ----------------- #
logger.info("Step 4 | Register uniprot IDs ...")
start_time = time.time()
snv1.register_uniprot()
elapsed_time = time.time() - start_time
len_uniprot = 13594
logger.info("... done in {:10.2f} min for # unique uniprot ids = {}".format(elapsed_time/60, 13594))

# # ----------------- Create GENE EXPRESSION ----------------- #
logger.info("Step 5 | Loading gene expressions ...")
start_time = time.time()
g1.fetch_gene_expr_and_create(settings_data['clinical_studies'])
len_epxr = GENE_EXPR.select().count()
elapsed_time = time.time() - start_time
logger.info("... done in {:10.2f} min for # gene expressions = {}".format(elapsed_time/60, len_epxr))

# # ----------------- END Creation ----------------- #
creation_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S") 
logger.info(f"Last version of {settings_data['db_name']} database created on {creation_time}")