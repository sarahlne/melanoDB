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

from db.mutations import Mutations

#import settings data
with open('settings.json') as setting_file:
    settings_data=json.load(setting_file)

db = SqliteDatabase(settings_data['db_dir']+settings_data['db_name'])
#db = SqliteDatabase(settings_data['db_dir']+settings_data['db_test_name'])

class SNV(Mutations):
    """
    This class represents SNVs (Single Nucleotide Polymorphism).

    :property patienr_id: id of the patients
    :type pt_id: CharField 
    :property sample_id: id of the biological sample
    """

    id = IntegerField(primary_key=True)
    creation_datetime = DateTimeField(default=datetime.now)
    patientID = CharField(null=True, index=True)
    sample_id = CharField(null=True, index=True)
    HGNC = CharField(null=True, index=True)
    mutated = CharField(null=True, index=True)
    HGVSp = CharField(null=True, index=True)
    HGVSp_short = CharField(null=True, index=True)
    uniprot_id = CharField(null=True, index=True)
    consequence = CharField(null=True, index=True)
    variant_classification = CharField(null=True, index=True)
    variant_type = CharField(null=True, index=True)
    chromosome = IntegerField(null=True, index=True)
    start_position = IntegerField(null=True, index=True)
    end_position = IntegerField(null=True, index=True)
    strand = IntegerField(null=True, index=True)
    ref_allele = CharField(null=True, index=True)
    tumor_allele_1 = CharField(null=True, index=True)
    tumor_allele_2 = CharField(null=True, index=True)
    temporality = CharField(null=True, index=True)
    other_prelevements = CharField(null=True, index=True)
    source = CharField(null=True, index=True)

    def fetch_snvs_and_create(cls, settings_data):
        from _helper.blateau import Blateau
        from _helper.catalanotti import Catalanotti
        from _helper.vanallen import VanAllen
        from _helper.louveau import Louveau
        from _helper.yan import Yan

        b = Blateau()
        list_snv_blateau = b.parse_xlsx_from_file(settings_data['blateau_file'], 'mutations')

        c = Catalanotti()
        list_snv_catalanotti = c.parse_mutations_from_file(settings_data['catalanotti_clinical_mutations_extd'], settings_data['catalanotti_clinical_sample'])
        
        v = VanAllen()
        list_snv_vanallen = v.parse_mutation_from_file(settings_data['van_allen_mutations_extd'])

        l = Louveau()
        list_snv_louveau = l.parse_mutations_from_file(settings_data['louveau_file'])
        

        cls.create_snv_table_from_list(list_snv_blateau)
        cls.create_snv_table_from_list(list_snv_catalanotti)
        cls.create_snv_table_from_list(list_snv_vanallen)
        cls.create_snv_table_from_list(list_snv_louveau)

    #### Method for Uniprot Database requests ############################################################

    POLLING_INTERVAL = 3
    API_URL = "https://rest.uniprot.org"

    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_url( url, **kwargs):
        import requests
        import sys
        response = requests.get(url, **kwargs)
    
        if not response.ok:
            print(response.text)
            response.raise_for_status()
            sys.exit()
        return response

    def check_response(cls, response):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json())
            raise
    
    def submit_id_mapping(cls, from_db, to_db, ids):
        request = requests.post(
            f"{SNV.API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids),"taxId":"9606"},
        )
        cls.check_response(request)
        return request.json()["jobId"]
    
    def get_id_mapping_results_link(cls, job_id):
        url = f"{SNV.API_URL}/idmapping/details/{job_id}"
        request = SNV.session.get(url)
        cls.check_response(request)
        return request.json()["redirectURL"]

    def get_id_mapping_results_search(cls, url):
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        if "size" in query:
            size = int(query["size"][0])
        else:
            size = 10
            query["size"] = size
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()
        request = SNV.session.get(url)
        cls.check_response(request)
        results = cls.decode_results(request, file_format, compressed)
        total = int(request.headers["x-total-results"])
        cls.print_progress_batches(0, size, total)
        for i, batch in enumerate(cls.get_batch(request, file_format, compressed), 1):
            results = cls.combine_batches(results, batch, file_format)
            cls.print_progress_batches(i, size, total)
        if file_format == "xml":
            return cls.merge_xml_results(results)
        return results
    
    def get_batch(cls, batch_response, file_format, compressed):
        batch_url = cls.get_next_link(batch_response.headers)
        while batch_url:
            batch_response = cls.session.get(batch_url)
            batch_response.raise_for_status()
            yield cls.decode_results(batch_response, file_format, compressed)
            batch_url = cls.get_next_link(batch_response.headers)
    
    def get_next_link(cls, headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)
    
    def decode_results(cls, response, file_format, compressed):
        if compressed:
            decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                j = json.loads(decompressed.decode("utf-8"))
                return j
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
            elif file_format == "xlsx":
                return [decompressed]
            elif file_format == "xml":
                return [decompressed.decode("utf-8")]
            else:
                return decompressed.decode("utf-8")
        elif file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        elif file_format == "xlsx":
            return [response.content]
        elif file_format == "xml":
            return [response.text]
        return response.text

    def combine_batches(cls, all_results, batch_results, file_format):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results and batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]
        else:
            return all_results + batch_results
        return all_results

    def print_progress_batches(cls, batch_index, size, total):
        n_fetched = min((batch_index + 1) * size, total)
        print(f"Fetched: {n_fetched} / {total}")

    
    def get_xml_namespace(cls, element):
        m = re.match(r"\{(.*)\}", element.tag)
        return m.groups()[0] if m else ""

    def merge_xml_results(cls, xml_results):
        merged_root = ElementTree.fromstring(xml_results[0])
        for result in xml_results[1:]:
            root = ElementTree.fromstring(result)
            for child in root.findall("{http://uniprot.org/uniprot}entry"):
                merged_root.insert(-1, child)
        ElementTree.register_namespace("", cls.get_xml_namespace(merged_root[0]))
        return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

    


    ########################################################################################################
    @classmethod
    def create_snv_table_from_list(cls, list_mut):
        with db.atomic():
            for batch in chunked(list_mut,100):
                cls.insert_many(batch).execute()

    def register_uniprot(cls):
        #https://rest.uniprot.org/beta/docs/
        WEBSITE_API = 'https://rest.uniprot.org/uniprotkb'
        PROTEINS_API = 'https://https://www.ebi.ac.uk/proteins/api'
        gene_list = []
        for snv in cls.select():
            gene_list.append(snv.HGNC)
        gene_list = list(dict.fromkeys(gene_list))  
        num_gene = len(gene_list)
        print(f'There is {num_gene} unique HGNC')
        job_id_test = cls.submit_id_mapping(from_db="Gene_Name", to_db="UniProtKB-Swiss-Prot", ids=gene_list)
        print(job_id_test)
        search_str = "https://rest.uniprot.org/idmapping/uniprotkb/results/"+job_id_test+"?format=json&size=500"
        print(search_str)
        resp = cls.get_id_mapping_results_search(search_str)
        uniprot_IDs = {}
        for elem in resp['results']:
            if (elem['from'] in uniprot_IDs.keys()):
                if (elem['from'] in elem['to']['uniProtkbId']):
                    uniprot_IDs[elem['from']] = elem['to']['primaryAccession']
            else:
                uniprot_IDs[elem['from']] = elem['to']['primaryAccession']
        
        print(f'Number of results: {len(resp["results"])}')
        print(f'Number of unfound ids: {len(resp["failedIds"])}')
        with db.atomic():
            for snv in cls.select():
                if(snv.HGNC in uniprot_IDs.keys()):
                    HGNC = snv.HGNC
                    snv.set_uniprot_id(uniprot_IDs[HGNC])
                    snv.save()
    
    def set_sample_ID(self, sample_ID):
        if(sample_ID):
            self.sample_id = sample_ID
    
    def set_patient_ID(self, internal_patientID):
        if(internal_patientID):
            self.patientID = internal_patientID
    def set_hgnc(self, hgnc):
        if(hgnc):
            self.HGNC = hgnc
    def set_mutated(self, mutated):
        if(mutated):
            self.mutated = mutated
    def set_consequence(self, consequence):
        if(consequence):
            self.consequence = consequence
    def set_HGVSp(self, HGVSp):
        if(HGVSp):
            if('p' in HGVSp):
                HGVSp = HGVSp.split('p.')[1]
            self.HGVSp = HGVSp
    def set_HGVSp_short(self, HGVSp_short):
        if(HGVSp_short):
            if('p' in HGVSp_short):
                HGVSp_short = HGVSp_short.split('p.')[1]
            self.HGVSp_short = HGVSp_short
    def set_uniprot_id(self, uniprot_id):
        if(uniprot_id):
            self.uniprot_id = uniprot_id
    def set_variant_class(self, variant_class):
            self.variant_classification = variant_class
    def set_variant_type(self, variant_type):
            self.variant_type = variant_type
    def set_chromosome(self, chr):
            self.chromosome = chr
    def set_chromosome(self, chr):
        if(chr):
            self.chromosome = chr
    def set_start_posisition(self, start):
        if(start):
            self.start_position = start
    def set_end_posisition(self, end):
        if(end):
            self.end_position = end
    def set_strand(self, strand):
            if(strand):
                self.strand = strand
    def set_ref_allele(self, ref_allele):
            if(ref_allele):
                self.ref_allele = ref_allele
    def set_tumor_allele1(self, tum_allele1):
            if(tum_allele1):
                self.tumor_allele_1 = tum_allele1
    def set_tumor_allele2(self, tum_allele2):
            if(tum_allele2):
                self.tumor_allele_2 = tum_allele2
    def set_temporality(self, temp):
            if(temp):
                self.temporality = temp
    def set_other_prelev(self, other_prelev):
        if(other_prelev):
            self.other_prelevements = other_prelev
    def set_source(self, source):
        if(source):
            self.source = source


    class Meta():
        database = db
        table_name = 'snvs'