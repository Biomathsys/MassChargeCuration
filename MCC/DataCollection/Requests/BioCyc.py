import re
import codecs
import json
import requests
import xml.etree.ElementTree as ET
import logging
import pandas as pd
from io import StringIO

from .databaseInterface import DatabaseInterface

replace_capital_ids = re.compile(r"([A-Z])([A-Z])")
remove_1 = re.compile(r"([A-Z][a-z]?)(1)([A-Z]|$)")

class BioCycInterface(DatabaseInterface):
    def __init__(self, data_path, no_local = False, biocyc_base_path = None):
        self.data_path = data_path
        self.no_local = no_local
        self.biocyc_base_path = biocyc_base_path
        self.load_biocyc_db()

    def load_biocyc_db(self):
        """
        Loads the BioCyc database files (and on demand consolidates them).
        """
        if self.no_local: 
            self.BioCyc_dict = {}
            return
        try:
            with open(f"{self.data_path}/BioCyc.json", "r") as f:
                self.BioCyc_dict = json.loads(f.read())
        except FileNotFoundError:
            if not self.biocyc_base_path is None:
                self.consolidate_BioCyc()
                with open(f"{self.data_path}/BioCyc.json", "r") as f:
                    self.BioCyc_dict = json.loads(f.read())
            else:
                self.BioCyc_dict = {}
            
    def get_formulae_by_id(self, meta_id):
        formula = self.BioCyc_dict.get(meta_id.replace("META:", ""), {}).get("formula", None)
        charge = self.BioCyc_dict.get(meta_id.replace("META:", ""), {}).get("charge", None)
        if formula is None:
            result = self.get_BioCyc_information(meta_id)
            if result is None: return result
            else: return [result]
        return [(formula, charge)]

    def get_BioCyc_information(self, meta_id):
        """
        Online information fetching as fallback if no offline data is found.

        Args:
            meta_id (str): Id of the metabolite for which we want to fetch information.
        """
        base_url = 'https://websvc.biocyc.org/getxml?{}&detail=low'
        try:
            answer = requests.get(base_url.format(meta_id))
            logging.debug(f"Could not find entry {meta_id} in local BioCyc file. Trying to fetch from {base_url.format(meta_id)}")
            xml_tree = ET.fromstring(answer.text)
            formula = xml_tree.find(".//formula")
            charge = xml_tree.find(".//Compound/cml/molecule")
            if not charge is None:
                charge = int(charge.attrib['formalCharge'])
            if not formula is None:
                formula = formula.attrib['concise'].replace(" ", "")
                formula = replace_capital_ids.sub(lambda pat: pat.group(1) + pat.group(2).lower(), formula)
                formula = remove_1.sub(r"\1\3", formula)
                return (formula, charge)
        except ET.ParseError: #Dead links return HTML files instead of XML
            return

    def search_identifier(self, names, other_ids):
        id_mapping = {"metanetx.chemical" : "MetaNetX",
                 "bigg.metabolite" : "BiGG",
                 "seed.compound" : "Seed",
                 "sabiork.compound": None,
                 "biocyc" : "BioCyc" # somewhat redundant
                 }
        id_list = [(db, identifier) for db, value in other_ids.items() for identifier in value["ids"]]
        query = "\n".join([*names, *[f"{id_mapping[key]}:{value} " for key, value in id_list]])
        if len(names) == 0: return None
        response = requests.post("https://biocyc.org/META/metabolite-translation-service?=file", files = {"file" : query})
        df = pd.read_csv(StringIO(response.text), sep = "\t")
        if not "Result" in df.columns:
            return []
        else:
            return [f"META:{biocyc_id}" for biocyc_id in df[df["Result"] == "success"].BioCyc.unique()]

    def get_other_references(self, id, relevant_dbs):
        id_mapping = {"METANETX" : "metanetx.chemical", # somewhat redundant
                    "BIGG" : "bigg.metabolite",
                    "SEED" : "seed.compound"
                    }
        result = {}
        for database_id in id_mapping:
            links = self.BioCyc_dict.get(id, {}).get("db_links", {}).get(database_id, None)
            if links is None: continue
            else:
                result[id_mapping[database_id]] = links.split(";")
        return result

    def consolidate_BioCyc(self):
        """
        Function to condense the information we need from the BioCyc database files.
        """
        compound_file_name = f"{self.biocyc_base_path}/compounds.dat"
        class_file_name = f"{self.biocyc_base_path}/classes.dat"
        with codecs.open(compound_file_name, 'r', encoding='utf-8', errors='ignore') as fdata:
            classes = fdata.read().split("//\n")[1:]
            classes_split = [c.split("\n") for c in classes]
            compound_dict = {}
            for split in classes_split:
                meta_id, information = parse_biocyc_compound(split)
                compound_dict[meta_id] = information
        with codecs.open(class_file_name, 'r', encoding='utf-8', errors='ignore') as fdata:
            classes = fdata.read().split("//\n")[1:]
            classes_split = [c.split("\n") for c in classes]
            for split in classes_split:
                meta_id, information = parse_biocyc_class(split)
                compound_dict[meta_id] = information

        with open(f"{self.data_path}/BioCyc.json" ,"w") as f:
            f.write(json.dumps(compound_dict))

find_id = re.compile(r"UNIQUE-ID - (.*)")
find_formula = re.compile(r"CHEMICAL-FORMULA - \((\S+) (\d+)\)")
find_charge = re.compile(r"ATOM-CHARGES - \((\d+) (\S+)\)")
find_name = re.compile(r"(?:COMMON-NAME|SYNONYMS) - (.*)")
find_link = re.compile(r'DBLINKS - \((\S+) "(.*)"')

def parse_biocyc_compound(compound_split):
    """
    Function to parse the biocyc database files for compounds.
    """
    elements = {}
    db_links = {}
    names = []
    charge = 0
    meta_id = ""
    for line in compound_split:
        if not (found := find_id.search(line)) is None:
            meta_id = found.groups()[0]
        if not (found := find_name.search(line)) is None:
            names.append(found.groups()[0].replace("<sup>", "").replace("</sup>", "").replace("<i>", "").replace("</i>", ""))
        if not (found := find_formula.search(line)) is None:
            elements[found.groups()[0]] = int(found.groups()[1])
        if not (found := find_charge.search(line)) is None:
            charge += int(found.groups()[1])
        if not (found := find_link.search(line)) is None:
            db_links[found.groups()[0]] = found.groups()[1]
    return meta_id, {"names" : names,
                    "formula" : dict_to_formula(elements),
                    "charge" : charge,
                    "db_links" : db_links,
                    "type" : "compound"}

def parse_biocyc_class(compound_split):
    """
    Function to parse the BioCyc class files.
    """

    names = []
    meta_id = ""
    for line in compound_split:
        if not (found := find_id.search(line)) is None:
            meta_id = found.groups()[0]
        if not (found := find_name.search(line)) is None:
            names.append(found.groups()[0].replace("<sup>", "").replace("</sup>", "").replace("<i>", "").replace("</i>", ""))
    return meta_id, {"names" : names,
                    "type" : "class"}

