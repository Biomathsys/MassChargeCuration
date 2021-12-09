import re
import codecs
import json
import requests
import cobra
import logging
import pandas as pd

from .databaseInterface import DatabaseInterface

class ModelSEEDInterface(DatabaseInterface):
    def __init__(self, data_path, no_local = False, biocyc_base_path = None):
        self.data_path = data_path
        self.no_local = no_local
        self.biocyc_base_path = biocyc_base_path
        self.load_db()
        self.ModelSEED_dict = {}
        def fill_dict(row):
            self.ModelSEED_dict[row['id']] = (row["formula"], row["charge"])
        self.df.apply(fill_dict, axis = 1) 

    def load_db(self):
        """
        Allows to download the modelSEED database for faster access.
        """
        if self.no_local:
            logging.warning("ModelSEED interface is currently only implemented to download the entire database.")
            self.df = pd.DataFrame()
        try:
            self.df = pd.read_csv(f"{self.data_path}/ModelSEED_compounds.tsv", sep = "\t", dtype=object)
        except FileNotFoundError:
            base_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv"
            result = requests.get(base_url)
            result.raise_for_status()
            with open(f"{self.data_path}/ModelSEED_compounds.tsv", "w") as f:
                f.write(result.text)
            self.df = pd.read_csv(f"{self.data_path}/ModelSEED_compounds.tsv", sep = "\t", dtype=object)

            
    def get_formulae_by_id(self, meta_id):
        result = self.ModelSEED_dict.get(meta_id, None)
        if result is None:
            return result
        else:
            return [result]

    def search_identifier_seed(self, names, other_ids):
        """
        Function to search for an identifier in ModelSEED with the given name and references in other databases.
        
        Args:
            names ([str]): List of known names of the metabolite.
            other_ids ({db_id : [meta_ids]}): Dictionary mapping database identifiers to lists of known identifiers for the metabolite in that database. 
        """
        other_ids = [meta_id.replace("META:", "") for meta_ids in other_ids.values() for meta_id in meta_ids["ids"]]
        found_aliases = self.df["aliases"].dropna()[self.df["aliases"].dropna().apply(lambda aliases: any([(True if re.search(f"(:|;) {re.escape(x)}(\||;|\Z)", aliases) else False) for x in [*names, *other_ids]]))]
        found_ids = self.df["id"][found_aliases.index]
        return list(found_ids)

    def get_other_references(self, id, relevant_dbs):
        id_mapping = {"BiGG" : "bigg.metabolite",
                    "MetaCyc" : "biocyc"  
                    }
        references = {}
        found = self.df[self.df["id"] == id]
        if len(found) == 0:
            logging.warning(f"ModelSEED id: {id} could not be found in the database.")
        out_refs = found["aliases"].values[0]
        if type(out_refs) != str: return references
        else: 
            out_refs = out_refs.split("|")
        def split_into_identifiers(s):
            colon_index = s.find(":")
            if colon_index > -1:
                db_identifier = s[:colon_index]
                meta_id = s[colon_index + 2:]
                if (db_identifier in id_mapping):
                    db_ref = references.get(db_identifier, set())
                    db_ref.update([mid.strip() for mid in meta_id.split(";")])
                    references[id_mapping[db_identifier]] = db_ref
        [split_into_identifiers(s) for s in out_refs]
        
        return references