import json
import requests
import os
from ...ModelInterface.ModelInterface import ModelInterface
from .databaseInterface import DatabaseInterface

class BiGGInterface(DatabaseInterface):
    def __init__(self, data_path, no_local = False, biocyc_base_path = None):
        self.data_path = data_path
        self.no_local = no_local
        self.biocyc_base_path = biocyc_base_path
        self.load_db()

    def load_db(self):
        """
        Loads in (or on demand downloads and creates) the BiGG database file.
        """
        if self.no_local:
            self.BiGG_dict = {}
        try:
            with open(f"{self.data_path}/BiGG_Database.json", "r") as f:
                self.BiGG_dict = json.loads(f.read())
        except FileNotFoundError:
            result = requests.get("http://bigg.ucsd.edu/api/v2/models") 
            if result.status_code != 200:
                result.raise_for_status()  # Will only raise for 4xx codes, so...
                raise RuntimeError(f"Request to http://bigg.ucsd.edu/api/v2/models returned status code {result.status_code}")
            base_url = "http://bigg.ucsd.edu/static/models/{}.xml"
            base_path = f"{self.data_path}/BiGG_models"
            try:
                os.makedirs(base_path)
            except OSError:
                pass
            for res in result.json()["results"]:
                model_name = res["bigg_id"]
                with open(f"{base_path}/{model_name}.xml", "w") as file_name:
                    file_name.write(requests.get(base_url.format(model_name)).text)
            self.consolidate_BiGG()

            with open(f"{self.data_path}/BiGG_Database.json", "r") as f:
                self.BiGG_dict = json.loads(f.read())

    def get_assignments_by_id(self, meta_id):
        if not meta_id.startswith("M_"):
            meta_id = f"M_{meta_id}"
        sub_dict = self.BiGG_dict.get(meta_id, {})
        formulae = set()
        for key, value in sub_dict.items():
            if key in ["names", "annotations"] or (value[0] == ""): continue
            else:
                formulae.add(tuple(value))
        if len(formulae) == 0:
            formulae = self.get_BiGG_information(meta_id)
        return formulae 
    
    def search_identifier(self, names, other_ids):
        found = set()
        for key, values in self.BiGG_dict.items():
            if any(name in values['names'] for name in names) or any([(db_id in values['annotations']) and values['annotations'][db_id] == other_ids[db_id] for db_id in other_ids]):
                found.add(key)
        return list(found)

    def get_other_references(self, id, relevant_dbs):
        return {database_id: self.BiGG_dict.get(id, {}).get("annotations", {}).get(database_id, None) for database_id in relevant_dbs}

    def consolidate_BiGG(self):
        """
        Condenses the information contained in all downloaded BiGG models into one file.
        """
        BiGG_Database = {}
        incomplete_models = {}
        for model_name in os.listdir(f"{self.data_path}/BiGG_models")[:3]:
            model_interface = ModelInterface(f"{self.data_path}/BiGG_models/{model_name}")
            for metabolite in model_interface.metabolites.values():
                meta_dict = BiGG_Database.get(metabolite.id[:-2], {})
                if (metabolite.formula is None) or (metabolite.charge is None):
                    model_dict = incomplete_models.get(model_name, {})
                    model_dict[metabolite.id] = (str(metabolite.formula), metabolite.charge)
                    incomplete_models[model_name] = model_dict
                    continue
                meta_dict[model_name[:-4]] = (str(metabolite.formula), metabolite.charge)
                annotations = meta_dict.get("annotations", {})
                names = meta_dict.get("names", set())
                names.add(metabolite.name)
                meta_dict["names"] = names
                sbml_annotations = metabolite.cv_terms
                for anno, value in sbml_annotations.items():
                    annotation = annotations.get(anno, [])
                    annotation.append(value)
                    annotations[anno] = annotation
                meta_dict["annotations"] = annotations
                BiGG_Database[metabolite.id[:-2]] = meta_dict

        for metabolite in BiGG_Database:
            for key, value in BiGG_Database[metabolite]["annotations"].items():
                BiGG_Database[metabolite]["annotations"][key] = list(value)
            BiGG_Database[metabolite]["names"] = list(BiGG_Database[metabolite]["names"])
        with open(f"{self.data_path}/BiGG_Database.json", "w") as db_file:
            db_file.write(json.dumps(BiGG_Database)) 

    

    def get_BiGG_information(self, meta_id):
        """
        Online fallback function if a meta_id cannot be found in the offline databases.
        """
        base_url = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'
        # since generally "M_" is not included in the meta_id in BiGG
        if meta_id.startswith("M_"):
            meta_id = meta_id[2:]
        result = requests.get(base_url.format(meta_id))

        if result.status_code != 200:
            result.raise_for_status()  # Will only raise for 4xx codes, so...
            raise RuntimeError(f"Request to {base_url.format(meta_id)} returned status code {result.status_code}")
        try:
            metabolite_json = result.json()
            charges = metabolite_json['charges']
            formulae = metabolite_json['formulae']
            if len(charges) == 0:
                charges = [None]
            if len(formulae) == 0:
                formulae = [None]
        except json.decoder.JSONDecodeError:
            charges = [None]
            formulae = [None]
        return set((formula, charge) for formula in formulae for charge in charges)
