import requests
import logging
import pandas as pd
from difflib import SequenceMatcher

from .databaseInterface import DatabaseInterface

class MetaNetXInterface(DatabaseInterface):
    def __init__(self, data_path, no_local = False):
        self.data_path = data_path
        self.no_local = no_local
        self.load_metanetx_db()
        self.prop_dict = {}
        def fill_dict(row):
            self.prop_dict[row['#ID']] = (row["formula"], row["charge"])
        self.prop_df.apply(fill_dict, axis = 1)


    def load_metanetx_db(self):
        if self.no_local:
            logging.warning("MetaNetX interface is currently only implemented to download the entire database.")
            self.xref_df = pd.DataFrame()
            self.depr_df = pd.DataFrame()
            self.prop_df = pd.DataFrame()

        # TODO: make skiprows dynamic
        base_url = "http://www.metanetx.org/cgi-bin/mnxget/mnxref/{}"
        try:
            self.xref_df = pd.read_csv(f"{self.data_path}/chem_xref.tsv", sep = "\t", skiprows = 351)
        except FileNotFoundError:
            result = requests.get(base_url.format("chem_xref.tsv"))
            result.raise_for_status()
            with open(f"{self.data_path}/chem_xref.tsv", "w") as f:
                f.write(result.text)
            self.xref_df = pd.read_csv(f"{self.data_path}/chem_xref.tsv", sep = "\t", skiprows = 351)

        try:
            self.depr_df = pd.read_csv(f"{self.data_path}/chem_depr.tsv", sep = "\t", skiprows = 351)
        except FileNotFoundError:
            result = requests.get(base_url.format("chem_depr.tsv"))
            result.raise_for_status()
            with open(f"{self.data_path}/chem_depr.tsv", "w") as f:
                f.write(result.text)
            self.depr_df = pd.read_csv(f"{self.data_path}/chem_depr.tsv", sep = "\t", skiprows = 351)

        try:
            self.prop_df = pd.read_csv(f"{self.data_path}/chem_prop.tsv", sep = "\t", skiprows = 351)
        except FileNotFoundError:
            result = requests.get(base_url.format("chem_prop.tsv"))
            result.raise_for_status()
            with open(f"{self.data_path}/chem_prop.tsv", "w") as f:
                f.write(result.text)
            self.prop_df = pd.read_csv(f"{self.data_path}/chem_prop.tsv", sep = "\t", skiprows = 351)

    def get_formulae_by_id(self, meta_id):
        result = self.prop_dict.get(meta_id, None)
        return [result]

    def search_identifier(self, names, other_ids):
        id_mapping = {"metanetx.chemical" : "mnx", # somewhat redundant
                    "bigg.metabolite" : "bigg.metabolite",
                    "seed.compound" : "seed.compound",
                    "sabiork.compound": "sabiork.compound",
                    "biocyc" : "metacyc.compound" 
                    }
        other_ids = [f'{id_mapping[db_id]}:{meta_id.replace("META:", "")}' for db_id, meta_ids in other_ids.items() for meta_id in meta_ids["ids"]]
        return list(self.xref_df["ID"][self.xref_df["#source"].apply(lambda x: x in other_ids)].unique())

    def update_id(self, id):
        """
        Returns the old and new ids for the given id in MetaNetX.
        Only works with the downloaded depr file.

        Args:
            id: Id to update.

        Returns:
            Tuple of deprecated and new ids.
                => ([deprecated ids], [new ids])
        """
        current_ids = set([id])
        old_ids = set()
        remove = None
        while (remove is None) or len(remove) > 0:
            remove = set()
            new = set()
            for cur_id in current_ids:
                depr_rows = self.depr_df[self.depr_df['#deprecated_ID'] == (cur_id)]
                if len(depr_rows) > 0:
                    remove.add(cur_id)
                    new.update([metabolite_id for metabolite_id in depr_rows["ID"]])
            current_ids.update(new)
            current_ids -= remove
            old_ids.update(remove)
        return old_ids, current_ids 

    def update_ids(self, ids, names):
        new_ids = set()
        names_and_ids = [*names, *ids]
        old_ids = set()
        for meta_id in ids:
            old, new = self.update_id(meta_id)
            # filter for most similar name
            filtered_new = []
            for mid in new:
                found_names = self.prop_df["name"][self.prop_df["#ID"] == mid]
                if len(found_names) == 0: continue
                max_sim = max(found_names.apply(lambda x : max([similar(x, name) for name in names_and_ids])))
                filtered_new.append((max_sim, mid))
            if len(filtered_new) == 0: continue
            max_similarity = max([scored[0] for scored in filtered_new])
            if (max_similarity > .8) and (len(new) > 1):
                filtered_meta_ids = [scored[1] for scored in filtered_new if scored[0] > max_similarity * .9]
                removed_meta_ids = [scored[1] for scored in filtered_new if scored[0] <= max_similarity * .9]
            else:
                filtered_meta_ids = []
                removed_meta_ids = [scored[1] for scored in filtered_new]
                #logging.warning(f"Max metanetX similarity for {metabolite.id} was less then .8 with {filtered_meta_ids} chosen.")
            new_ids.update(filtered_meta_ids)
            new_ids.difference_update(old)
            new_ids.difference_update(removed_meta_ids)
            old_ids.update(removed_meta_ids)
            old_ids.update(old)
        return old_ids, new_ids
            
        ids["metanetx.chemical"]["ids"] = new_ids
        return old, new_ids
        

    def get_other_references(self, id, relevant_dbs):
        id_mapping = {"mnx" : "metanetx.chemical", # somewhat redundant
                 "bigg.metabolite" : "bigg.metabolite",
                 "seed.compound" : "seed.compound",
                 "sabiork.compound": "sabiork.compound",
                 "metacyc.compound" : "biocyc"  
                 }
        references = {}
        out_refs = self.xref_df["#source"][self.xref_df["ID"] == id]
        def split_into_identifiers(s):
            colon_index = s.find(":")
            if colon_index > -1:
                db_identifier = s[:colon_index]
                meta_id = s[colon_index + 1:]
                if (db_identifier in id_mapping):
                    db_ref = references.get(db_identifier, set())
                    db_ref.add(meta_id.replace("META:", ""))
                    references[id_mapping[db_identifier]] = db_ref
        out_refs.apply(split_into_identifiers)
        return references

def similar(a, b):
    """
    Determines similarity between two strings. If a starts with L, than b must also start with L.
    Otherwise we return 0.

    Args:
        a (str): first string to compare.
        b (str): Second string to compare.

    Returns:
        String similarity.
    """
    if a.startswith("L ") or a.startswith("L-"):
        if not(b.startswith("L ") or b.startswith("L-")): return 0
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()