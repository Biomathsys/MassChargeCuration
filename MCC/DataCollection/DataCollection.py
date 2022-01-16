from __future__ import annotations
from MCC.util import get_sbml_annotations
import numpy as np
import os
from ..util import *
import logging
import time

from .Requests.BiGG import BiGGInterface
from .Requests.BioCyc import BioCycInterface
from .Requests.MetaNetX import MetaNetXInterface
from .Requests.ModelSEED import ModelSEEDInterface


default_interfaces = {"metanetx.chemical" : MetaNetXInterface,
                    "bigg.metabolite" : BiGGInterface,
                    "seed.compound" : ModelSEEDInterface,
                    "biocyc" : BioCycInterface
                    }


class DataCollector:
    """
    Class intended to be used by the ModelBalancer to gather the information for the metabolites from different databases.
    Uses further specified interfaces to different databases to gather the data, cleans the gathered formulae and provides a function to 
    gather the ids for different databases for the metabolites based on the present information.

    The class can be used with online and offline resources, however, full functionality is only given if both can be used.

    Args:
        model (cobrapy.Model): Model for which the DataCollector will gather the information and takes information to gather ids.
        data_path (string): Optional; Path where offline data is searched for/downloaded to.
        update_ids (bool): Optional; Whether or not to gather/update the database identifiers for the used databases.
        gather_information (bool): Optional; If False, will not actually gather any information. Useful if you like 
            to register other interfaces first, or only want to update ids.
        used_annotations ([str]): Optional; Annotations to use, defaults to all. If you distrust a certain database or only want to use a specific one
            you can specify which databases to use here.
        no_local (bool): Optional; Whether or not to use local data/ download database files.
        biocyc_path (str): Optional; Path to biocyc database files.
    """

    def __init__(self, model, data_path = "./data", update_ids = False, gather_information = True, used_annotations = None, no_local = False, biocyc_path = None):
        self.model = model
        self.no_local = no_local
        self.interfaces = {}
        self.used_annotations = list(default_interfaces.keys()) if used_annotations is None else used_annotations
        self.data_path = data_path
        self.strict_linkback = True
        try:
            os.makedirs(self.data_path)
        except PermissionError:
            raise PermissionError
        except OSError:
            pass
        self._load_default_interfaces(data_path, biocyc_path)
        self.assignments = {}
        self.allow_undefined_charge = True
        if update_ids:
            self.get_all_ids()
        if gather_information:
            self.gather_info()

    def _load_default_interfaces(self, data_path, biocyc_path = None):
        """
        Function to load all the default interfaces, specified in default_interfaces at the top of the file.
        Currently compromised of BiGG, MetaNetX, BioCyc and ModelSEED.

        Args:
            data_path: Path for the offline database files.
            biocyc_path: Optinoal; Path for the BioCyc offline database file.
        """
        for identifier, constructor in default_interfaces.items():
            if not (self.used_annotations is None or (identifier in self.used_annotations)): continue
            if identifier == "biocyc":
                interface = constructor(data_path, no_local = self.no_local, biocyc_base_path = biocyc_path)
            else:
                interface = constructor(data_path, no_local = self.no_local)
            self.register_interface(identifier, interface)


    def register_interface(self, identifier, interface):
        """
        Function to register a database interface.

        Args:
            identifier (str): Name of the identifiers.org identifier. Should be the same as in metabolite.annotation.
            interface (MCC.DatabaseInterface): Database interface. Should inherit from MCC.DatabaseInterface.
        """
        self.interfaces[identifier] = interface

    

    def get_formulae(self, metabolite):
        """
        Function to gather all available formulae/charges with the given interfaces.

        Args:
            metabolite (cobrapy.Metabolite): Metabolite for which to gather the formulae/charges.

        Returns:
            Dictionary mapping all found formulae/charges to the containing databases.
                => {(formula, charge): set(database_identifiers)}

        """
        formulae = {}
        annotations = get_sbml_annotations(metabolite)
        for db_id, interface in self.interfaces.items():
            if db_id in annotations:
                ids = [annotations[db_id]] if type(annotations[db_id]) != list else annotations[db_id]
                for identifier in ids:
                    try:
                        if not ((cur_formulae := interface.get_formulae_by_id(identifier)) is None):
                            logging.debug(f"{db_id}, {identifier}")
                            for formula in cur_formulae:
                                if formula is None: continue
                                if (type(formula[0]) == float) or (formula[0] is None): continue
                                if (type(formula[1]) ==float) and np.isnan(formula[1]) or (formula[1] is None):
                                    if self.allow_undefined_charge:
                                        formula = (formula[0], None)
                                    else:
                                        continue
                                formula = (clean_formula(formula[0]), int(formula[1]) if not formula[1] is None else None)
                                cur_db = formulae.get(formula, set())
                                cur_db.add((db_id, identifier))
                                formulae[formula] = cur_db
                    except KeyboardInterrupt:
                        raise
                    except Exception as e:
                        logging.exception(f"Error getting formula for {identifier} in {db_id}:")


        return formulae


    def gather_info(self):
        """
        Gathers formulae and charges for all metabolites in self.model.
        """
        for i, metabolite in enumerate(self.model.getListOfSpecies()):
            logging.info(f"{i + 1}/{self.model.getNumSpecies()}: Getting information for {metabolite.id[2:]}")
            self.assignments[metabolite.id[2:]] = self.get_formulae(metabolite)

    def get_assignments(self, metabolite, clean = True, partial = True, database_seperated = False):
        """
        Function to return all assignments for the given metabolite that were found using all registered interfaces.
        
        Args:
            metabolite (cobrapy.Metabolite): Metabolite for which to gather information.
            clean (bool): Whether or not to also return formulae which were not cleaned (False currently not implemented).
            partial (bool): Whether or not to return formulae containing wildcard symbols (False currently not implemented).
            database_seperated (bool): Whether or not to return formulae mapping to the databases that contain them.
        
        Returns:
            If database_seperated = False: Set of formula/charge combinations that could be found for this metabolite.
            Otherwise: Dictionary mapping (formula, charge) to set of database identifiers where it could be found.

        """
        if len(self.assignments) == 0:
            logging.warn("Tried to get assignments with no gathered information. Try calling gather_info before.")
            return None
        assignments = self.assignments.get(metabolite.id[2:], {})
        if get_sbml_notes(metabolite).get("type", "metabolite") != "class":
            filtered_assignments = {assignment: databases for assignment, databases in assignments.items() if not ("R" in assignment[0])} 
            if len(filtered_assignments) > 0:
                assignments = filtered_assignments
        if (clean == False) or (partial == False):
            raise NotImplementedError
        
        if database_seperated:
            return assignments
        else:
            return set(assignments.keys())
            

    def get_all_ids(self):
        """
        Updates/gathers all ids for all metabolites in the currently registered database interfaces.
        """
        now = time.process_time()
        for i, metabolite in enumerate(self.model.getListOfSpecies()):
            logging.info(f"{i + 1}/{self.model.getNumSpecies()}: {metabolite.id[2:]}")
            ids = self.get_ids(metabolite)
            logging.debug(f"Ids were {ids}")
            annotations = get_sbml_annotations(metabolite)
            for db in self.used_annotations: 
                if db == "biocyc":
                    annotations[db] = [f"META:{entry}" for entry in ids[0][db]["ids"]]
                else:
                    annotations[db] = list(ids[0][db]["ids"])
            logging.debug(f"Updated metabolite {metabolite.id[2:]} annotations to {annotations}")
        logging.info(f"{time.process_time() - now}")

    def get_ids(self, metabolite):
        """
        Updates/gathers all ids for the given metabolite in the currently registered database interfaces.

        Args:
            metaoblite (cobrapy.Metabolite): Metabolite for which to gather the ids.

        Returns:
            Dictionary mapping database identifiers to outdated and current ids.
                => {db_identifer : {"old_ids" : [outdated ids], "ids" : set(current ids)}}
        """

        names = [metabolite.name] if metabolite.name else ([metabolite.id[2:]] if metabolite.id else [])
        DB_ids = {db_identifier : {"old_ids": [], "ids" : set()} for db_identifier in self.used_annotations}
        missing_links = set(self.used_annotations)
        check_list = set()
        # taking ids from annotations
        annotations = get_sbml_annotations(metabolite)
        for db_identifier in annotations:
            if db_identifier in self.used_annotations:
                if type(annotations[db_identifier]) is list:
                    DB_ids[db_identifier]["ids"].update([meta_id.replace("META:", "") for meta_id in annotations[db_identifier]])
                else:
                    DB_ids[db_identifier]["ids"].add(annotations[db_identifier].replace("META:", ""))
                check_list.update([(db_identifier, meta_id.replace("META:", "")) for meta_id in DB_ids[db_identifier]["ids"]])
                missing_links.remove(db_identifier)
        
        # fetching missing ids from other information
        for db_identifier in missing_links:
            try:
                if (not (found := self.interfaces[db_identifier].search_identifier(names, DB_ids)) is None) and (len(found) > 0):
                    DB_ids[db_identifier]["ids"].update(found)
                    logging.info(f"Found new ids {found} in {db_identifier} via id & name based search for {metabolite.id[2:]}")
                    check_list.update([(db_identifier, meta_id.replace("META:", "")) for meta_id in DB_ids[db_identifier]["ids"]])
            except KeyboardInterrupt:
                raise
            except Exception as e:
                logging.exception(f"Error searching for identifier in {db_identifier}:")

        #update identifiers
        flattened_ids = [meta_id.replace("META:", "") for meta_ids in DB_ids.values() for meta_id in meta_ids["ids"]]
        for db_identifier in self.interfaces:
            try:
                old, new = self.interfaces[db_identifier].update_ids(flattened_ids, names)
                check_list.difference_update([(db_identifier, meta_id) for meta_id in old])
                DB_ids[db_identifier]["ids"].difference_update(old)
                check_list.update([(db_identifier, meta_id) for meta_id in new])
                DB_ids[db_identifier]["ids"].update(new)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                logging.exception(f"Error updating identifier in {db_identifier}:")

        # iteratively trying to gather more ids and pruning potential metaNetX ids
        while len(check_list) > 0:
            next_id = check_list.pop()
            try:
                result = self.interfaces[next_id[0]].get_other_references(next_id[1], self.used_annotations)
                if result is None: continue
                for database_id, meta_ids in result.items():
                    if (meta_ids is None) or (database_id not in self.used_annotations): continue
                    for meta_id in meta_ids:
                        meta_id = meta_id.replace("META:", "")
                        if (meta_id not in DB_ids[database_id]["ids"]):
                            if (meta_id in DB_ids[database_id]["old_ids"]):
                                _, new_ids = self.interfaces[database_id].update_ids(meta_id)
                            else:
                                new_ids = [meta_id]
                            for meta_id in new_ids:
                                back_references = self.interfaces[database_id].get_other_references(meta_id, self.used_annotations)
                                back_ref = back_references.get(next_id[0], [])
                                strict_condition = self.strict_linkback or len(back_references) != 0
                                if (not strict_condition) or ((not back_ref is None) and (next_id[1] in back_ref)):
                                    DB_ids[database_id]["ids"].add(meta_id)
                                    check_list.add((database_id, meta_id))
                                    logging.info(f"Found new id {meta_id} in {database_id} from {next_id} for {metabolite.id}")
            except KeyboardInterrupt:
                raise
            except Exception as e:
                logging.exception(f"Error getting other identifiers:")
        return DB_ids, names
