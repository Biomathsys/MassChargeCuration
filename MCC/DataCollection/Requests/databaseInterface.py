import json
import requests
import logging
import os

class DatabaseInterface():
    """
    Class representing an Interface to a database. Could use online and offline resources.
    Provides functions to access the database.
    """
    def __init__(self):
        pass

    def get_assignment_by_id(self, meta_id):
        """
        Returns a list of all formulae/charge combinations which can be found for the given id in this database.

        Args:
            meta_id (str): Metabolite id for which to get the information.

        Returns:
            List of formulae/charge combinations which can be found with the given id in this database.
        """
        return None
    
    def search_identifier(self, names, other_ids):
        """
        Returns a list of ids for a metabolite with the given name or cross references.

        Args:
            names ([str]): Available names of the metabolite to search for.
            other_ids ({identifier: id}): Dictionary mapping identifier.org identifers to respective metabolite ids.
                Can be used to search for entries cross referencing.

        Returns:
            List of identifiers which can be found with the given information.
        """
        return []

    def get_other_references(self, id, relevant_dbs):
        """
        Gets references to other databases found using this interface and the given id.

        Args:
            id (str): Metabolite id for which to find other databases ids.
            relevant_dbs ([str]): List of database identifiers of interest.

        Returns:
            Dictionary mapping identifers.org database identifiers to the respective metabolite ids.
                => {database_id : metabolite_id}
        """
        return {}

    def update_ids(self, ids, names):
        """
        Updates the given ids for the database of this interface.
        Returns a tuple of old and new ids.

        Args:
            ids ([str]): List of currently known ids.
            names ([str]): List of currently known names. 
                If multiple new ids are found, we try to filter them out by name.

        Returns:
            Tuple of outdated and new ids.
                ([outdated_ids], [new_ids])
        """
        return ([], [])