from copy import deepcopy
import logging
from . LibSBMLInterface import LibSBMLInterface
from ..core import *

interfaces = {"libsbml" : LibSBMLInterface}

class ModelInterface():
    def __init__(self, model = None):
        if model is None:
            return
        self.interface = self._read_model(model)
        self.metabolites = self._read_metabolites(self.interface)
        self.reactions = self._read_reactions(self.interface)

    def _read_metabolites(self, interface):
        metabolites = {}
        for metabolite_id in interface.get_metabolite_ids():
            name = interface.get_metabolite_name(metabolite_id)
            formula_string = interface.get_metabolite_formula_by_id(metabolite_id)
            charge = interface.get_metabolite_charge_by_id(metabolite_id)
            sbo = interface.get_metabolite_sbo(metabolite_id)
            cv_terms = interface.get_metabolite_cv_terms(metabolite_id)
            notes = interface.get_metabolite_notes(metabolite_id)
            metabolites[metabolite_id] = Metabolite(metabolite_id, name, Formula(formula_string), charge, sbo, set(), cv_terms, notes)
        return metabolites

    def _read_reactions(self, interface):
        reactions = {}
        for reaction_id in interface.get_reaction_ids():
            name = interface.get_reaction_name(reaction_id)
            metabolite_ids = interface.get_reaction_metabolite_ids(reaction_id)
            metabolites = {self.metabolites[metabolite_id] : entry for metabolite_id, entry in metabolite_ids.items()}
            sbo = interface.get_reaction_sbo(reaction_id)
            cv_terms = interface.get_reaction_cv_terms(reaction_id)
            notes = interface.get_reaction_notes(reaction_id)
            reaction = Reaction(reaction_id, name, metabolites, sbo, cv_terms, notes)
            for metabolite in metabolites:
                metabolite.reactions.add(reaction)
            reactions[reaction_id] = reaction
        return reactions

    def write_model(self, filename):
        for metabolite in self.metabolites.values():
            self.interface.write_metabolite(metabolite.id, metabolite.name, str(metabolite.formula), metabolite.charge, metabolite.SBO, metabolite.cv_terms, metabolite.notes)
        
        for reaction in self.reactions.values():
            meta_dict = {metabolite.id : count for metabolite, count in reaction.metabolites.items()}
            self.interface.write_reaction(reaction.id, reaction.name, meta_dict, reaction.SBO, reaction.cv_terms, reaction.notes)

        self.interface.write_model(filename)

    def get_pseudo_reactions(self):
        pseudo_reactions = set()
        for reaction in self.reactions.values():
            # filter out exchange / sinks
            if all(coeff <= 0 for coeff in reaction.metabolites.values()):
                pseudo_reactions.add(reaction)
            if all(coeff >= 0 for coeff in reaction.metabolites.values()):
                pseudo_reactions.add(reaction)
            # filter out growth
            if (int(reaction.sbo) == 629) or ("growth" in reaction.name): 
                pseudo_reactions.add(reaction)
        return pseudo_reactions
            
    def _read_model(self, model):
        interface = None
        # if string is passed we try to read the model using libsbml
        if type(model) == str:
            try:
                import libsbml
                reader = libsbml.SBMLReader()
                document = reader.readSBML(model)
                model = document.getModel()
                interface = interfaces["libsbml"](model)
            except:
                logging.exception("Could not load model.")
        else:
            try:
                import libsbml
                if type(model) == libsbml.Model:
                    model = model
                    interface = interfaces["libsbml"](model)
            except:
                try:
                    import cobra
                    if type(model) == cobra.core.model.Model:
                        model = model
                        interface = interfaces["cobra"](model)
                except:
                    logging.exception("Could not read model.")
        return interface

    def copy(self):
        new_model_interface = ModelInterface()
        new_interface = self.interface.copy()
        new_model_interface.interface = new_interface
        new_model_interface.metabolites = new_model_interface._read_metabolites(new_interface)
        new_model_interface.reactions = new_model_interface._read_reactions(new_interface)
        return new_model_interface
