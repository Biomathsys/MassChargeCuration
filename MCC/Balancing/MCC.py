from ..ReportGeneration.reaction_report import reaction_report
from ..ReportGeneration.metabolite_report import metabolite_report
from ..ReportGeneration.visual_report import visual_report
from ..core import Formula

from ..ModelInterface.ModelInterface import ModelInterface
from .satCore import SatCore
from .formulaOptimizer import FormulaOptimizer
from .adherenceOptimizer import AdherenceOptimizer
from ..DataCollection.DataCollection import DataCollector
from ..util import  adjust_proton_count
import logging
import time
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import re

remove_H = re.compile(r"H\d*([^gfe]|$)") # remember to not also remove H from Hg

class MassChargeCuration:
    """
    Class to balance a model system. This class will take a model balance the model system.

    Args:
        model (str or libsbml.Model or cobra.core.model.Model): Path to the model to read in or model that was read in using libsbml or cobrapy.
        data_collector (DataCollector): DataCollector object to use for data collection. If None a new DataCollector will be created.
        data_path (str): Path to the data folder.
        fixed_assignments (dict): Dictionary of fixed assignments for metabolites.
        fixed_reactions (dict): Dictionary of fixed reactions.
        run_optimization (bool): If True the optimization will be run.
        cache_ids (list): List of ids to cache.

    """
    def __init__(self, model, data_collector = None, data_path = "/data", fixed_assignments = None, fixed_reactions = None, run_optimization = True, cache_ids = None, **kw):
        total_time = time.process_time()
        logging.info(f"Setting up model system.")
        self.model_interface = ModelInterface(model)
        self.pseudo_reactions = self.model_interface.get_pseudo_reactions()
        self.original_model_interface = self.model_interface.copy()
        logging.info(f"Collecting Data...")
        t = time.process_time()
        self.data_collector = DataCollector(model, data_path, cache_ids = cache_ids, **kw) if data_collector is None else data_collector
        logging.info(f"[{time.process_time() - t:.3f} s] Collected Data.")
        self.fixed_assignments = fixed_assignments
        self.proton_adjusted_reactions = {}
        logging.info(f"[{time.process_time() - total_time:.3f} s] Finished setting up model system.")

        # finding unbalancable reactions
        t = time.process_time()
        logging.info(f"Started constructing SMT model.")
        self.balancer = SatCore(self.model_interface, self.data_collector, fixed_assignments, fixed_reactions)
        logging.info(f"[{time.process_time() - t:.3f} s] Finished constructing SMT model.")
        self.balancer.balance()
        # as the model simplification already excludes all inherently unbalancable reactions, we must see which were unbalancable from the beginning
        logging.info(f"[{time.process_time() - t:.3f} s] Finished balancibility check. {len([r.id for r in self.balancer.unbalancable_reactions.difference(self.pseudo_reactions)])} non-pseudo reactions were unbalancable.")

        if run_optimization:
            logging.info(f"Started optimizing model.")
            # optimize formula selections
            t = time.process_time()
            self.optimizer = AdherenceOptimizer(self.balancer, self.original_model_interface)
            self.optimizer.balance()
            logging.info(f"[{time.process_time() - t:.3f} s] Finished adherence optimization.")

           # optimize formula selections
            t = time.process_time()
            self.optimizer = FormulaOptimizer(self.balancer, self.original_model_interface)
            self.optimizer.balance()
            logging.info(f"[{time.process_time() - t:.3f} s] Finished formula optimization.")

        # add wildcards back in for unconstrained formulae
        self.reintroduce_wildcards()

        # try to take representations from original model
        self.fit_to_original()

        # adjusting protons for reactions
        self.adjust_protons()
        self.total_time = time.process_time() - total_time

    def reintroduce_wildcards(self):
        # get any unknown metabolites
        wildcard_metabolites = set()
        for metabolite in self.model_interface.metabolites.values():
            assignments = self.assignments[metabolite.id]
            if any(("R" in formula) for formula, charge in assignments):
                matched_clean = any(metabolite.formula == formula for formula, charge in assignments)
                if not matched_clean:
                    wildcard_metabolites.add(metabolite)
            if (len(assignments) == 0):
                original_metabolite = self.original_model_interface.metabolites[metabolite.id]
                if (not metabolite.formula == original_metabolite.formula) or (not metabolite.charge == original_metabolite.charge):
                        wildcard_metabolites.add(metabolite)
    
        # check for inferred formulae
        unchecked_metabolites = set(wildcard_metabolites)
        inferred = set()
        while len(unchecked_metabolites):
            unchecked_metabolite = unchecked_metabolites.pop()
            reactions_counts = {}
            for reaction in unchecked_metabolite.reactions:
                reactions_counts[reaction] = len([m for m in reaction.metabolites if (m in wildcard_metabolites)])
            if any(count == 1 for count in reactions_counts.values()):
                wildcard_metabolites.remove(unchecked_metabolite)
                inferred.add(unchecked_metabolite)
                for reaction, count in reactions_counts.items():
                    if count > 1:
                        unchecked_metabolites.update(wildcard_metabolites.intersection(reaction.metabolites))
                        unchecked_metabolites.discard(unchecked_metabolite)
        # for inferred formulae we add ECO terms
        # FIXME: reactivate
        """
        for inferred_metabolite in inferred:
            inferred_metabolite.annotation["eco"] = "ECO:0000305"
            inferred_metabolite.notes["inferrence"] = "Inferred formulae"
        """
        # add wildcard symbol
        for metabolite in wildcard_metabolites:
            assignments = self.assignments[metabolite.id]
            matched_clean = False
            for formula, charge in assignments:
                formula_no_R = formula.copy()
                formula_no_R["R"] = 0
                metabolite_formula_no_R = metabolite.formula.copy()
                metabolite_formula_no_R["R"] = 0
                if formula_no_R == metabolite_formula_no_R:
                    metabolite.formula = formula
                    matched_clean = True
            if not matched_clean:
                metabolite.formula["R"] = 1
                

    def fit_to_original(self):
        assignments = self.balancer._calculate_cH_equivalents(reduce = False)
        for metabolite in self.model_interface.metabolites.values():
            original_metabolite = self.original_model_interface.metabolites[metabolite.id]
            nonH_formula = metabolite.formula.copy()
            nonH_formula["H"] = 0
            nonH_original = original_metabolite.formula.copy()
            nonH_original["H"] = 0
            if nonH_formula == nonH_original:
                diff_seperated = assignments[metabolite.id].get(nonH_formula, {})
                if not (metabolite.charge is None):
                    diff = metabolite.formula["H"] - metabolite.charge
                    equivalent_assignments = diff_seperated.get(diff, [])
                else:
                    equivalent_assignments = diff_seperated.get(None, [])
                for formula, charge in equivalent_assignments:
                    if formula == original_metabolite.formula and (charge == original_metabolite.charge):
                        metabolite.formula = formula
                        metabolite.charge = charge

    def adjust_protons(self):
        for reaction in self.model_interface.reactions.values():
            if reaction in self.pseudo_reactions:
                continue
            self.proton_adjusted_reactions[reaction.id] = adjust_proton_count(reaction, self.model_interface)
        
    @property
    def unbalancable_reactions(self):
        return self.balancer.unbalancable_reactions

    @property
    def unknown_metabolites(self):
        return self.balancer.unknown_metabolite_ids

    @property
    def reaction_reasons(self):
        return self.balancer.reaction_reasons

    @property
    def assignments(self):
        return self.balancer.assignments

    def add_unbalancable_reaction(self, reaction_id):
        self.balancer.unbalancable_reactions.add(self.model_interface.reactions[reaction_id])

    
    def generate_visual_report(self, filename=None):
        return visual_report(self, filename)


    def generate_reaction_report(self, filename=None, proton_threshold = 7):
        return reaction_report(self, filename, proton_threshold)
        

    def generate_metabolite_report(self, filename = None, original_model = None, target_model = None):
        return metabolite_report(self, filename, original_model, target_model)