from .balancer import Balancer
from ..util import get_pseudo_reactions
import logging
from math import prod

class ModelBalancer(Balancer):
    """
    Abstract Class providing functionality to load in a given cobrapy model and associated information 
    for it to be balanced.

    Args:
        model (cobrapy.Model): Model which we intend to balance.
        data_collector: Object which provides a get_assignments function.
        fixed_assignments ({metabolite_id: (formula, charge)}): Optional, Dictionary mapping metabolite ids
            to their respective fixed formula and charge. Will overwrite any information from the data_collector.
        include partial_information (bool): Optional, Whether or not information without charge should be used. Defaults to True.

    """
    def __init__(self, model, data_collector = None, fixed_assignments = None, include_partial_information = True):
        self.model = model
        self.include_partial_information = include_partial_information
        self.data_collector = data_collector
        self.fixed_assignments = {} if fixed_assignments is None else fixed_assignments
        self.reaction_reasons = {}
        self.incomplete_formulae = set()
        self.unbalancable_reactions = get_pseudo_reactions(model)
        self.assignments = self._get_assignments()
        self.unconstrained_formulae = self.find_unconstrained_metabolites()
        super().__init__()
    
    def _get_assignments(self):
        """
        Gets assignment for all metabolites of the models from the given data_collector.
        Checks for incomplete formulae (i.e. formulae containing an R.)

        Returns:
            Dictionary mapping metabolite ids to their found assignments.
        """
        all_assignments = self.fixed_assignments.copy()
        for metabolite in self.model.metabolites:
            if metabolite.id in self.fixed_assignments:
                continue
            elif not (assignments := self.data_collector.get_assignments(metabolite, partial=self.include_partial_information)) is None:
                if any(("R" in assignment[0]) for assignment in assignments):
                    self.incomplete_formulae.add(metabolite)
                all_assignments[metabolite.id] = assignments
            else:
                logging.error(f"No information for {metabolite.id} in DataCollector! You might want to gather information first or set a fixed assignment.")
                raise RuntimeError
        return all_assignments

    def find_unconstrained_metabolites(self):
        # TODO: check if this is still in use?
        unconstrained_formulae = self.incomplete_formulae.copy()
        to_check = self.incomplete_formulae.copy()
        while len(to_check) != 0:
            cur_metabolite = to_check.pop()
            is_constrained = False
            for reaction in cur_metabolite.reactions:
                if len(set(reaction.metabolites.keys()).intersection(unconstrained_formulae)) == 1:
                    is_constrained = True
                    # TODO: detect multimetabolites
            if is_constrained:
                unconstrained_formulae.remove(cur_metabolite)
                for reaction in cur_metabolite.reactions:
                    to_check.update(unconstrained_formulae.intersection(set(reaction.metabolites.keys())))
        return unconstrained_formulae

            


    def _get_balanced_combinations(self, reaction, possible_assignments, is_balanced, accept_any = True):
        """
        Function to find all balanced combinations for a given reaction based on this instances
        assignments attribute and the given "is balanced" function.

        Args:
            reaction (cobrapy reaction): Reaction for which to try the combinations.
            possible_assignments ({metabolite : [assignment]}): Dictionary mapping metabolites to lists of possible assignments to them.
            is_balanced (function(reaction, metabolite assignments) => bool): Function returning true if the reaction can be considered balanced wrt to the given assignment.
            accept any (bool) (Optional): Determines whether unknown information is considered balanced. Defaults to True.

        Returns:
            List of balanced assignments to this reactions metabolites.
        """
        # CURRENTLY ONLY USED IN REACTION SCORING; Could be moved to satCore?
        balanced_combinations = []
        
        def try_combinations(fixed, variable):
            '''
            Recursive function to try every possible combination of formulae. 
            it makes use of balanced_combinations as a closure variable to communicate the results.
            
            '''
            if not variable:
                if "any" in fixed.values() and accept_any:
                    balanced_combinations.append(fixed)
                elif is_balanced(reaction, fixed, accept_any):
                    balanced_combinations.append(fixed)
                    
            else:
                options = variable.pop()
                for option in options["assignment"]:
                    next_fixed = fixed.copy()
                    next_fixed[options['id']] = option
                    try_combinations(next_fixed, variable.copy())
                    
        options = []
        for metabolite in reaction.metabolites:
            if metabolite.id in possible_assignments:
                if (len(possible_assignments[metabolite.id]) == 0) or any(("R" in assignment[0]) for assignment in self.assignments[metabolite.id]):
                    options.append({"id" : metabolite.id, "assignment" : ["any"]})
                else:
                    options.append({"id" : metabolite.id, "assignment" : possible_assignments[metabolite.id]})
            else:
                options.append({"id" : metabolite.id, "assignment" : ["any"]})
        if prod([len(s["assignment"]) for s in  options]) > 1e6:
            logging.debug(f"For {reaction.id} there were too many combinations: {prod([len(s['assignment']) for s in options])}. Not checking for direct balancability.")
            return []
        try_combinations({}, options)

        return balanced_combinations
