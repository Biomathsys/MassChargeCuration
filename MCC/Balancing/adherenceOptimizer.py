from .fullBalancer import FullBalancer
import logging
import z3
from ..util import is_cH_balanced 


class AdherenceOptimizer(FullBalancer):
    """
    Class to optimize the found balancing to adhere as much as possible to the assignments in the given target model.

    Args:
        balancer (MCC.FullBalancer): Balancer whose solution the optimization will be based upon. 
        target_model (cobrapy.Model): Model to adhere to.
    """
    
    def __init__(self, balancer, target_model):
        self.balancer = balancer
        super().__init__(balancer.model, balancer.data_collector, balancer.fixed_assignments, target_model=target_model)
        

    def generate_assertions(self):
        """
        Generates all assertions for the z3 solver. Same as for the balancer but
        extended by adding the optimization (soft) constraints.

        Rechecks for any unbalancable reaction.
        """
        self.relevant_elements = self.balancer.relevant_elements
        self.unbalancable_reactions = self.balancer.unbalancable_reactions.copy()
        for reaction in self.model.reactions:
            if not is_cH_balanced(reaction):
                self.unbalancable_reactions.add(reaction.id)
        self._generate_metabolite_assertions()
        self._generate_reaction_assertions()
        self._add_soft_constraints()

    def _setup_z3(self):
        """
        Function to setup the solver, deviates slighty from the balancer setup, since we cannot minimize the unsat core when optimizing.
        """
        # use optimizer instead of solver
        self.simplifying_tactic = z3.Then("simplify", "solve-eqs")
        self.solver = z3.Optimize()
        z3.set_option("parallel.enable", True)

    def _add_soft_constraints(self):
        """
        Adds optimiziation (soft) constraints to the solver. Specifically adds for every metabolite the soft constraints to have the same
        assignment as in the self.target_model.
        """
        for metabolite in self.model.metabolites:
            original_metabolite = self.target_model.metabolites.get_by_id(metabolite.id)

            # if we can adhere to the already given formula, we will try to
            constraints = []
            for element in self.relevant_elements:
                constraints.append(self.metabolite_symbols[metabolite.id][element] == original_metabolite.elements.get(element, 0))
            constraints.append(self.charge_symbols[metabolite.id] == original_metabolite.charge)
            if metabolite.id == "mn2_c":
                print(constraints)
            self.solver.add_soft(z3.And(constraints))
