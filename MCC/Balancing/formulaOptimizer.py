from .fullBalancer import FullBalancer
import logging
import z3
import re

remove_H = re.compile(r"H\d*([^gfe]|$)") # remember to not also remove H from Hg
remove_H = re.compile(r"H\d*([^gfe]|$)") # remember to not also remove H from Hg

class FormulaOptimizer(FullBalancer):
    """
    Optimizer to try to choose assignments based on the following criteria:

    - highest possible atom count

    - if unconstrained:

        -> adhering to an unconstrained database representation

        -> lowest atom count possible

    Only assignments which are not the same as in the target model are further optimized.

    Args:
        balancer (MCC.FullBalancer): Balancer whose solution the optimization will be based upon. 
        target_model (core.ModelInterface): ModelInterface of the model to adhere to.
    """
    
    def __init__(self, balancer, target_model_interface):
        self.balancer = balancer
        super().__init__(balancer.model_interface, balancer.data_collector, balancer.fixed_assignments, target_model_interface=target_model_interface)
        

    def generate_assertions(self):
        """
        Generates all assertions for the z3 solver. Same as for the balancer but
        extended by adding the optimization (soft) constraints.

        Rechecks for any unbalancable reaction.
        """
        self.relevant_elements = self.balancer.relevant_elements
        self.unbalancable_reactions = self.balancer.unbalancable_reactions.copy()
        self.unknown_metabolite_ids = self.balancer.unknown_metabolite_ids.copy()
        for reaction in self.balancer.model_interface.reactions.values():
            if not reaction.is_balanced():
                self.unbalancable_reactions.add(reaction)
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

    def _generate_metabolite_assertion(self, metabolite):
        """
        Function the generate the assertion for a metabolite. Extended by adhering to a formula if it is the same formula as in self.target_model_interface,
        as we do not want to optimize these further.
        """
        element_symbols = {}
        constraints = []
        self.charge_symbols[metabolite.id] = z3.Int(f"charge_{metabolite.id}")
        for element in self.relevant_elements:
            element_symbols[element] = z3.Int(f"{element}_{metabolite.id}")
        self.metabolite_symbols[metabolite.id] = element_symbols
        original_metabolite = self.target_model_interface.metabolites[metabolite.id]
        formula_no_R = metabolite.formula.copy()
        formula_no_R["R"] = 0
        original_no_R = original_metabolite.formula.copy()
        original_no_R["R"] = 0

        if (formula_no_R == original_no_R) and metabolite.charge == original_metabolite.charge:
            charge_constraint = self.charge_symbols[metabolite.id] == metabolite.charge
            constraints.append(z3.And(*[element_symbols[element] == metabolite.formula[element] for element in self.relevant_elements], charge_constraint))
        else:
            for formula, charge in self.assignments[metabolite.id]:
                if (not charge is None):
                    charge_constraint = self.charge_symbols[metabolite.id] == charge
                else:
                    charge_constraint = True
                if "R" in formula:
                    constraints.append(z3.And(*[element_symbols[element] >= formula[element] for element in self.relevant_elements], charge_constraint))
                else:
                    constraints.append(z3.And(*[element_symbols[element] == formula[element] for element in self.relevant_elements], charge_constraint))
        if len(constraints) > 0:
            return z3.Or(constraints)
        else:
            self.unknown_metabolite_ids.add(metabolite.id)
            logging.debug(f"No assignments for {metabolite.id} found.")
            return z3.And(*[element_symbols[element] >= 0 for element in self.relevant_elements])

    def _add_soft_constraints(self):
        """
        Adds optimiziation (soft) constraints to the solver. Specifically adds soft constraints to try to choose formulae based on the following criteria:
        - highest possible atom count

        - if unconstrained:

            1. adhering to an unconstrained database representation

            2. no unnecessary atoms

        """
        for metabolite in self.model_interface.metabolites.values():
            # if there is no constraint on a metabolite formula, we want it to become 0
            if metabolite.id in self.unknown_metabolite_ids:
                for element in self.relevant_elements:
                    self.solver.add_soft(self.metabolite_symbols[metabolite.id][element] == 0)

            # if there is multiple formula, we prefer larger ones
            assignments = set(f for f in self.assignments[metabolite.id])
            if len(assignments) > 1:
                length_sorted_formulae = list(sorted(assignments, key = lambda f: sum(f[0].elements.values())))
                for i in range(len(length_sorted_formulae)):
                    constraints = []
                    for element in length_sorted_formulae[i][0]:
                        if (element == "R"): continue
                        constraints.append(self.metabolite_symbols[metabolite.id][element] == length_sorted_formulae[i][0][element])
                    if not (length_sorted_formulae[i][1] is None):
                        constraints.append(self.charge_symbols[metabolite.id] == length_sorted_formulae[i][1])
                    self.solver.add_soft(z3.And(constraints), weight = 10 * (i + 1))
            