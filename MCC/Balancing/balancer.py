import logging
import z3

class Balancer():
    """
    Abstract Class used to balance a given model. Represents the basic structure of 
    1. setting up the solver => _setup_z3,
    2. generating and adding the assertions => generate_assertions,
    3. simplify assertions => self.simplifying_tactic,
    4. (trying to) find a SAT solution => balance,
    5. resolving UNSAT solutions => resolve_unsat,
    6. assigning from solved model => assign_from solver
    """
    def __init__(self):
        # setting up z3
        self._setup_z3()
        self._prepare_solver()

    def _setup_z3(self):
        """
        Sets up the Z3 solver, generates the simplifying tactics and enables unsat core minimization.
        """
        self.reaction_answer_literals = {}
        self.simplifying_tactic = z3.Then("simplify", "solve-eqs")
        self.solver = z3.Solver()
        self.solver.set(':core.minimize', True)
        z3.set_option("parallel.enable", True) 

    def _prepare_solver(self):
        """
        Prepares the solver by adding all relevant assertions and applying simplifying tactics.
        """
        self.goal = z3.Goal()
        self.generate_assertions()
        self.subgoal = self.simplifying_tactic(self.goal)[0]
        self.solver.add(self.subgoal)

    def balance(self, answer_literals):
        """
        Tries to find a SAT solution assuming the given literals. 
        If a SAT solution is found, the solution is passed to this classes assign_from_solver. 
        Otherwise the UNSAT core is passed to this classes resolve_unsat.

        Args:
            answer_literals([z3.Literal]): Literals which are assumed to be true. The unsat core will be a subset of 
                these literals which cannot be true at the same time.

        """
        if self.solver.check(answer_literals) == z3.sat:
            self.assign_from_solver(self.subgoal.convert_model(self.solver.model()))
        else:
            self.resolve_unsat(self.solver.unsat_core())

    def generate_assertions(self):
        """
        Abstract function - should be implemented by subclasses. This function is intended to add all
        relevant constraints to self.solver.
        """
        raise(NotImplementedError)

    def assign_from_solver(self, sat_model):
        """
        Abstract function - should be implemented by subclasses. This function is intended to take a
        satisfying assignment to all literals and apply them to some actual model.
        """
        raise(NotImplementedError)

    def resolve_unsat(self, unsat_core):
        """
        Function to resolve any unsatisfying solution. Base function only reports the unsat core and returns.
        """
        logging.warn(f"Was unsat with {unsat_core}.")