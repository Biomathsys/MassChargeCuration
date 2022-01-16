from .fullBalancer import FullBalancer
import time
import logging
from z3 import sat
from ..util import is_balanced, get_sbml_metabolites

class SatCore(FullBalancer):
    """
    Class to resolve unsat cores. The strategy is to remove any unsat core from the used literals until the model is sat.
    Then for each unsat core, we try to add as many literals back to the model without it becoming unsat again.
    The literals are added back based on the following heuristic:

        - Smaller unsat cores are considered first.

        - Literals within an unsat cores are scored, based on the reactions which they represent.

            > A reaction is scored higher the closer it is to a fixed metabolite assignment. If there is no fixed assignment, the reactions will be scored
                heuristically based on the metabolite assignments which would balance them and the score of each of these assignments.

            > The score for a metabolite assignment is based on the number of reactions which can be balanced using this assignment.
    """
    
    def __init__(self, model, data_collector, fixed_assignments = None, fixed_reactions = None, *args, **kwargs):
        self.fixed_reactions = set() if fixed_reactions is None else fixed_reactions
        super().__init__(model, data_collector, fixed_assignments, *args, **kwargs)
        self.balance_function = is_balanced
        self.reaction_scores = self.score_reactions()
        self.unsat_cores = []

    def resolve_unsat(self, unsat_core):
        """
        Function to resolve an unsat core. Gets called when the first balancing fails.
        The unsat core is solved via the following strategy:
        Remove any unsat core from the used literals until the model is sat.
        Then for each unsat core, we try to add as many literals back to the model without it becoming unsat again.
        The literals are added back based on the following heuristic:

        - Smaller unsat cores are considered first.

        - Literals within an unsat cores are scored, based on the reactions which they represent.

            > A reaction is scored higher the closer it is to a fixed metabolite assignment. If there is no fixed assignment, the reactions will be scored
                heuristically based on the metabolite assignments which would balance them and the score of each of these assignments.

            > The score for a metabolite assignment is based on the number of reactions which can be balanced using this assignment.

        Returns:
            z3.sat if unsat core could be resolved, z3.unsat otherwise (ideally shouldn't happen).
        """

        # first we collect all unsat cores and determine the sat core
        unsat_cores = []
        used_literal_ids = set(self.answer_literals.keys()).difference(self.unbalancable_reactions)
        sat_core = set([self.answer_literals[id] for id in used_literal_ids])
        t = time.process_time()
        while len(unsat_core) > 0:
            unsat_cores.append(unsat_core)
            sat_core = sat_core.difference(unsat_core)
            self.solver.check(sat_core)
            unsat_core = self.solver.unsat_core()
        self.unsat_cores = sorted(unsat_cores, key=len)
        logging.debug(f"[{time.process_time() - t:.3f} s] unsat cores were: {unsat_cores}")
        t = time.process_time()

        # trying to recover as many reactions as possible from unsat cores
        for unsat_core in self.unsat_cores:
            if len(unsat_core) == 1: 
                reaction_id = get_rid(unsat_core[0])
                if reaction_id in self.fixed_reactions:
                    logging.error(f"Could not balance fixed reaction {reaction_id}.")
                self.unbalancable_reactions.add(reaction_id)
                self.reaction_reasons[reaction_id] = [get_rid(literal) for literal in unsat_core]
                continue # reaction is unbalancable by itself
            else:
                # first determine fixed scores
                fixed_scores = {}
                reactions_ids = [get_rid(literal) for literal in unsat_core]
                for answer_literal in unsat_core:
                    reaction_id = get_rid(answer_literal)
                    # fixed reactions get the highest scores
                    if reaction_id in self.fixed_reactions:
                        fixed_scores[reaction_id] = 100 * len(unsat_core)
                        continue
                    # otherwise we determine the scores based on fixed metabolites
                    for metabolite_id in get_sbml_metabolites(self.model.getReaction(reaction_id)):
                        if metabolite_id in self.fixed_assignments:
                            distances = self._get_reaction_distances(reaction_id, reactions_ids)
                            for answer_literal in unsat_core:
                                score = fixed_scores.get(answer_literal, 0) + 1/distances.get(get_rid(answer_literal), len(unsat_core) + 1)
                                fixed_scores[answer_literal] = score
                            pass
                grouped_literals = {}

                # finally we bucket by score and then
                for answer_literal in unsat_core:
                    score = fixed_scores.get(answer_literal, 0)
                    cur_group = grouped_literals.get(score, set())
                    cur_group.add(answer_literal)
                    grouped_literals[score] = cur_group

                # sort each bucket by the reaction qualities
                for score in sorted(grouped_literals.keys(), reverse=True):
                    cur_group = grouped_literals[score]
                    for answer_literal in sorted(cur_group, key = self._get_reaction_score, reverse=True):
                        sat_core.add(answer_literal)
                        if self.solver.check(sat_core) != sat:
                            reaction_id = get_rid(answer_literal)
                            if reaction_id in self.fixed_reactions:
                                logging.error(f"Could not balance fixed reaction {reaction_id}.")
                            sat_core.remove(answer_literal)
                            self.unbalancable_reactions.add(reaction_id)
                            self.reaction_reasons[reaction_id] = [get_rid(literal) for literal in unsat_core]
        return self.balance()

    def _get_reaction_distances(self, reaction_id, reactions_ids):
        """
        Function to get the closest distance between two reactions using BFS. 
        Used to determine the score when based on fixed assignments.

        Args:
            reaction_id (str): Id of the reaction from which to determine the distances. 
            reactions_ids ([str]): Ids of all reactions for which to determine the distance.

        Returns:
            Dictionary mapping the reactions_ids of reactions to their respective distance from the reaction with reaction_id.
            => {reaction_id : int}
        """
        distances = {}
        covered = set()
        worklist = [(reaction_id, 1)]
        while len(worklist) != 0:
            reaction_id, distance = worklist.pop(0)
            if reaction_id in covered:
                continue
            distances[reaction_id] = distance
            covered.add(reaction_id)
            for metabolite in self.model.reactions.get_by_id(reaction_id).metabolites:
                for reaction in metabolite.reactions:
                    if (not reaction.id in covered) and (reaction.id in reactions_ids):
                        worklist.append((reaction.id, distance  + 1))
        return distances

            


    def _get_reaction_score(self, answer_literal):
        """
        Function to return the score for the given answer_literal, based on the corresponding reaction.

        Args:
            answer_literal (z3.Literal): Answer literal for which to determine the score.
        """
        return self.reaction_scores[get_rid(answer_literal)]

    def score_reactions(self):
        """
        Function to score all reactions.
        Reactions are scored based on: 

            - A reaction is scored higher the closer it is to a fixed metabolite assignment. If there is no fixed assignment, the reactions will be scored
                heuristically based on the metabolite assignments which would balance them and the score of each of these assignments.

            - The score for a metabolite assignment is based on the number of reactions which can be balanced using this assignment.

        Returns:
            Dictionary mapping reaction ids to scores.
            => {reaction_id : float}
        """
        # first round collecting the assignment votes
        assignment_votes = {metabolite.id[2:] : {} for metabolite in self.model.getListOfSpecies()}
        balanced_combinations = {}
        for reaction in self.model.getListOfReactions():
            balanced_combinations[reaction.id] = self._get_balanced_combinations(reaction, self.assignments, self.balance_function)
            for assignments in balanced_combinations[reaction.id]:
                for metabolite_id, assignment in assignments.items():
                    assignment_votes[metabolite_id][assignment] = assignment_votes[metabolite_id].get(assignment, 0) + 1

        # normalize votes
        for metabolite_id in assignment_votes:
            summed_votes = sum(assignment_votes[metabolite_id].values())
            for assignment, votes in assignment_votes[metabolite_id].items():
                assignment_votes[metabolite_id][assignment] = votes/summed_votes

        reaction_scores = {}
        # score reactions based on their metabolites, the combination with the best score determines the score of the reaction
        for reaction in self.model.getListOfReactions():
            scores = []
            for combination in balanced_combinations[reaction.id]:
                combination_score = 0
                for metabolite_id, assignment in combination.items():
                    if (assignment is None) or (assignment == 'any'): continue
                    else:
                        combination_score += 2**assignment_votes[metabolite_id][assignment]
                scores.append((combination_score/len(combination)) if combination else 0)
            if len(scores) == 0: reaction_scores[reaction.id] = 0
            else:
                reaction_scores[reaction.id] = max(scores)
        return reaction_scores


def get_rid(literal):
    """
    Returns the corresponding reaction for a given literal.
    """
    return literal.decl().name()