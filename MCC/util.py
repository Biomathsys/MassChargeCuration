import logging
import sys
from turtle import pos
import numpy as np
import re
import z3
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


logging_is_setup = False
def logging_setup(loglevel = "info"): 
    """
    Sets up logging, deprecated.
    """
    global logging_is_setup
    if logging_is_setup: return
    # setting logging stuff
    loglevel = loglevel.lower()
    logging_map = {"debug" : logging.DEBUG, "info" : logging.INFO, "warning" : logging.WARNING, "error" : logging.ERROR, "critical" : logging.CRITICAL}
    level = logging_map.get(loglevel, logging.INFO)
    logging.basicConfig(format='%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s', level = level)
    formatter = logging.Formatter('%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s')
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(formatter)
    logging_is_setup = True


def get_assertion_leafs(assertion):
    """
    Finds the leafs of a z3 assertion.
    """
    if assertion == False:
        print("found false in assertions")
    if len(assertion.children()) == 0:
        if type(assertion) != z3.z3.IntNumRef:
            return [assertion]
        else:
            return []
    else:
        assertions = []
        for child in assertion.children():
            assertions.extend(get_assertion_leafs(child))
        return assertions

def get_pseudo_reactions(model):
    """
    Gives all pseudo reactions ids for a given model.

    Args:
        model (cobrapy.Model): Model for which to determine the pseudoreactions.

    Returns:
        List of reaction ids which are pseudo reactions.
            => [pseudoreaction_ids]
    """
    pass
def adjust_proton_count(reaction, model_interface):
    """
    Adds protons to a given reaction to fully balance it out.
    Tries to use protons which are already present in the reaction. If none can be found,
    protons of a compartment of an aribitrary reactants are chosen.

    Args:
        reaction (cobrapy.Reaction): Reaction which we want to add protons to. 
    
    Returns:
        Old charge balance.
    """
    charge_balance = 0
    h_balance = 0
    for metabolite, coeff in reaction.metabolites.items():
        charge_balance += metabolite.charge * coeff
        h_balance += metabolite.formula["H"] * coeff
    charge_balance = np.round(charge_balance)
    h_balance = np.round(h_balance)
    if not charge_balance == h_balance:
            logging.warn(f"Adding {h_balance} protons to charge unbalanced reaction {reaction.id} to maintain stoichiometrtic consistency.")
    elif charge_balance > 10:
        logging.info(f"added {charge_balance} protons to reaction {reaction.id}")
    h_id = None
    for metabolite in reaction.metabolites:
        if str(metabolite.formula) == "H":
            h_id = metabolite
            break
    if h_id is None:
        possible_h = set()
        for metabolite in model_interface.metabolites.values():
            if str(metabolite.formula) == "H":
                possible_h.add(metabolite)
        for metabolite in reaction.metabolites:
            for hydrogen in possible_h:
                if hydrogen.id[-1] == metabolite.id[-1]:
                    h_id = hydrogen
    if h_id is None:
        if len(possible_h) > 0:
            h_id = possible_h.pop()
            logging.warn(f"Could not find appropriate hydrogen to balance reaction {reaction.id}. Chose {h_id}.")
        else:
            logging.error(f"Could not find appropriate hydrogen to balance reaction {reaction.id}.")
            return 0
    if not h_id.id.lower().startswith("h"):
        logging.warn(f"Found {h_id} as proton metabolite. If this is not intended, try to fix its formula using the fixed_assignments argument.")
    reaction.metabolites[h_id] =  reaction.metabolites.get(h_id, 0) - h_balance
    # since note appending does not work / documentation is unclear, we use a workaround
    reaction.notes["Inferred"] =  f"{h_balance} protons added to {'reactants' if h_balance > 0 else 'products'} to balance equation."
    return h_balance

def subset_formula(subset_f, superset_f):
    """
    Function to determine if one mass formula string represents a subset of atoms of the other
    mas formula string.

    Args:
        subset_f (str): Mass formula string for which we want to determine if it is a subset of the other.
        superset_f (str): Mass formula string for which we want to determine if it is a superset of the other.


    Return:
        True if the atoms represented by subset_f are a subset of the atoms represented by superset_f, otherwise False.
    """
    if (subset_f is None) or (superset_f is None):
        return subset_f == superset_f
    same = True
    for key in subset_f:
        if key == "R": continue
        if subset_f[key] > superset_f[key]: same = False
    for key in superset_f:
        if key == "R": continue
        if superset_f[key] < subset_f[key]: same = False
    return same 

def get_integer_coefficients(reaction):
    """
    Function to try to find integer coefficents for a given reaction.
    Works by repeatedly multiplying with the reciprocal of a non-integer coefficent. 
    If this does not work after 3 times, tries to multipy by 10, if that doesnt work, gives up.

    Args:
        reaction (cobrapy.Reaction): Reaction for which we want to determine the interger coefficients.

    Returns:
        Factor with which the stochiometric coeficients need to be multiplied to gain integer coefficients.
    """
    non_int_found = None
    factor = 1
    counter = 0
    while not (non_int_found is False):
        counter += 1
        if counter == 3:
            factor = 10
        if counter == 5:
            logging.warning(f"Could not compute integer coefficients for reaction {reaction.id}... skipping it for now")
            return None
        non_int_found = False
        for _, coeff in reaction.metabolites.items():
            coeff *= factor
            if not coeff.is_integer(): 
                non_int_found = True
                factor = 1/np.absolute(coeff)
                break
    return factor
