import cobra
import logging
import sys
import numpy as np
import z3

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

def formula_to_dict(string):
    '''
    Function to convert a formula given as a string to a dictionary representing the formula.
    
    Args:
        string: The formula as a string.
        
    Returns:
        The formula represented as a dict. For example for C48H67N14O38P6S2:
        {'C': 48, 'H': 67, 'N': 14, 'O': 38, 'P': 6, 'S': 2}
        
    '''
    formula = cobra.core.formula.Formula(string)
    formula.parse_composition()
    return formula.elements

def dict_to_formula(formula_dict):
    '''
    Function to turn a formula represented as a dictionary to a formula represented by string.
    
    Args:
        formula_dict: formula represented as a dictionary mapping elements to the number of their occurence.
        
    Returns:
        Given formula as a string.
        
    '''
    sorted_keys = sorted(formula_dict.keys())
    result = []
    for key in sorted_keys:
        result.append(key)
        if formula_dict[key] > 1:
            result.append(str(formula_dict[key]))
    return "".join(result)

def non_hydrogen_balanced(reaction, assignments = None, unknown_is_balanced = False):
    """
    Checks if the given reaction is balanced wrt to any element but hydrogen with the given assignments.
    
    Args:
        reaction (cobrapy.Reaction): Reaction for which to check the balance.
        assignments ({metabolite : formula}): Optional; Dictionary containing the assignments for every metabolite.
            If none is passed, the current formula is used.
        unknown_is_balanced (bool): Deprecated.

    Returns:
        True if the reaction is non-hydrogen balanced else false.
    """

    if assignments is None: assignments = {}
    mass_dict = {}
    for metabolite, coeff in reaction.metabolites.items():
        if metabolite.id not in assignments:
            mass = metabolite.elements
        elif assignments[metabolite.id] == "any": 
                mass = {}
        else:
            mass = formula_to_dict(assignments[metabolite.id])
        for atom, count in mass.items():
            if atom == "H": continue
            mass_dict[atom] = mass_dict.get(atom, 0) + (count * coeff)
    return all(val == 0 for val in mass_dict.values())

def is_cH_balanced(reaction, assignments = None, unknown_is_balanced = False):
    """
    Function the check whether or not a reaction is balanced wrt its hydrogen atoms and charges.

    Args:
        reaction (cobrapy.Reaction): Reaction for which to check the balance.
        assignments ({metabolite : (H count, charge)}): Optional; Dictionary containing the assignments for every metabolite.
            If none is passed, the current values is used.
        unknown_is_balanced (bool): If True, if any charge is unknown the reaction is considered balanced.

    Returns:
        True if the reaction is cH balanced, else False.
    """
    h_sum = 0
    charge_sum = 0
    for metabolite, coeff in reaction.metabolites.items():
        if metabolite.charge is None:
            return unknown_is_balanced
        if assignments is None:
            h_sum += metabolite.elements.get("H", 0) * coeff
            charge_sum += metabolite.charge * coeff
        else:
            if ((assignments[metabolite.id][0] is None) or (assignments[metabolite.id][1] is None) or (assignments[metabolite.id] == "any")):
                return unknown_is_balanced
            h_sum += assignments[metabolite.id][0] * coeff
            charge_sum += int(assignments[metabolite.id][1]) * coeff
        
    return h_sum == charge_sum # protons can be added to balance
        

def get_pseudo_reactions(model):
    """
    Gives all pseudo reaction for a given model.

    Args:
        model (cobrapy.Model): Model for which to determine the pseudoreactions.

    Returns:
        List of reaction ids which are pseudo reactions.
            => [pseudoreaction_ids]
    """
    pseudo_reactions = set()
    for reaction in model.reactions:
        if len(reaction.products) == 0 or len(reaction.reactants) == 0:
            pseudo_reactions.add(reaction.id)
        if ('sbo' in reaction.annotation and reaction.annotation['sbo'] == 'SBO:0000629') or ("growth" in reaction.id.lower()):
            pseudo_reactions.add(reaction.id)
    return pseudo_reactions


def lowest_charge_preferred(balance_dict):
    """
    For a given dictionary mapping charge/hydrogen differences determines the least absolute charge
    representative for every charge/hydrogen difference.

    Args:
        balance_dict ({difference : [(formula, charge)]})

    Returns:
        Dictionary mapping each charge/hydrogen difference to one representative.
            => {difference : (formula, charge)}
    """
    new_dict = {}
    for balance, assignments in balance_dict.items():
        new_dict[balance] = min(list(assignments), key = lambda x : 1e9 if x[1] is None else np.abs(x[1]))
    return new_dict
        
def adjust_proton_count(reaction):
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
        h_balance += metabolite.elements.get("H", 0) * coeff
    charge_balance = np.round(charge_balance)
    h_balance = np.round(h_balance)
    if charge_balance == h_balance:
        if charge_balance > 10:
            logging.info(f"added {charge_balance} protons to reaction {reaction.id}")
        h_id = None
        for metabolite in reaction.metabolites:
            if metabolite.id.startswith("h_"):
                h_id = metabolite.id
                break
        if h_id is None:
            for metabolite in reaction.metabolites:
                h_id = "h_{}".format(metabolite.id[-1])
                try:
                    reaction.subtract_metabolites({h_id : 0})
                    break
                except KeyError as e:
                    h_id = None
        if h_id is None:
            return 0
        reaction.subtract_metabolites({h_id : charge_balance})
        old_notes = reaction.notes
        old_notes["Inferred"] = "Protons added to balance equation"
        reaction.notes = old_notes
        return charge_balance
    return 0


def calculate_balance(reaction, assignments = None):
    """
    Calculates the mass and charge balance for a given reaction.

    Args:
        reaction (cobrapy.Reaction): Reaction for which to determine the balance.
        assignments ({metabolite : (formula, charge)}): Optional; Assignments to use for every metabolite.
            If none is given, defaults to the current values for the metabolites.

    Returns:
        Dictionary containing mass and charge difference for this reaction.
            => {"mass": {atom : difference}, "charge" : charge_difference}
    """
    charge_balance = 0
    mass_dict = {}
    for metabolite, coeff in reaction.metabolites.items():
        if assignments:
            mass = formula_to_dict(assignments[metabolite.id][0])
            charge = assignments[metabolite.id][1]
        else:
            mass = metabolite.elements
            charge = metabolite.charge
        charge_balance += charge * coeff
        for atom, count in mass.items():
            mass_dict[atom] = mass_dict.get(atom, 0) + (count * coeff)
    return {"mass" : mass_dict, "charge" : charge_balance}

def is_balanced(reaction, assignments = None, unknown_is_balanced = False):
    """
    Returns True if the given reaction is balanced, otherwise False.

    Args:
        reaction (cobrapy.Reaction): Reaction for which to determine whether it is balanced.
        assignments ({metabolite : (formula, charge)}): Optional; Assignments which to use for the metabolites.
            If none is given, the current values are used.
        unknown_is_balanced (bool): Optional; If True if any metabolites charge or mass is unknown the reaction is assumed to be balanced.
        TODO: unkown is balanced makes currently no sense
    """
    if unknown_is_balanced:
        raise NotImplementedError
    balance = calculate_balance(reaction, assignments)
    return all(np.isclose(val, 0) for val in balance["mass"].values()) and (balance["charge"] == 0)
    



def score_nonH_assignment(model, assignments = None):
    """
    Function to count the number of non-hydrogen imbalanced reactions in a model.

    Args:
        model (cobrapy.Model): Model which to evaluate.
        assignments ({metabolite : (formula, charge)}): Assignments to use.

    Returns:
        Number of non-hydrogen unbalanced reactions in the given model.
    """
    score = 0
    pseudo_reactions = get_pseudo_reactions(model)
    for reaction in model.reactions:
        if reaction in pseudo_reactions: continue
        if not non_hydrogen_balanced(reaction, assignments):
            score += 1
        
    return score

def apply_assignment(assignment, model):
    """
    Function to write a set of given assignments to the given model.

    Args:
        assignments ({metabolite : (formula, charge)}): Assignments to use.
        model (cobrapy.Model): Model to which we want to apply the assignments to.  

    """
    for metabolite_id, a in assignment.items():
        m = model.metabolites.get_by_id(metabolite_id)
        e_dict = m.elements
        e_dict.update(formula_to_dict(a))
        m.elements = e_dict 

def same_formula(f1, f2, ignore_rest = False):
    """
    Function to evaluate if two mass formula strings represent the same atom counts.

    Args:
        f1 (str): First mass formula string to compare.
        f2 (str): Second mass formula string to compare.
        ignore_rest (bool): Optional; If True, wildcards are ignored in the comparison.

    Returns:
        True if the two strings represent the same atom counts, otherwise False.
    """
    if (f1 is None) or (f2 is None):
        return f1 == f2
    dict1 = formula_to_dict(f1) 
    dict2 = formula_to_dict(f2)
    same = True
    for key in dict1:
        if key == "R" and ignore_rest: continue
        if dict1[key] != dict2.get(key, 0): same = False
    for key in dict2:
        if key == "R" and ignore_rest: continue
        if dict2[key] != dict1.get(key, 0): same = False
    return same 

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
    subset_dict = formula_to_dict(subset_f) 
    superset_dict = formula_to_dict(superset_f)
    same = True
    for key in subset_dict:
        if key == "R": continue
        if subset_dict[key] > superset_dict.get(key, 0): same = False
    for key in superset_dict:
        if key == "R": continue
        if superset_dict[key] < subset_dict.get(key, 0): same = False
    return same 

def get_proton_count(reaction):
    """
    Function to determine the proton balance in a reaction.

    Args:
        reaction (cobrapy.Reaction): Reaction for which we want to determine the proton balance.
    """
    total_count = 0
    for metabolite, count in reaction.metabolites.items():
        if metabolite.formula == "H":
            total_count += np.absolute(count)
    return total_count


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
        for metabolite, coeff in reaction.metabolites.items():
            coeff *= factor
            if not coeff.is_integer(): 
                non_int_found = True
                factor = 1/np.absolute(coeff)
                break
    return factor