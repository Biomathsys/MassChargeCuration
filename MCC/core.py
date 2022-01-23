from copy import deepcopy
from dataclasses import dataclass, field
from typing import Dict, Set
import logging
import re
import numpy as np


element_re = re.compile(r"([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")
replace_capital_ids = re.compile(r"([A-Z])([A-Z])(\d+)")
replace_deuterium = re.compile(r"D([^ybs]|$)")
replace_tritium = re.compile(r"T([^abceslhmi]|$)")
replace_rest_R = re.compile(r"R\d*([^abefghnu]|$)")
replace_rest_X = re.compile(r"X\d*([^e]|$)")
replace_placeholders = re.compile(r"[*\.]\d*")
remove_1 = re.compile(r"([A-Z][a-z]?)(1)([A-Z]|$)")
remove_isotope_notation = re.compile(r"\[\d+([A-Z][a-z]?)\]")

    

class Formula():

    def __init__(self, formula) -> None:
        self.elements = {}
        if isinstance(formula, str):
            self.elements = Formula._to_dict(formula)
        elif isinstance(formula, dict):
            self.elements = formula
        elif isinstance(formula, Formula):
            self.elements = Formula._to_dict(str(Formula))
        else:
            logging.error(f"Expected argument of type str or dict, not {type(formula)}.")
            raise ValueError
    
    def __repr__(self) -> str:
        return Formula._from_dict(self.elements)

    def __getitem__(self, key):
        return self.elements.get(key, 0)
            
    def __setitem__(self, key, value):
        self.elements[key] = value

    def __iter__(self):
        return iter(self.elements)

    def __hash__(self):
        return hash(repr(self))

    def copy(self):
        return Formula(self)


    @staticmethod
    def _from_dict(formula_dict):
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
            if formula_dict[key] > 0:
                result.append(key)
            if formula_dict[key] > 1:
                result.append(str(formula_dict[key]))
        return "".join(result)

    @staticmethod
    def _to_dict(formula):
        formula = Formula.clean(formula)
        element_dict = {}
        parsed = element_re.findall(formula)
        for (element, count) in parsed:
            if count == "":
                count = 1
            cur_count = element_dict.get(element, 0)
            cur_count += int(count)
            element_dict[element] = cur_count
        return element_dict

    @staticmethod
    def clean(formula):
        """
        Function to clean a given formula of possible artifacts.
        
            - Replaces Deuterium and Tritium with Hydrogen

            - Removes any 1s from atom counts

            - Removes isotope notations

            - Replaces any wildcard symbol with "R"


        Returns:
            This cleaned formula.
        """
        # replace different hydrogen symbols
        formula = replace_deuterium.sub("H", formula)
        formula = replace_tritium.sub("H", formula)
        # replace 1s
        formula = remove_1.sub(r"\1\3", formula)
        # remove isotope notation
        formula = remove_isotope_notation.sub(r"\1", formula)
        # gather rest symbols
        if not ((found_R := replace_rest_R.search(formula)) is None):
            formula = replace_rest_R.sub(r"\1", formula)
        if not ((found_X := replace_rest_X.search(formula)) is None):
            formula = replace_rest_X.sub(r"\1", formula)
        if not ((found_other := replace_placeholders.search(formula)) is None):
            formula = replace_placeholders.sub("", formula)
        if any([found_R, found_X, found_other]):
            formula += "R"
        return formula

    def __eq__(self, __x: object) -> bool:
        if not isinstance(__x, Formula):
            return False 
        else:
            dict1 = self.elements
            dict2 = __x.elements
            same = True
            for key in dict1:
                if dict1[key] != dict2.get(key, 0): same = False
            for key in dict2:
                if dict2[key] != dict1.get(key, 0): same = False
            return same 

@dataclass
class Metabolite():
    id: str
    name: str
    formula: Formula
    charge: int
    SBO: str = field(default = None) 
    reactions: Set = field(repr=False, default_factory = lambda: set())
    cv_terms: Dict = field(repr=False, default_factory = lambda: {})
    notes: Dict = field(repr=False, default_factory = lambda: {})

    def __hash__(self) -> int:
        return hash(self.id)
        
@dataclass
class Reaction():
    id: str
    name: str
    metabolites: Dict = field(default_factory = lambda: {})
    SBO: str = field(default = None) 
    cv_terms: Dict = field(repr=False, default_factory = lambda: {})
    notes: Dict = field(repr=False, default_factory = lambda: {})

    def __hash__(self) -> int:
        return hash(self.id)
 
    def mass_balance(self, assignments = None):
        if assignments is None: assignments = {}
        mass_dict = {}
        for metabolite, coeff in self.metabolites.items():
            if metabolite not in assignments:
                formula = metabolite.formula
            else:
                formula = assignments[metabolite][0]
            for element, count in formula.elements.items():
                mass_dict[element] = mass_dict.get(element, 0) + (count * coeff)
        return mass_dict

    def charge_balance(self, assignments = None):
        charge_balance = 0
        for metabolite, coeff in self.metabolites.items():
            if metabolite not in assignments:
                charge = metabolite.charge
            else:
                charge = assignments[metabolite][1]
            charge_balance += (charge * coeff)
        return charge_balance

    def is_balanced(self, assignments = None, just_charge = False, proton_balance = False):
        if just_charge:
            mass_balance = self.mass_balance(assignments)
            return np.isclose(self.charge_balance(assignments), mass_balance.get("H", 0))
        if not proton_balance:
            mass_balance = self.mass_balance(assignments)
            return all(np.isclose(value, 0) for element, value in mass_balance.items() if element != "H") and np.isclose(self.charge_balance(assignments), mass_balance.get("H", 0))
        else:
            return all(np.isclose(value, 0) for value in self.mass_balance(assignments).values()) and np.isclose(self.charge_balance(assignments), 0)
                
    @property
    def sbo(self):
        return self.SBO

    @sbo.setter
    def sbo(self, value):
        self.SBO = value

    def copy(self):
        return deepcopy(self)

