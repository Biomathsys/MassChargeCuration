import logging

from .ReaderInterface import ReaderInterface


class CobraPyInterface(ReaderInterface):

    def __init__(self, model):
        super().__init__(model)

    def get_model_id(self):
        return self.model.id

    def write_metabolite(self, metabolite_id : str, name : str, formula : str, charge : int, sbo, cv_terms, notes):
        metabolite = self.model.metabolites.get_by_id(metabolite_id)
        metabolite.name = name
        metabolite.formula = formula
        metabolite.charge = charge
        metabolite.annotation = cv_terms
        if sbo is not None:
            metabolite.annotation["sbo"] = f"SBO:{sbo}"
        metabolite.notes = notes

    def write_reaction(self, reaction_id : str, name : str, metabolites, sbo, cv_terms, notes):
        reaction = self.model.reactions.get_by_id(reaction_id)
        reaction.name = name
        to_remove = {}
        for metabolite, coeff in reaction.metabolites.items():
            if metabolite.id not in metabolites:
                to_remove[metabolite.id] = -coeff
        reaction.add_metabolites(to_remove)
        reaction.add_metabolites({meta_id : coeff for meta_id, coeff in metabolites.items()}, combine = False)
        reaction.annotation = cv_terms
        if sbo is not None:
            reaction.annotation["sbo"] = f"SBO:{sbo}"
        reaction.notes = notes

    def write_model(self, filename):
        try:
            import cobra
            cobra.io.write_sbml_model(self.model, filename)
        except ImportError:
            logging.exception("Tried to use cobraPy to write the model but failed to import it.")

    def get_metabolite_ids(self):
        return [metabolite.id for metabolite in self.model.metabolites]

    def get_reaction_ids(self):
        return [reaction.id for reaction in self.model.reactions]

    def get_metabolite_formula_by_id(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id).formula
        
    def get_metabolite_charge_by_id(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id).charge
        
    def get_reaction_metabolite_ids(self, reaction_id):
        return {m.id : coeff for m, coeff in self.model.reactions.get_by_id(reaction_id).metabolites.items()}
       
    def get_metabolite_name(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id).name
      
    def get_reaction_name(self, reaction_id):
        return self.model.reactions.get_by_id(reaction_id).name
     
    def get_metabolite_cv_terms(self, metabolite_id):
        anno = self.model.metabolites.get_by_id(metabolite_id).annotation.copy()
        for key, value in anno.items():
            if type(value) is str:
                anno[key] = [value]
        if "sbo" in anno: del anno["sbo"]
        return anno
    
    def get_reaction_cv_terms(self, reaction_id):
        anno = self.model.reactions.get_by_id(reaction_id).annotation.copy()
        for key, value in anno.items():
            if type(value) is str:
                anno[key] = [value]
        if "sbo" in anno: del anno["sbo"]
        return anno

    def get_metabolite_notes(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id).notes.copy()
   
    def get_reaction_notes(self, reaction_id):
        return self.model.reactions.get_by_id(reaction_id).notes.copy()

    def get_metabolite_sbo(self, metabolite_id):
        anno = self.model.metabolites.get_by_id(metabolite_id).annotation.copy()
        sbo = anno.get("sbo", None)
        if sbo is None: return None
        else: return sbo[4:]
  
    def get_reaction_sbo(self, reaction_id):
        anno = self.model.reactions.get_by_id(reaction_id).annotation.copy()
        sbo = anno.get("sbo", None)
        if sbo is None: return None
        else: return sbo[4:]
 
    def copy(self):
        new_model = self.model.copy()
        return CobraPyInterface(new_model)
        