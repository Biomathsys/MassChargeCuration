

class ReaderInterface():

    def __init__(self, model):
        self.model = model
    
    def get_model_id(self):
        """
        Returns this models id.
        """
        raise NotImplementedError

    #FIXME: current implementation will most likely delete any qualifiers which are not IS
    def write_metabolite(self, metabolite_id : str, name : str, formula : str, charge : int, sbo : int, cv_terms, notes):
        """
        Writes the given metabolite information to the model.

        Args:
            metabolite_id (str): Id of the metabolite for which we write the information.
            name (str): Name of the metabolite.
            formula (str): Formula of the metabolite.
            charge (int): Charge of the metabolite.
            sbo (int): SBO term describing this metabolite.
            cv_terms ({identifier : [terms]}): Dictionary mapping identifiers to lists of terms.
            notes ({note_id : note}): Dictionary mapping note identifiers to the respective note.
        """
        raise NotImplementedError

    def write_reaction(self, reaction_id : str, name : str, metabolites, sbo, cv_terms, notes):
        """
        Writes the given reaction information to the model.

        Args:
            reaction_id (str): Id of the reaction for which we write the information.
            name (str): Name of the reaction.
            metabolites ({metabolite_id (str) : coeff (int)}): Dictionary mapping metabolite id to the respective stoichiometric coefficients.
            sbo (int): SBO term describing this reaction.
            cv_terms ({identifier : [terms]}): Dictionary mapping identifiers to lists of terms.
            notes ({note_id : note}): Dictionary mapping note identifiers to the respective note.
        """
        raise NotImplementedError

    def write_model(self, filename):
        """
        Writes the current model to the given filename.

        Args:
            filename (str): Path to write the model to.
        """
        raise NotImplementedError

    def get_metabolite_ids(self):
        """
        Gets ids for all metabolties in this interfaces model.

        Returns:
            List of metabolite ids.
        """
        raise NotImplementedError

    def get_reaction_ids(self):
        """
        Gets ids for all reactions in this interfaces model.

        Returns:
            List of reaction ids.
        """
        raise NotImplementedError

    def get_metabolite_formula_by_id(self, metabolite_id):
        """
        Gets the formula as a string for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites formula as a string.
        """
        raise NotImplementedError
        
    def get_metabolite_charge_by_id(self, metabolite_id):
        """
        Gets the charge for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites charge.
        """
        raise NotImplementedError
        
    def get_reaction_metabolite_ids(self, reaction_id):
        """
        Gets the ids and stoichiometric coefficents for all metabolites participating in the reaction with the given id.

        Args:
            reaction_id (str): Id of the reaction for which to fetch the information.

        Returns:
            Dictionary mapping metabolite ids (str) to their stoichiometric coefficents.
        """
        raise NotImplementedError
       
    def get_metabolite_name(self, metabolite_id):
        """
        Gets the name for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites name.
        """
        raise NotImplementedError
      
    def get_reaction_name(self, reaction_id):
        raise NotImplementedError
     
    def get_metabolite_cv_terms(self, metabolite_id):
        """
        Gets the cv terms for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites cv terms in form of a dictionary mapping the identifiers to the terms. Only terms with the qualifier IS should be considered.
        """
        raise NotImplementedError
    
    def get_reaction_cv_terms(self, reaction_id):
        """
        Gets the cv terms for the reaction with the given id.
        
        Args:
            reaction_id (str): Id of the reaction for which to fetch the information.

        Returns:
            The reactions cv terms in form of a dictionary mapping the identifiers to the terms. Only terms with the qualifier IS should be considered.
        """
        raise NotImplementedError

    def get_metabolite_notes(self, metabolite_id):
        """
        Gets the notes for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites notes in form of a dictionary mapping the identifiers to the notes.
        """
        raise NotImplementedError
   
    def get_reaction_notes(self, reaction_id):
        """
        Gets the notes for the reaction with the given id.
        
        Args:
            reaction_id (str): Id of the reaction for which to fetch the information.

        Returns:
            The reactions notes in form of a dictionary mapping the identifiers to the notes.
        """
        raise NotImplementedError

    def get_metabolite_sbo(self, metabolite_id):
        """
        Gets the sbo value for the metabolite with the given id.
        
        Args:
            metabolite_id (str): Id of the metabolite for which to fetch the information.

        Returns:
            The metabolites SBO term.
        """
        raise NotImplementedError
  
    def get_reaction_sbo(self, reaction_id):
        """
        Gets the sbo value for the reaction with the given id.
        
        Args:
            reaction_id (str): Id of the reaction for which to fetch the information.

        Returns:
            The reactions SBO term.
        """
        raise NotImplementedError
 
    def copy(self):
        raise NotImplementedError
        