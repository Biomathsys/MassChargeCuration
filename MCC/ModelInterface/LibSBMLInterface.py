import libsbml
import re
import logging

identifiers_pattern = re.compile(r"https?://identifiers.org/(.+?)[:/](.+)")
identifiers_url = "https://identifiers.org"

class LibSBMLInterface():
    def __init__(self, model):
        self.model = model

    def write_metabolite(self, metabolite_id : str, name : str, formula : str, charge : int, sbo, cv_terms, notes):
        metabolite = self.model.getSpecies(metabolite_id)
        self._set_sbml_name(metabolite, name)
        plugin = self._get_fbc_plugin(metabolite)
        plugin.setChemicalFormula(formula)
        plugin.setCharge(charge)
        self._set_sbml_SBO(metabolite, sbo)
        self._set_sbml_cv_terms(metabolite, cv_terms)
        self._set_sbml_notes(metabolite, notes)

    def write_reaction(self, reaction_id : str, name : str, metabolites, sbo, cv_terms, notes):
        reaction = self.model.getReaction(reaction_id)
        self._set_sbml_name(reaction, name)

        # setting metabolites
        new_reactants = {metabolite_id : - count for metabolite_id, count in metabolites.items() if count < 0}
        new_products = {metabolite_id : count for metabolite_id, count in metabolites.items() if count > 0}
        for metabolite in reaction.getListOfReactants():
            if metabolite.species in new_reactants:
                metabolite.setStoichiometry(new_reactants[metabolite.species])
                new_reactants.pop(metabolite.species)
            else:
                reaction.removeReactant(metabolite.species)
                
        for metabolite in reaction.getListOfProducts():
            if metabolite.species in new_products:
                metabolite.setStoichiometry(new_products[metabolite.species])
                new_products.pop(metabolite.species)
            else:
                reaction.removeProduct(metabolite.species)
        for metabolite_id, count in new_reactants.items():
            reaction.addReactant(self.model.getSpecies(metabolite_id), count)
        for metabolite_id, count in new_products.items():
            reaction.addProduct(self.model.getSpecies(metabolite_id), count)
        # setting other values
        self._set_sbml_SBO(reaction, sbo)
        self._set_sbml_cv_terms(reaction, cv_terms)
        self._set_sbml_notes(reaction, notes)

    def write_model(self, filename):
        libsbml.writeSBMLToFile(self.model.getSBMLDocument(), filename)


    def get_metabolite_ids(self):
        return [species.id for species in self.model.getListOfSpecies()]

    def get_reaction_ids(self):
        return [reaction.id for reaction in self.model.getListOfReactions()]

    def get_metabolite_formula_by_id(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        plugin = self._get_fbc_plugin(metabolite)
        return plugin.chemical_formula

    def get_metabolite_charge_by_id(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        plugin = self._get_fbc_plugin(metabolite)
        return plugin.charge

    def get_reaction_metabolite_ids(self, reaction_id):
        reaction = self.model.getReaction(reaction_id)
        metabolites = {}
        for metabolite in reaction.getListOfProducts():
            metabolites[metabolite.species] = metabolite.stoichiometry
        for metabolite in reaction.getListOfReactants():
            metabolites[metabolite.species] = -metabolite.stoichiometry
        return metabolites

    def get_metabolite_name(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        return metabolite.name

    def get_reaction_name(self, reaction_id):
        reaction = self.model.getReaction(reaction_id)
        return reaction.name


    def get_metabolite_cv_terms(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        return self._get_sbml_cv_terms(metabolite)
        
    def get_reaction_cv_terms(self, reaction_id):
        reaction = self.model.getReaction(reaction_id)
        return self._get_sbml_cv_terms(reaction)

    def get_metabolite_notes(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        return self._get_sbml_notes(metabolite)
        
    def get_reaction_notes(self, reaction_id):
        reaction = self.model.getReaction(reaction_id)
        return self._get_sbml_notes(reaction)

    def get_metabolite_sbo(self, metabolite_id):
        metabolite = self.model.getSpecies(metabolite_id)
        return self._get_sbml_SBO(metabolite)

    def get_reaction_sbo(self, reaction_id):
        reaction = self.model.getReaction(reaction_id)
        return self._get_sbml_SBO(reaction)

    def _get_sbml_SBO(self, sbml_object):
        return sbml_object.getSBOTerm()

    def _set_sbml_SBO(self, sbml_object, term):
       sbml_object.setSBOTerm(term) 

    def _set_sbml_name(self, sbml_object, name):
        sbml_object.setName(name)

    def _set_sbml_notes(self, sbml_object, notes):
        notes_str = "<html xmlns='http://www.w3.org/1999/xhtml'>\n<p>" + "</p>\n<p>".join(f"{identifier} : {note}" for identifier, note in notes.items()) + "</p>\n</html>"
        sbml_object.setNotes(notes_str)

    def _get_sbml_notes(self, sbml_object):
        notes_dict = {}
        notes = sbml_object.getNotes()
        for i in range(notes.getNumChildren()):
            html = notes.getChild(i)
            if html.getName() == "html":
                for i in range(html.getNumChildren()):
                    note = html.getChild(i)
                    if (note.getName() == "p") and (note.getNumChildren() >= 1) and (note.getChild(0).isText()):
                        base_text = note.getChild(0).getCharacters()
                        splitted = base_text.split(":")
                        identifier = splitted[0]
                        value = "".join(splitted[1:]).strip()
                        notes_dict[identifier] = value
        return notes_dict

    def _set_sbml_cv_terms(self, sbml_object, cv_terms):
        resources = set([f"{identifiers_url}/{db_id}/{identifier}"for db_id, identifiers in cv_terms.items() for identifier in identifiers])
        cv_terms = sbml_object.getCVTerms()
        if cv_terms is None:
            cv_term = libsbml.CVTerm()
            cv_term.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
            cv_term.setBiologicalQualifierType(libsbml.BQB_IS)
            for resource in resources:
                cv_term.addResource(resource)
            sbml_object.addCVTerm(cv_term)
            return
        #FIXME: Do we need to treat nested CV Terms?
        found = False
        for cv_term in cv_terms:
            # we only care for is relationships
            if cv_term.getBiologicalQualifierType() == libsbml.BQB_IS:
                if found == True:
                    logging.error(f"Found multiple CV-Terms with qualifiertype 'is' this is not implemented.")
                # identify those to remove and those to add
                to_remove = set()
                for resource_id in range(cv_term.getNumResources()):
                    resource = cv_term.getResourceURI(resource_id)
                    if not resource in resources:
                        to_remove.add(resource) # keep track to remove the resource
                    else: 
                        resources.remove(resource) # the resource is already present
                for missing_resource in resources:
                    cv_term.addResource(missing_resource)
                for resource in to_remove:
                    cv_term.removeResource(resource)

    def _get_sbml_cv_terms(self, sbml_object):
        terms = {}
        cv_terms = sbml_object.getCVTerms()
        if cv_terms is None:
            return {}
        #FIXME: Do we need to treat nested CV Terms?
        for cv_term in cv_terms:
            # we only care for is relationships
            if cv_term.getBiologicalQualifierType() == libsbml.BQB_IS:
                for resource_id in range(cv_term.getNumResources()):
                    matches = identifiers_pattern.match(cv_term.getResourceURI(resource_id))
                    db_id, identifier = matches.group(1), matches.group(2)
                    cur_identifiers = terms.get(db_id, [])
                    cur_identifiers.append(identifier)
                    terms[db_id] = cur_identifiers
        return terms

    def _get_sbml_annotations(self, sbml_object):
        identifiers = {}
        annotations_node = sbml_object.getAnnotation() 
        if annotations_node is None: return identifiers
        for i in range(annotations_node.getNumChildren()):
            child_rdf = annotations_node.getChild(i)
            if child_rdf.getName() == "RDF":
                for i in range(child_rdf.getNumChildren()):
                    child_dscr = child_rdf.getChild(i)
                    if child_dscr.getName() == "Description":
                        for i in range(child_dscr.getNumChildren()):
                            child_is = child_dscr.getChild(i)
                            if child_is.getName() == "is":
                                for i in range(child_is.getNumChildren()):
                                    child_bag = child_is.getChild(i)
                                    if child_bag.getName() == "Bag":
                                        for i in range(child_bag.getNumChildren()):
                                            child_anno = child_bag.getChild(i)
                                            attributes = child_anno.getAttributes()
                                            for i in range(attributes.getNumAttributes()):
                                                if attributes.getName(i) == "resource":
                                                    attribute = attributes.getValue(i)
                                                    splitted = attribute.split("/")
                                                    identifier = splitted[-1]
                                                    db_id = splitted[-2]
                                                    identifiers[db_id] = identifier
        return identifiers


    def _get_fbc_plugin(self, sbml_object):
        for i in range(sbml_object.num_plugins):
            plugin = sbml_object.getPlugin(i)
            if plugin.package_name == "fbc":
                return plugin

    def copy(self):
        new_model = self.model.clone()
        return LibSBMLInterface(new_model)