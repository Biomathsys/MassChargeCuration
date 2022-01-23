from MCC.util import get_fbc_plugin, get_sbml_metabolites

from ..ModelInterface.ModelInterface import ModelInterface
from .satCore import SatCore
from .formulaOptimizer import FormulaOptimizer
from .adherenceOptimizer import AdherenceOptimizer
from ..DataCollection.DataCollection import DataCollector
from ..util import  adjust_proton_count, formula_to_dict, get_pseudo_reactions, is_balanced, same_formula, clean_formula, get_fbc_plugin, get_sbml_metabolites
import logging
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

remove_H = re.compile(r"H\d*([^gfe]|$)") # remember to not also remove H from Hg

class MassChargeCuration:
    
    def __init__(self, model, data_collector = None, data_path = "/data", fixed_assignments = None, fixed_reactions = None, run_optimization = True, **kw):
        total_time = time.process_time()
        self.model_interface = ModelInterface(model)
        self.pseudo_reactions = self.model_interface.get_pseudo_reactions()
        self.original_model_interface = self.model_interface.copy()
        self.data_collector = DataCollector(model, data_path, **kw) if data_collector is None else data_collector
        self.fixed_assignments = fixed_assignments
        self.proton_adjusted_reactions = {}

        # finding unbalancable reactions
        t = time.process_time()
        self.balancer = SatCore(self.model_interface, self.data_collector, fixed_assignments, fixed_reactions)
        logging.info(f"[{time.process_time() - t:.3f} s] Finished constructing model.")
        self.balancer.balance()
        logging.info(f"[{time.process_time() - t:.3f} s] Finished balancibility check. {len(self.balancer.unbalancable_reactions.difference(self.pseudo_reactions))} non-pseudo reactions were unbalancable.")
        for reaction in self.model_interface.reactions:
            if not reaction.is_balanced():
                self.balancer.unbalancable_reactions.add(reaction.id)
        logging.info(f"[{time.process_time() - t:.3f} s] Finished balancibility check. {len(self.balancer.unbalancable_reactions.difference(self.pseudo_reactions))} non-pseudo reactions were unbalancable.")

        if run_optimization:
            # optimize formula selections
            t = time.process_time()
            self.optimizer = AdherenceOptimizer(self.balancer, self.original_model_interface)
            self.optimizer.balance()
            logging.info(f"[{time.process_time() - t:.3f} s] Finished adherence optimization.")

            # optimize formula selections
            t = time.process_time()
            self.optimizer = FormulaOptimizer(self.balancer, self.original_model_interface)
            self.optimizer.balance()
            logging.info(f"[{time.process_time() - t:.3f} s] Finished formula optimization.")

        # add wildcards back in for unconstrained formulae
        self.reintroduce_wildcards()

        # try to take representations from original model
        self.fit_to_original()

        # remove any 0 entries from metabolite formulae
        self.clear_formulae()

        # adjusting protons for reactions
        self.adjust_protons()
        self.total_time = time.process_time() - total_time

    def reintroduce_wildcards(self):
        # get any unkown metabolites
        wildcard_metabolites = set()
        for metabolite in self.model_interface.metabolites():
            assignments = self.assignments[metabolite.id]
            if any(("R" in formula) for formula, charge in assignments):
                matched_clean = any(metabolite.formula == formula for formula, charge in assignments)
                if not matched_clean:
                    wildcard_metabolites.add(metabolite)
            if (len(assignments) == 0):
                original_metabolite = self.original_model_interface[metabolite.id]
                if (not metabolite.formula == original_metabolite.formula) or (not metabolite.charge == original_metabolite.charge):
                        wildcard_metabolites.add(metabolite)
    
        # check for inferred formulae
        unchecked_metabolites = set(wildcard_metabolites)
        inferred = set()
        while len(unchecked_metabolites):
            unchecked_metabolite = unchecked_metabolites.pop()
            reactions_counts = {}
            for reaction in unchecked_metabolite.reactions:
                reactions_counts[reaction] = len([m for m in reaction.metabolites.values() if (m in wildcard_metabolites)])
            if any(count == 1 for count in reactions_counts.values()):
                wildcard_metabolites.remove(unchecked_metabolite)
                inferred.add(unchecked_metabolite)
                for reaction, count in reactions_counts.items():
                    if count > 1:
                        unchecked_metabolites.update(wildcard_metabolites.intersection(reaction.metabolites.values()))
                        unchecked_metabolites.discard(unchecked_metabolite)
        # for inferred formulae we add ECO terms
        # FIXME: reactivate
        """
        for inferred_metabolite in inferred:
            inferred_metabolite.annotation["eco"] = "ECO:0000305"
            inferred_metabolite.notes["inferrence"] = "Inferred formulae"
        """
        # add wildcard symbol
        for metabolite in wildcard_metabolites:
            assignments = self.assignments[metabolite.id]
            matched_clean = False
            for formula, charge in assignments:
                formula_no_R = formula.copy()
                formula_no_R["R"] = 0
                metabolite_formula_no_R = metabolite.formula.copy()
                metabolite_formula_no_R["R"] = 0
                if formula_no_R == metabolite_formula_no_R:
                    metabolite.formula = formula
                    matched_clean = True
            if not matched_clean:
                formula["R"] = 1
                metabolite.formula = formula
                

    def fit_to_original(self):
        assignments = self.balancer._calculate_cH_equivalents(reduce = False)
        for metabolite in self.model_interface.metabolites.values():
            original_metabolite = self.original_model_interface.metabolites[metabolite.id]
            nonH_formula = metabolite.formula.copy()
            nonH_formula["H"] = 0
            nonH_original = original_metabolite.formula.copy()
            nonH_original["H"] = 0
            if nonH_formula == nonH_original:
                diff_seperated = assignments[metabolite.id].get(nonH_formula, {})
                if not (metabolite.charge is None):
                    diff = metabolite.formula["H"] - metabolite.charge
                    equivalent_assignments = diff_seperated.get(diff, [])
                else:
                    equivalent_assignments = diff_seperated.get(None, [])
                for formula, charge in equivalent_assignments:
                    if formula == original_metabolite.formula and (charge == original_metabolite.charge):
                        metabolite.formula = formula
                        metabolite.charge = charge

    def adjust_protons(self):
        for reaction in self.model.reactions.values():
            if reaction in self.pseudo_reactions:
                continue
            self.proton_adjusted_reactions[reaction.id] = adjust_proton_count(reaction)
        
    @property
    def unbalancable_reactions(self):
        return self.balancer.unbalancable_reactions

    @property
    def unknown_metabolites(self):
        return self.balancer.unknown_metabolite_ids

    @property
    def reaction_reasons(self):
        return self.balancer.reaction_reasons

    @property
    def assignments(self):
        return self.balancer.assignments

    def add_unbalancable_reaction(self, reaction_id):
        self.balancer.unbalancable_reactions.add(self.model_interface.reactions[reaction_id])

    
    def generate_visual_report(self, filename=None):
        fig, image_ax = plt.subplots()
        fig.set_size_inches(6,6)
        size = 0.3

        cur_total = 0
        for reaction in self.model.reactions:
            if reaction.id in self.pseudo_reactions: continue
            if not is_balanced(reaction):
                cur_total += 1
        old_total = 0
        for reaction in self.original_model.reactions:
            if reaction.id in self.pseudo_reactions: continue
            if not is_balanced(reaction):
                old_total += 1

        metabolite_df = self.generate_metabolite_report()
        metabolite_df["Inferrence Type"] = metabolite_df.apply(lambda row: row["Inferrence Type"] if row["Inferrence Type"] != "Unconstrained" else f"Unconstrained {'with DB' if row['Used Databases'] != '' else 'without DB'}", axis=1)
        grouped_report = metabolite_df.groupby(["Inferrence Type", "Similarity"]).size().to_frame().reset_index().rename(columns={0 : "Count"})

        #values:
        outer_vals = np.array([grouped_report[grouped_report["Inferrence Type"] == i_type]["Count"].values for i_type in grouped_report["Inferrence Type"].unique()], dtype=object)
        inner_vals = grouped_report.groupby("Inferrence Type")["Count"].sum()

        # 4 colors or less:
        cmap = sns.color_palette("tab20b", as_cmap=True)
        inner_color_pos = np.array([i * 4 for i in range(len(inner_vals))]) + 4
        outer_color_pos = np.array([(i * 4) + np.arange(1, len(out_vals) + 1)  + 4 for i, out_vals in enumerate(outer_vals)], dtype=object)

        inner_colors = cmap(inner_color_pos)
        outer_colors = cmap(np.hstack(outer_color_pos).astype(int))

        #labels:
        inner_labels = ["" if val < sum(inner_vals)*.01 else val for val in inner_vals]
        outer_labels = ["" if val < sum(inner_vals)*.01 else val for val in np.hstack(outer_vals)]

        #drawing:
        outer_wedges, outer_texts = image_ax.pie(np.hstack(outer_vals), radius=1, colors=outer_colors, labels=outer_labels, labeldistance=1.1,
            wedgeprops=dict(width=size, edgecolor='w'))

        inner_wedges, inner_texts = image_ax.pie(inner_vals, radius=1-size, colors=inner_colors, labels=inner_labels, labeldistance=.75,
            wedgeprops=dict(width=size, edgecolor='w'))

        #add legend:
        legend_labels = [f"{i_type} / {same}" for i_type, same in zip(grouped_report["Inferrence Type"], grouped_report["Similarity"])]
        image_ax.legend(outer_wedges, [f"{val} - {label}" for val, label in zip(np.hstack(outer_vals), legend_labels)],
                title="Inferrence type / comparison with original model",
                loc="upper left",
                bbox_to_anchor=(1, 0, 1, 1))

        #set title
        image_ax.set(aspect="equal", title=f'Metabolite Formulae for {self.model.id}\nin comparison with original formulae\n/w {cur_total} from {old_total} originally unbalanced\n reactions within {self.total_time:.1f} seconds.')
        if not filename is None:
            plt.savefig(f'{filename}.png', dpi=400,bbox_inches='tight')
        return fig
#FIXME: adjust to libsbml
    def generate_reaction_report(self, filename=None, proton_threshold = 7):

        def get_shared_metabolites(reaction_ids):
            metabolites = []
            for reaction_id in reaction_ids:
                metabolites.append(set(self.model.reactions.get_by_id(reaction_id).metabolites.keys()))
            if len(metabolites) == 0: return set()
            else: return set.intersection(*metabolites)
        
        def difference_string(mass_dict):
            return ", ".join(f"{element}: {diff}" for element, diff in mass_dict.items() if diff != 0)

        reaction_report = []
        pseudo_reactions = get_pseudo_reactions(self.model)
        for reaction in self.model.reactions:
            if reaction.id in pseudo_reactions: continue
            balance = calculate_balance(reaction)
            unbalance_type = []
            if not all(np.isclose(val, 0) for val in balance["mass"].values()):
                unbalance_type.append("Mass") 
            if not (np.isclose(balance["charge"],0)):
                unbalance_type.append("Charge")
            if len(unbalance_type) != 0:
                reason = self.reaction_reasons.get(reaction.id, [])
                if len(reason) == 0:
                    reason == [reaction.id]
                reaction_report.append({"Id" : reaction.id,
                                        "Unbalanced Reaction" : str(reaction),
                                        "Unbalanced Type" : ", ".join(unbalance_type),
                                        "Reason" : ", ".join(reason),
                                        "Shared Metabolites" : ", ".join([m.id for m in get_shared_metabolites(reason)]),
                                        "Mass Difference": difference_string(balance["mass"]),
                                        "Charge Difference" : balance["charge"]
                })
            elif self.proton_adjusted_reactions[reaction.id] > proton_threshold:
                reaction_report.append({"Id" : reaction.id,
                                        "Unbalanced Reaction" : str(reaction),
                                        "Unbalanced Type" : "High Proton Count",
                                        "Reason" : f"Added {self.proton_adjusted_reactions[reaction.id]} protons.",
                                        "Shared Metabolites" : "",
                                        "Mass Difference": "",
                                        "Charge Difference" : 0
                })
        information_df = pd.DataFrame(reaction_report)
        def type_order_func(series):
            order = {"Mass, Charge": 0, "Mass": 3, "Charge": 6, "High Proton Count" : 9}[series["Unbalanced Type"]]
            return order
        #sort frame by most interesting information
        information_df["type_order"] = information_df.apply(type_order_func, axis = 1)
        information_df["num_reasons"] = information_df["Reason"].apply(lambda x: len(x.split(", ")))
        information_df = information_df.sort_values(["type_order", "num_reasons"]).reset_index(drop=True)
        information_df = information_df.drop(columns = ["type_order", "num_reasons"])
        if not filename is None: information_df.to_csv(f"{filename}.csv")
        return information_df
        

    def generate_metabolite_report(self, filename = None, target_model = None):
        other_model = self.original_model

        def other_metabolite(metabolite_id):
            return other_model.getSpecies(metabolite_id)

        def generate_metabolite_information(metabolite):
            plugin = get_fbc_plugin(metabolite)
            other = other_metabolite(metabolite.id)
            other_plugin = get_fbc_plugin(other)
            this_databases = set()
            for assignment, dbs in self.data_collector.get_assignments(metabolite, database_seperated = True).items():
                if same_formula(assignment[0], plugin.chemical_formula) and ((plugin.charge == assignment[1]) or (assignment[1] is None)):
                    this_databases.update([f"{db[0]}:{db[1]}" for db in dbs])
            other_databases = set()
            for assignment, dbs in self.data_collector.get_assignments(metabolite, database_seperated = True).items():
                if same_formula(assignment[0], other_plugin.chemical_formula) and ((other_plugin.charge == assignment[1]) or (assignment[1] is None)):
                    other_databases.update([f"{db[0]}:{db[1]}" for db in dbs])

            inferrence_type = "Clean"
            if (len(this_databases) == 0):
                inferrence_type = "Inferred"
            if "R" in plugin.chemical_formula:
                inferrence_type = "Unconstrained"
            
            if same_formula(clean_formula(plugin.chemical_formula), clean_formula(other_plugin.chemical_formula)) and (plugin.charge == other_plugin.charge):
                if len(this_databases) > 0:
                    target = "Target & "
                else:
                    target = "Target"
            else:
                target = ""

            result ={"Id" : metabolite.id,
                    "Name" : metabolite.name,
                    "Determined Formula" : plugin.chemical_formula,
                    "Determined Charge" : plugin.charge,
                    "Previous Formula" : other_plugin.chemical_formula,
                    "Previous Charge" : other_plugin.charge,
                    "Inferrence Type" : inferrence_type,
                    "Reasoning" : target + ', '.join(this_databases),
                    "Used Databases" : ', '.join(this_databases),
                    "Previous Databases" : ', '.join(other_databases)
            }
 
            if not (target_model is None):
                target = target_model.getSpecies(metabolite.id)
                target_plugin = get_fbc_plugin(target)
                target_databases = set()
                for assignment, dbs in self.data_collector.get_assignments(metabolite, database_seperated = True).items():
                    if same_formula(assignment[0], target_plugin.chemical_formula) and ((target_plugin.charge == assignment[1]) or (assignment[1] is None)):
                        target_databases.update([f"{db[0]}:{db[1]}" for db in dbs])
                result.update({
                    "Target Formula" : target_plugin.chemical_formula,
                    "Target Charge" : target_plugin.charge,
                    "Target Databases" : ', '.join(target_databases)
                })

            return result
        information = {}
        for metabolite in self.model.getListOfSpecies():
            information[metabolite.id] = generate_metabolite_information(metabolite)

        
        fixing_reactions = set()
        for reaction in self.model.getListOfReactions():
            if reaction.id in self.pseudo_reactions: continue
            unknown_count = 0
            for metabolite_id in get_sbml_metabolites(reaction):
                if information[metabolite_id]["Reasoning"] == "":
                    unknown_count += 1
            if unknown_count == 1:
                fixing_reactions.add(reaction.id)
        while len(fixing_reactions) > 0:
            fixing_reaction_id = fixing_reactions.pop()
            reasons = []
            fixed_metabolite_id = None
            for metabolite_id in get_sbml_metabolites(self.model.getReaction(fixing_reaction_id)):
                reason = information[metabolite_id]["Reasoning"]
                if reason == "": fixed_metabolite_id = metabolite_id
                else:
                    reasons.append(f"{metabolite_id} -> {reason}")
            if fixed_metabolite_id is not None: #otherwise we already fixed it
                reason = f"{fixing_reaction_id}: {'(' if len(reasons) > 1 else ''}{'; '.join(reasons)}{')' if len(reasons) > 1 else ''}"
                information[fixed_metabolite_id]["Reasoning"] = reason
                for reaction_id in self.balancer.metabolite_reactions[fixed_metabolite_id]:
                    if reaction_id in self.pseudo_reactions: continue
                    unknown_count = 0
                    for metabolite_id in get_sbml_metabolites(self.model.getReaction(reaction_id)):
                        if information[metabolite_id] == "": unknown_count += 1
                    if unknown_count == 1:
                        fixing_reactions.add(reaction.id)
        
        for metabolite_id, info in information.items():
            if not info["Reasoning"].startswith("Target"):
                if same_formula(info["Determined Formula"], info["Previous Formula"], ignore_rest = True) and (info["Determined Charge"] == info["Previous Charge"]):
                    info["Reasoning"] = f"unconstrained Target{' & ' if info['Reasoning'] != '' else ''}{info['Reasoning']}"


        information_df = pd.DataFrame(list(information.values()))
        def similarity(row, columns):
            base = columns[0]
            target = columns[1]
            same_f = same_formula(row[base[0]], row[target[0]])
            same_charge = row[base[1]] == row[target[1]]
            
            hydrogen_difference = formula_to_dict(row[base[0]]).get("H", 0) - formula_to_dict(row[target[0]]).get("H", 0)
            charge_difference = row[base[1]] - row[target[1]]
            same_f_nonH = same_formula(remove_H.sub(r"\1", row[base[0]]), remove_H.sub(r"\1", row[target[0]]))
            hc_same = hydrogen_difference == charge_difference
            if same_f and same_charge: return "Same"
            if same_f_nonH and hc_same: return "Proton Diff"
            return "Diff"

        information_df["Similarity"] = information_df.apply(similarity, columns = (("Determined Formula", "Determined Charge"), ("Previous Formula", "Previous Charge")), axis = 1)
        if not (target_model is None):
                information_df["Target Similarity"] = information_df.apply(similarity, columns = (("Determined Formula", "Determined Charge"), ("Target Formula", "Target Charge")), axis = 1)
                information_df["Target Change"] = information_df.apply(similarity, columns = (("Previous Formula", "Previous Charge"), ("Target Formula", "Target Charge")), axis = 1)
                information_df = information_df[["Id", "Name", "Determined Formula", "Determined Charge", "Target Formula", "Target Charge", "Reasoning", "Target Databases", "Target Similarity", "Inferrence Type", "Target Change", "Previous Formula", "Previous Charge", "Previous Databases", "Similarity"]]
        def type_order_func(series):
            order = {"Inferred": 0, "Unconstrained": 3, "Clean": 6}[series["Inferrence Type"]]
            order = order if series["Determined Formula"] == series["Previous Formula"] else order + 1.5
            order = order if series["Determined Charge"] == series["Previous Charge"] else order + .5
            return order
        #sort frame by most interesting information
        information_df["type_order"] = information_df.apply(type_order_func, axis = 1)
        information_df["num_db"] = information_df["Reasoning"].apply(lambda x: len(x.split(", ")))
        information_df = information_df.sort_values(["type_order", "num_db"]).reset_index(drop=True)
        information_df = information_df.drop(columns = ["type_order", "num_db"])
        if not filename is None: information_df.to_csv(f"{filename}.csv")
        return information_df
