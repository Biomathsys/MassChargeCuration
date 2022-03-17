from ..ModelInterface.ModelInterface import *
import pandas as pd

def metabolite_report(curator, filename = None, original_model = None, target_model = None):

        if not (target_model is None):
            target_model_interface = ModelInterface(target_model)

        if not (original_model is None):
            original_model_interface = ModelInterface(original_model)
        else:
            original_model_interface = curator.original_model_interface

        def generate_metabolite_information(metabolite):
            original = original_model_interface.metabolites[metabolite.id]
            this_databases = set()
            for (formula, charge), dbs in curator.data_collector.get_assignments(metabolite, database_seperated = True).items():
                if (metabolite.formula == formula) and ((metabolite.charge == charge) or (charge is None)):
                    this_databases.update([f"{db[0]}:{db[1]}" for db in dbs])
            original_databases = set()
            for (formula, charge), dbs in curator.data_collector.get_assignments(metabolite, database_seperated = True).items():
                if (original.formula == formula) and ((original.charge == charge) or (charge is None)):
                    original_databases.update([f"{db[0]}:{db[1]}" for db in dbs])

            inferrence_type = "Clean"
            if (len(this_databases) == 0):
                inferrence_type = "Inferred"
            if "R" in metabolite.formula:
                inferrence_type = "Unconstrained"
            
            if (metabolite.formula == original.formula) and (metabolite.charge == original.charge):
                if len(this_databases) > 0:
                    target = "Adherence & "
                else:
                    target = "Adherence"
            else:
                target = ""

            result ={"Id" : metabolite.id,
                    "Name" : metabolite.name,
                    "Determined Formula" : str(metabolite.formula),
                    "Determined Charge" : metabolite.charge,
                    "Previous Formula" : str(original.formula),
                    "Previous Charge" : original.charge,
                    "Inferrence Type" : inferrence_type,
                    "Reasoning" : target + ', '.join(this_databases),
                    "Used Databases" : ', '.join(this_databases),
                    "Previous Databases" : ', '.join(original_databases)
            }
 
            if not (target_model is None):
                target = target_model_interface.metabolites[metabolite.id]
                target_databases = set()
                for (formula, charge), dbs in curator.data_collector.get_assignments(metabolite, database_seperated = True).items():
                    if (formula == target.formula) and ((target.charge == charge) or (charge is None)):
                        target_databases.update([f"{db[0]}:{db[1]}" for db in dbs])
                result.update({
                    "Target Formula" : str(target.formula),
                    "Target Charge" : target.charge,
                    "Target Databases" : ', '.join(target_databases)
                })

            return result
        information = {}
        for metabolite in curator.model_interface.metabolites.values():
            information[metabolite.id] = generate_metabolite_information(metabolite)

        
        # more complicated reasoning
        fixing_reactions = set()
        for reaction in curator.model_interface.reactions.values():
            if reaction in curator.pseudo_reactions: continue
            unknown_count = 0
            for metabolite in reaction.metabolites:
                if information[metabolite.id]["Reasoning"] == "":
                    unknown_count += 1
            if unknown_count == 1:
                fixing_reactions.add(reaction)
        while len(fixing_reactions) > 0:
            fixing_reaction = fixing_reactions.pop()
            reasons = []
            fixed_metabolite = None
            for metabolite in fixing_reaction.metabolites:
                reason = information[metabolite.id]["Reasoning"]
                if reason == "": fixed_metabolite = metabolite
                else:
                    reasons.append(f"{metabolite.id} -> {reason}")
            if fixed_metabolite is not None: #originalwise we already fixed it
                reason = f"{fixing_reaction.id}: {'(' if len(reasons) > 1 else ''}{'; '.join(reasons)}{')' if len(reasons) > 1 else ''}"
                information[fixed_metabolite.id]["Reasoning"] = reason
                for reaction in fixed_metabolite.reactions:
                    if reaction in curator.pseudo_reactions: continue
                    unknown_count = 0
                    for metabolite in reaction.metabolites:
                        if information[metabolite.id] == "": unknown_count += 1
                    if unknown_count == 1:
                        fixing_reactions.add(reaction)
        
        for metabolite_id, info in information.items():
            if not info["Reasoning"].startswith("Target"):
                determined_formula_no_R = Formula(info["Determined Formula"])
                determined_formula_no_R["R"] = 0
                previous_formula_no_R = Formula(info["Previous Formula"])
                previous_formula_no_R["R"] = 0
                if (determined_formula_no_R == previous_formula_no_R) and (info["Determined Charge"] == info["Previous Charge"]):
                    info["Reasoning"] = f"unconstrained Target{' & ' if info['Reasoning'] != '' else ''}{info['Reasoning']}"

        # writing to the dataframe and determine order
        information_df = pd.DataFrame(list(information.values()))
        def similarity(row, columns):
            base = columns[0]
            target = columns[1]
            f1 = Formula(row[base[0]])
            f2 = Formula(row[target[0]])
            same_f = f1 == f2 
            same_charge = row[base[1]] == row[target[1]]
            
            hydrogen_difference = f1["H"] - f2["H"]
            charge_difference = row[base[1]] - row[target[1]]
            f1["H"] = 0
            f2["H"] = 0
            same_f_nonH = f1 == f2
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
