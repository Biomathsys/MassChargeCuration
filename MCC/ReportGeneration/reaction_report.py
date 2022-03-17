import numpy as np
import pandas as pd

def reaction_report(curator, filename = None, proton_threshold = 7):

    def get_shared_metabolites(reaction_ids):
        metabolites = []
        for reaction_id in reaction_ids:
            metabolites.append(set(curator.model_interface.reactions[reaction_id].metabolites))
        if len(metabolites) == 0: return set()
        else: return set.intersection(*metabolites)
    
    def difference_string(mass_dict):
        return ", ".join(f"{element}: {diff}" for element, diff in mass_dict.items() if diff != 0)

    reaction_report = []
    for reaction in curator.model_interface.reactions.values():
        if reaction in curator.pseudo_reactions: continue
        mass_balance = reaction.mass_balance()
        charge_balance = reaction.charge_balance()
        unbalance_type = []
        if not all(np.isclose(val, 0) for val in mass_balance.values()):
            unbalance_type.append("Mass") 
        if not (np.isclose(charge_balance, 0)):
            unbalance_type.append("Charge")
        if len(unbalance_type) != 0:
            reason = curator.reaction_reasons.get(reaction.id, [])
            if len(reason) == 0:
                reason == [reaction.id]
            reaction_report.append({"Id" : reaction.id,
                                    "Unbalanced Reaction" : str(reaction),
                                    "Unbalanced Type" : ", ".join(unbalance_type),
                                    "Reason" : ", ".join(reason),
                                    "Shared Metabolites" : ", ".join([m.id for m in get_shared_metabolites(reason)]),
                                    "Mass Difference": difference_string(mass_balance),
                                    "Charge Difference" : charge_balance
            })
        elif curator.proton_adjusted_reactions[reaction.id] > proton_threshold:
            reaction_report.append({"Id" : reaction.id,
                                    "Unbalanced Reaction" : str(reaction),
                                    "Unbalanced Type" : "High Proton Count",
                                    "Reason" : f"Added {curator.proton_adjusted_reactions[reaction.id]} protons.",
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