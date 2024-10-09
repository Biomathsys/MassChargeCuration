import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

def visual_report(curator, filename=None, ax = None, size = None, dpi = 400, **kwargs):
    if ax is None:
        fig, image_ax = plt.subplots()
        fig.set_size_inches(6,6)
    else:
        image_ax = ax
        fig = ax.get_figure()

    if size is None:
        size = 0.3

    cur_total = 0
    for reaction in curator.model_interface.reactions.values():
        if reaction in curator.pseudo_reactions: continue
        if not reaction.is_balanced():
            cur_total += 1
    old_total = 0
    for reaction in curator.original_model_interface.reactions.values():
        if reaction in curator.pseudo_reactions: continue
        if not reaction.is_balanced():
            old_total += 1

    metabolite_df = curator.generate_metabolite_report()
    metabolite_df["Inferrence Type"] = metabolite_df.apply(lambda row: row["Inferrence Type"] if row["Inferrence Type"] != "Unconstrained" else f"Unconstrained {'with DB' if row['Used Databases'] != '' else 'without DB'}", axis=1)
    grouped_report = metabolite_df.groupby(["Inferrence Type", "Similarity"]).size().to_frame().reset_index().rename(columns={0 : "Count"})

    #values:
    outer_vals = np.array([grouped_report[grouped_report["Inferrence Type"] == i_type]["Count"].values for i_type in grouped_report["Inferrence Type"].unique()], dtype=object)
    inner_vals = grouped_report.groupby("Inferrence Type")["Count"].sum()

    # 4 colors or less:
    cmap = cm.get_cmap("tab20b")
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
    image_ax.set(aspect="equal", title=f'Metabolite Formulae for {curator.model_interface.get_model_id()}\nin comparison with original formulae\n/w {cur_total} from {old_total} originally unbalanced\n reactions within {curator.total_time:.1f} seconds.')
    if not filename is None:
        plt.savefig(f'{filename}.png', dpi=dpi, bbox_inches='tight', **kwargs)
    return fig