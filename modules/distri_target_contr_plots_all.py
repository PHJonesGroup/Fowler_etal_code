import matplotlib.pyplot as plt
import numpy as np
import os

def distri_target_contr_plots_all(binn, perc_z, perc_i, perc_t, perc_n, cond1, output_dir):
    """
    Plot LFC distribution for zGE, intergenic, target, and NT control.
    Any category whose array is None is skipped (e.g. zGE when absent).
    """
    width = np.diff(binn)[0]

    # ---- Figure 1: separate panels (only those with data) ----
    panels = [
        (perc_z, 'Zero Expressed Genes',   'c'),
        (perc_i, 'Intergenic Genes',       'm'),
        (perc_n, 'Non-Targetting Genes',   'orange'),
        (perc_t, 'Target Genes',           'k'),
    ]
    panels = [p for p in panels if p[0] is not None]

    plt.figure(figsize=(10, 2 * len(panels)))
    plt.suptitle(f"Condition: {cond1}", fontsize=16)
    for k, (perc, title, color) in enumerate(panels, start=1):
        plt.subplot(len(panels), 1, k)
        plt.bar(binn, perc, width=width, color=color)
        plt.grid(True)
        plt.title(title, fontsize=12)
        plt.ylabel('%', fontsize=12)
        if k == len(panels):                 # xlabel only on bottom panel
            plt.xlabel('LFC', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"distri_separate_target_controls_{cond1}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # ---- Figure 2: overlaid (only categories with data) ----
    plt.figure(figsize=(10, 5))
    if perc_z is not None:
        plt.bar(binn, perc_z, width=width, color='c',      label='zGE')
    if perc_i is not None:
        plt.bar(binn, perc_i, width=width, color='m', alpha=0.7, label='intergenic')
    if perc_t is not None:
        plt.bar(binn, perc_t, width=width, color='k', alpha=0.5, label='target')
    if perc_n is not None:
        plt.bar(binn, perc_n, width=width, color='g', alpha=0.3, label='NT control')
    plt.grid(True)
    plt.title(f'Condition : {cond1}', fontsize=14)
    plt.xlabel('LFC', fontsize=12)
    plt.ylabel('%', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"distri_target_controls_{cond1}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return 1