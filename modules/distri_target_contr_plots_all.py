import matplotlib.pyplot as plt
import numpy as np
import os

def distri_target_contr_plots_all(binn, perc_z, perc_i, perc_t, perc_n, cond1, output_dir):
    """
    Plot LFC distribution for zGE, intergenic, target, and NT control.
    """
    # Plot 1: Separate subplots
    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.bar(binn, perc_z, width=np.diff(binn)[0], color='c')
    plt.grid(True)
    plt.title(f'LFC distribution zGE, one replica: {cond1}', fontsize=14)
    
    plt.subplot(3, 1, 2)
    plt.bar(binn, perc_i, width=np.diff(binn)[0], color='m')
    plt.grid(True)
    plt.title('LFC distribution intergenic', fontsize=12)
    
    plt.subplot(3, 1, 3)
    plt.bar(binn, perc_t, width=np.diff(binn)[0], color='k')
    plt.grid(True)
    plt.title('LFC distribution target genes', fontsize=12)
    plt.xlabel('LFC', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"distri_separate_target_controls_{cond1}_WTval.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Plot 2: Overlaid plot
    
    plt.figure(figsize=(10, 5))
    plt.bar(binn, perc_z, width=np.diff(binn)[0], color='c', label='zGE')
    plt.bar(binn, perc_i, width=np.diff(binn)[0], color='m', alpha=0.7, label='intergenic')
    plt.bar(binn, perc_t, width=np.diff(binn)[0], color='k', alpha=0.5, label='target')
    plt.bar(binn, perc_n, width=np.diff(binn)[0], color='g', alpha=0.3, label='NT control')

    plt.grid(True)
    plt.title(f'LFC distribution - Target Controls (one replica): {cond1}', fontsize=14)
    plt.xlabel('LFC', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,f"distri_target_controls_{cond1}_WTval.png"), dpi=300, bbox_inches="tight")
    plt.close()

    return 1