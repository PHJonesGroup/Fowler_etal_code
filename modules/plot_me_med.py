import numpy as np
import matplotlib.pyplot as plt
import os 

def plot_me_med(tab_norm_T0, tab_norm_T1, output_dir):

    T0 = tab_norm_T0.iloc[:, 2:4].to_numpy()
    T1 = tab_norm_T1.iloc[:, 2:4].to_numpy()

    me0 = np.mean(T0, axis=0)
    med0 = np.median(T0, axis=0)

    me1 = np.mean(T1, axis=0)
    med1 = np.median(T1, axis=0)

    me_med_T0_T1 = np.array([
        [me0[0], me0[1], me1[0], me1[1]],
        [med0[0], med0[1], med1[0], med1[1]]
    ])

    labels = ['T0F', 'T0M', 'T1F', 'T1M']
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8, 5))
    # Plot mean bars first
    ax.bar(x, me_med_T0_T1[0], label='Mean', color='green')
    # Overlay median bars on top
    ax.bar(x, me_med_T0_T1[1], label='Median', color='magenta')

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Normalized counts')
    ax.set_title('Mean and Median of normalized counts')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "normalised_counts_mean_med"), dpi=300, bbox_inches="tight")
    plt.close()

    return me_med_T0_T1