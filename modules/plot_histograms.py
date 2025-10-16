import matplotlib.pyplot as plt
import os

def plot_histograms(binn, perc_t1, perc_zt1, perc_t2, perc_zt2, condition_label, output_dir):

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 2, 1)
    plt.bar(binn, perc_t1, color='k')
    plt.grid(True)
    plt.title('LFC distri rep1')
    plt.xlim([-16, 8])
    plt.xlabel('LFC')

    plt.subplot(2, 2, 3)
    plt.bar(binn, perc_zt1, color='g')
    plt.grid(True)
    plt.title('Z LFC distri rep1')
    plt.xlim([-16, 8])
    plt.xlabel('LFC Z-normalised')

    plt.subplot(2, 2, 2)
    plt.bar(binn, perc_t2, color='k')
    plt.grid(True)
    plt.title('LFC distri rep2')
    plt.xlim([-16, 8])
    plt.xlabel('LFC')

    plt.subplot(2, 2, 4)
    plt.bar(binn, perc_zt2, color='g')
    plt.grid(True)
    plt.title('Z LFC distri rep2')
    plt.xlim([-16, 8])
    plt.xlabel('LFC Z-normalised')

    plt.suptitle(f'LFC gRNA distri: {condition_label}', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f"target_distri_LFC_Z_corr_{condition_label}"), dpi=300, bbox_inches="tight")
    plt.close()