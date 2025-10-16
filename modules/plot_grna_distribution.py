import matplotlib.pyplot as plt
import os

def plot_grna_distribution(gg, output_dir):
    plt.figure()
    plt.hist(gg, bins=300, edgecolor='black')
    plt.grid(True)
    plt.xlim(0, 6)
    plt.xlabel('number of gRNA per gene')
    plt.title('number of gRNA per gene distribution')
    plt.savefig(os.path.join(output_dir, "gene_distribution_gRNA"), dpi=300, bbox_inches="tight")
    plt.close()
    