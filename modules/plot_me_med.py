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

    return me_med_T0_T1