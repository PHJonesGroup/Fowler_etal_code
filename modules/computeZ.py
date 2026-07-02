import numpy as np
import matplotlib.pyplot as plt
from .make_histo_LFC import make_histo_LFC

def computeZ(st, en, step, LFC_t, me_sd_z):
    """
    Compute Z scores for targets based on mu and sd from control (CTR).

    Inputs:
    - LFC_t1: numpy arrays of LFC values
    - me_sd_z12: 2x2 array with mean and sd (rows)
    - st, en, step: float, histogram bin parameters

    Outputs:
    - Z_t1: Z-normalized LFC vectors 
    - binn: bin centers for histograms
    - perc_t, perc_zt: percentage histograms for LFC and Z-LFC
    """
    mu1, sd1 = me_sd_z[0], me_sd_z[1]

    binn, _, perc_t1 = make_histo_LFC(step, LFC_t, st, en)
    Z_t1 = (LFC_t - mu1) / sd1
    _, _, perc_zt1 = make_histo_LFC(step, Z_t1, st, en)

    return Z_t1, binn, perc_t1, perc_zt1
