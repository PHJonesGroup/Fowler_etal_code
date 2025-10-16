import numpy as np
import matplotlib.pyplot as plt
from .make_histo_LFC import make_histo_LFC

def computeZ(st, en, step, LFC_t1, LFC_t2, me_sd_z12):
    """
    Compute Z scores for targets based on mu and sd from control (CTR).

    Inputs:
    - LFC_t1, LFC_t2: numpy arrays of LFC values for replicate 1 and 2
    - me_sd_z12: 2x2 array with mean and sd for replicate 1 and 2 (rows)
    - st, en, step: float, histogram bin parameters

    Outputs:
    - Z_t1, Z_t2: Z-normalized LFC vectors for replicate 1 and 2
    - binn: bin centers for histograms
    - perc_t1, perc_zt1, perc_t2, perc_zt2: percentage histograms for LFC and Z-LFC
    """
    mu1, sd1 = me_sd_z12[0, 0], me_sd_z12[0, 1]
    mu2, sd2 = me_sd_z12[1, 0], me_sd_z12[1, 1]

    binn, _, perc_t1 = make_histo_LFC(step, LFC_t1, st, en)
    Z_t1 = (LFC_t1 - mu1) / sd1
    _, _, perc_zt1 = make_histo_LFC(step, Z_t1, st, en)

    _, _, perc_t2 = make_histo_LFC(step, LFC_t2, st, en)
    Z_t2 = (LFC_t2 - mu2) / sd2
    _, _, perc_zt2 = make_histo_LFC(step, Z_t2, st, en)

    return Z_t1, Z_t2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2
