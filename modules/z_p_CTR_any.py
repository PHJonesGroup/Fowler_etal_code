import pandas as pd
import numpy as np
from .p_control_any_implement import p_control_any_implement
from .make_tables_Z_two import make_tables_Z_two
from .computeZ import computeZ

def z_p_CTR_any(
    st: float,
    en: float,
    step: float,
    T_any: pd.DataFrame,
    binn: np.ndarray,
    p_cont12: np.ndarray,
    me_sd_12: np.ndarray,
    cond1: str,
    cond2: str,
):
    """
    Calibrate any control set (e.g. intergenic or NT) against zGE controls,
    assign p/q values, compute Z‑scores, and return per‑replicate tables + histos.
    """
    T_ij_q_verti_any = p_control_any_implement(T_any, binn, p_cont12)

    LFC_t1 = T_ij_q_verti_any.iloc[:, 2].to_numpy()
    LFC_t2 = T_ij_q_verti_any.iloc[:, 3].to_numpy()

    Z_t1, Z_t2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2 = computeZ(
        st, en, step, LFC_t1, LFC_t2, me_sd_12
    )

    T_z_q_any1, T_z_q_any2 = make_tables_Z_two(
        T_ij_q_verti_any, Z_t1, Z_t2, cond1, cond2
    )

    return T_z_q_any1, T_z_q_any2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2
