import numpy as np

def make_histo_LFC(step, vec_num, st, en):
    """
    Compute histogram counts and percentages with fixed bin size.

    Parameters:
    - step: bin width
    - vec_num: 1D array-like numeric values (e.g., LFC values)
    - st: minimum value (start of bins)
    - en: maximum value (end of bins)

    Returns:
    - bin: array of bin edges excluding the last edge 
    - his: counts of values in each bin
    - perc: percentage of counts in each bin
    """

    # Create bins from st to en with step size
    bins = np.arange(st, en + step, step)

    # Initialize counts array (length = number of bins - 1)
    val = np.zeros(len(bins) - 1, dtype=int)

    # Count values in each bin (vectorized approach)
    for i in range(len(bins) - 1):
        val[i] = np.sum((vec_num >= bins[i]) & (vec_num < bins[i + 1]))

    # Compute percentage
    total = val.sum()
    perc = 100 * val / total if total > 0 else np.zeros_like(val)

    # Return bins excluding the last edge
    bin_return = bins[:-1]

    return bin_return, val, perc

