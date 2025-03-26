import numpy as np

def apply_PBC(M, L):
    if not isinstance(M, np.ndarray):
        M = np.array(M)
    if not isinstance(L, np.ndarray):
        L = np.array(L)
    M = M - L*np.round(M/L)
    return M