import numpy as np
from mgs.tools.pbc import apply_PBC

def distence_3d(x, y, pbc = False, Ls = None):

    if not isinstance(x, (list, np.ndarray)) or not isinstance(y, (list, np.ndarray)):
        return 'must be list or numpy array'
    if len(x) != 3 or len(y) != 3:
        return 'must be 3d'
    
    dx = x[0] - y[0]
    dy = x[1] - y[1]
    dz = x[2] - y[2]

    if pbc:
        if Ls == None:
            return 'Must specify the cell matrixs!'
        if len(Ls) != 3:
            return 'Must be 3d cell matrixs!'

        dx = apply_PBC(dx, Ls[0])
        dy = apply_PBC(dy, Ls[1])
        dz = apply_PBC(dz, Ls[2])
    
    return np.sqrt(dx**2 + dy**2 + dz**2)


def distence_2d(x, y, pbc = False, Ls = None):

    if not isinstance(x, (list, np.ndarray)) or not isinstance(y, (list, np.ndarray)):
        return 'must be list or numpy array'
    if len(x) != 2 or len(y) != 2:
        return 'must be 2d'

    dx = x[0] - y[0]
    dy = x[1] - y[1]

    if pbc:
        if Ls == None:
            return 'Must specify the cell matrixs!'
        if len(Ls) != 2:
            return 'Must be 2d cell matrixs!'
        
        dx = apply_PBC(dx, Ls[0])
        dy = apply_PBC(dy, Ls[1])

    return np.sqrt(dx**2 + dy**2)


