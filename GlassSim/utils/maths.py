import numpy as np


def Sin(x, A0, A1, s):
    return A0+A1*np.sin(2*np.pi*x+s)