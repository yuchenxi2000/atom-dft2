import numpy as np


def interp_poly(z: float, y: np.ndarray, x: np.ndarray):
    N = x.shape[0]
    s = 0.0
    for i in range(N):
        prod = 1.0
        for j in range(N):
            if j != i:
                prod *= (z - x[j]) / (x[i] - x[j])
        s += y[i] * prod
    return s
