import numpy as np
import pylibxc


def get_Vxc_exc(density_s: np.ndarray, x_func: pylibxc.LibXCFunctional, c_func: pylibxc.LibXCFunctional) -> (np.ndarray, np.ndarray):
    inp = {'rho': density_s}
    res_x = x_func.compute(inp)
    res_c = c_func.compute(inp)
    Vxc = (res_x['vrho'] + res_c['vrho'])
    exc = (res_x['zk'] + res_c['zk'])
    return Vxc, exc


# TODO: begin test


def get_Y(y, b, c):
    return y**2 + b*y + c


def getvxc_scalar(n):
    y0 = -0.10498
    b = 3.72744
    c = 12.9352
    A = 0.0621814

    if n == 0:
        return 0.0, 0.0

    Q = np.sqrt(4 * c - b ** 2)
    rs = (3 / (4 * np.pi * n)) ** (1.0 / 3)
    y = np.sqrt(rs)
    ec = A / 2 * (np.log(y ** 2 / get_Y(y, b, c)) + 2 * b / Q * np.arctan(Q / (2 * y + b))- b * y0 / get_Y(y0, b, c) * (np.log((y-y0) ** 2 / get_Y(y, b, c)) + 2 * (b + 2 * y0) / Q * np.arctan(Q / (2 * y + b))))
    Vc = ec - A / 6 * (c * (y - y0) - b * y0 * y) / ((y - y0) * get_Y(y, b, c))
    ex = -3 / (4 * np.pi) * (3 * np.pi ** 2 * n) ** (1.0 / 3)
    Vx = 4 * ex / 3

    exc = ex + ec
    Vxc = Vx + Vc
    return Vxc, exc


def get_Vxc_exc_test(density_s: np.ndarray, x_func: pylibxc.LibXCFunctional, c_func: pylibxc.LibXCFunctional) -> (np.ndarray, np.ndarray):
    N = density_s.shape[0]
    Vxc = np.zeros([N, 1])
    exc = np.zeros(N)
    for i in range(N):
        V, e = getvxc_scalar(density_s[i, 0])
        Vxc[i, 0] = V
        exc[i] = e
    return Vxc, exc


# TODO: end test
