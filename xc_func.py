import numpy as np
import pylibxc


def get_Vxc_exc_lda(density_s: np.ndarray, x_func: pylibxc.LibXCFunctional, c_func: pylibxc.LibXCFunctional) -> (np.ndarray, np.ndarray):
    inp = {'rho': density_s}
    res_x = x_func.compute(inp)
    res_c = c_func.compute(inp)
    Vxc = (res_x['vrho'] + res_c['vrho'])
    exc = (res_x['zk'] + res_c['zk'])
    return Vxc, exc
