import scipy
import numpy as np
import interp
import radial_wave
import ode
import cbinding
import ctypes


def solve_poisson(density_: np.ndarray, grid: radial_wave.radial_grid, use_c: bool = False):
    if use_c:
        c_r0 = ctypes.c_double(grid.r0)
        c_rc = ctypes.c_double(grid.rc)
        c_N = ctypes.c_int(grid.N)
        VH = np.zeros(grid.N)
        dVH = np.zeros(grid.N)
        density_r = density_ * grid.r
        cbinding.libatomdft.solve_poisson(
            VH,
            dVH,
            np.ascontiguousarray(density_r),
            c_r0,
            c_rc,
            c_N,
        )
        return VH, dVH
    else:
        density_r = density_ * grid.r

        def func(i, y):
            return np.array([
                grid.r[i] * y[1], -2.0 * y[1] - 4.0 * np.pi * density_r[i]
            ])

        def func_mid(i, y):
            return np.array([
                grid.r0 * np.exp(grid.u[i]+0.5*grid.u_step) * y[1], -2.0 * y[1] - 4.0 * np.pi * density_r_mid[i]
            ])

        VH = np.zeros([grid.N, 2])
        # first three points
        VH[0, 0] = 4.0 * np.pi * scipy.integrate.simpson(density_r * grid.r) * grid.u_step
        VH[0, 1] = 0.0
        # get mid array
        density_r_mid = np.zeros(3)
        for i in range(3):
            density_r_mid[i] = interp.interp_poly(grid.u[i] + 0.5 * grid.u_step, density_r[0:4], grid.u[0:4])
        # runge kutta 4
        ode.runge_kutta_4(func, func_mid, grid.u_step, 4, VH)
        # adam 4
        ode.adam_4(func, grid.u_step, grid.N, VH)
        return VH[:, 0], VH[:, 1]
