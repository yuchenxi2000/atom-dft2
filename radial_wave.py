import matplotlib.pyplot as plt
import numpy as np
import dataclasses
import scipy
import cbinding
import ctypes
import ode


@dataclasses.dataclass
class radial_grid:
    N: int
    u: np.ndarray
    r: np.ndarray
    u_step: float
    r0: float
    rc: float

    @classmethod
    def exp_grid(cls, N, r0, rc):
        u_step = np.log(rc/r0) / (N-1)
        u = np.linspace(0.0, np.log(rc/r0), num=N)
        r = r0 * np.exp(u)
        return cls(N, u, r, u_step, r0, rc)


def solve_radial_from_zero(E: float, l: int, Z: float, Vext: np.ndarray, grid: radial_grid, use_c: bool = False):
    if use_c:
        chi = np.zeros(grid.N, dtype=np.float64, order="C")
        dchi = np.zeros(grid.N, dtype=np.float64, order="C")
        num_nodes_ptr = ctypes.pointer(ctypes.c_int(0))
        status = cbinding.libatomdft.solve_radial_from_zero(
            chi,
            dchi,
            num_nodes_ptr,
            ctypes.c_double(E),
            ctypes.c_int(l),
            ctypes.c_double(Z),
            np.ascontiguousarray(Vext),
            ctypes.c_double(grid.r0),
            ctypes.c_double(grid.rc),
            ctypes.c_int(grid.N),
        )
        if status != 0:
            raise Exception(f'call to C function solve_radial_from_zero failed with status {status}!')
        wfn = np.zeros([grid.N, 2])
        wfn[:, 0] = chi
        wfn[:, 1] = dchi
        return wfn, num_nodes_ptr.contents.value
    else:
        r = grid.r
        p2r = 2.0 * (E * r - 0.5 * l * (l + 1) / r + Z - Vext * r)

        def coeff_func(i):
            mat = np.array([
                [0.0, r[i]],
                [-p2r[i], 0.0],
            ])
            vec = np.zeros(2)
            return mat, vec

        wfn = np.zeros([grid.N, 2])
        for i in range(3):
            wfn[i, 0] = r[i] ** (l+1)
            wfn[i, 1] = (l + 1) * r[i] ** l
        ode.implicit_adam_4(coeff_func, grid.u_step, grid.N, wfn)

        num_nodes = count_nodes(wfn, grid.N)

        return wfn, num_nodes


def solve_radial_from_zero_and_inf(E: float, l: int, Z: float, Vext: np.ndarray, grid: radial_grid, use_c: bool = False):
    if use_c:
        chi = np.zeros(grid.N, dtype=np.float64, order="C")
        dchi = np.zeros(grid.N, dtype=np.float64, order="C")
        dE = np.zeros(1, dtype=np.float64, order="C")
        num_nodes_ptr = ctypes.pointer(ctypes.c_int(0))
        status = cbinding.libatomdft.solve_radial_from_zero_and_inf(
            chi,
            dchi,
            dE,
            num_nodes_ptr,
            ctypes.c_double(E),
            ctypes.c_int(l),
            ctypes.c_double(Z),
            np.ascontiguousarray(Vext),
            ctypes.c_double(grid.r0),
            ctypes.c_double(grid.rc),
            ctypes.c_int(grid.N)
        )
        wfn = np.zeros([grid.N, 2])
        wfn[:, 0] = chi
        wfn[:, 1] = dchi
        return status == 0, wfn, dE[0], num_nodes_ptr.contents
    else:
        if E >= 0.0:
            return False, None, 0.0, 0

        r = grid.r
        p2r = 2.0 * (E * r - 0.5 * l * (l + 1) / r + Z - Vext * r)

        # find last zero
        idx_last_zero = grid.N-1
        for idx_last_zero in range(grid.N-1, 0, -1):
            if p2r[idx_last_zero] == 0.0:
                break
            elif p2r[idx_last_zero] < 0.0 < p2r[idx_last_zero-1] or p2r[idx_last_zero] > 0.0 > p2r[idx_last_zero-1]:
                break

        if idx_last_zero == grid.N-1 or idx_last_zero == 1:
            return False, None, 0.0, 0

        # find inf
        pseudo_phase = 36.0
        idx_inf = idx_last_zero
        for idx_inf in range(idx_last_zero, grid.N):
            if np.sqrt(-2.0 * E) * (r[idx_inf] - r[idx_last_zero]) >= pseudo_phase:
                break

        if idx_inf == grid.N-1:
            return False, None, 0.0, 0

        if idx_inf - idx_last_zero < 10:
            return False, None, 0.0, 0

        if r[idx_last_zero] / r[idx_inf] > 0.5:
            return False, None, 0.0, 0

        def coeff_func(i):
            mat = np.array([
                [0.0, r[i]],
                [-p2r[i], 0.0],
            ])
            vec = np.zeros(2)
            return mat, vec

        def coeff_func_reverse(i):
            mat = np.array([
                [0.0, r[i + idx_last_zero]],
                [-p2r[i + idx_last_zero], 0.0],
            ])
            vec = np.zeros(2)
            return mat, vec

        wfn = np.zeros([grid.N, 2])
        # solve from zero
        for i in range(3):
            wfn[i, 0] = r[i] ** (l+1)
            wfn[i, 1] = (l+1) * r[i] ** l
        ode.implicit_adam_4(coeff_func, grid.u_step, idx_last_zero+1, wfn)
        # backup
        wfn_last_zero = np.copy(wfn[idx_last_zero])
        # solve from inf
        for i in range(idx_inf-3, idx_inf):
            wfn[i, 0] = np.exp(-np.sqrt(-2.0 * E) * (r[i] - r[idx_last_zero]))
            wfn[i, 1] = -np.sqrt(-2.0 * E) * wfn[i, 0]
        ode.implicit_adam_4_reversed(coeff_func_reverse, grid.u_step, idx_inf - idx_last_zero + 1, wfn[idx_last_zero:, :])
        factor = wfn_last_zero[0] / wfn[idx_last_zero, 0]
        wfn[idx_last_zero:idx_inf+1, :] *= factor
        # count number of nodes
        num_nodes = count_nodes(wfn, idx_inf+1)
        # TODO: should we fill tail?
        S = scipy.integrate.simpson(wfn[:, 0] ** 2 * r) * grid.u_step
        dE = wfn_last_zero[0] * (wfn_last_zero[1] - wfn[idx_last_zero, 1]) / (2 * S)
        return True, wfn, dE, num_nodes


def count_nodes(wfn: np.ndarray, N: int) -> int:
    num_nodes = 0
    for i in range(N-1):
        if abs(wfn[i, 0]) > 1e10:
            break
        if wfn[i, 0] == 0.0:
            num_nodes += 1
        elif wfn[i, 0] > 0.0 > wfn[i + 1, 0] or wfn[i, 0] < 0.0 < wfn[i + 1, 0]:
            num_nodes += 1
    return num_nodes


def get_eigen_radial_bisect(Emin: float, Emax: float, l: int, nr: int, Z: float, Vext: np.ndarray, grid: radial_grid, use_c: bool = False) -> float:
    while Emax - Emin >= 1e-9:
        E = 0.5 * (Emin + Emax)
        wfn, num_nodes = solve_radial_from_zero(E, l, Z, Vext, grid, use_c=use_c)
        if num_nodes <= nr:
            Emin = E
        else:
            Emax = E
    return Emin


def get_eigen_radial_combined(Emin: float, Emax: float, l: int, nr: int, Z: float, Vext: np.ndarray, grid: radial_grid, use_c: bool = False):
    converged = False
    E = 0.5 * (Emin + Emax)
    use_perturb = False
    while not converged:
        # print(f'E = {E}')
        if Emax - Emin < 1e-2:
            # try perturbation method
            success, wfn, dE, num_nodes = solve_radial_from_zero_and_inf(E, l, Z, Vext, grid, use_c=use_c)
            if success and num_nodes == nr:
                # print(f'dE = {dE}')
                if dE < 0.0:
                    Emax = E
                else:
                    Emin = E
                if E + dE > Emax or E + dE < Emin:
                    pass
                else:
                    E = E + dE
                    converged = abs(dE) < 1e-9
                    use_perturb = True
                    continue
        E = 0.5 * (Emin + Emax)
        wfn, num_nodes = solve_radial_from_zero(E, l, Z, Vext, grid, use_c=use_c)
        use_perturb = False
        if num_nodes <= nr:
            Emin = E
        else:
            Emax = E
        converged = Emax - Emin < 1e-9
    if not use_perturb:
        fix_wave_tail(wfn, grid.N)
    return E, wfn


def fix_wave_tail(wfn: np.ndarray, N: int) -> None:
    i = N-1
    for i in range(N-1, 0, -1):
        if abs(wfn[i-1, 0]) > abs(wfn[i, 0]):
            break
    if i <= 1:
        plt.plot(wfn)
        plt.ylim([-2, 2])
        plt.show()
        raise Exception('this wave function does not have peaks!')
    for j in range(i, N):
        wfn[j, :] = 0.0


def normalize_wave(wfn: np.ndarray, grid: radial_grid):
    norm = 4.0 * np.pi * scipy.integrate.simpson(wfn[:, 0] ** 2 * grid.r) * grid.u_step
    wfn /= np.sqrt(norm)
