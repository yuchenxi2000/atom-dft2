import configuration
import numpy as np
import scipy
import radial_wave
import dataclasses


@dataclasses.dataclass
class total_energy:
    Etot: float
    Ekin: float
    EH: float
    Eenuc: float
    Exc: float
    eig_sum: float


def get_total_energy(levels: list[configuration.Level], exc: np.ndarray, Vxc: np.ndarray, VH: np.ndarray, density_s: np.ndarray, Z: float, grid: radial_wave.radial_grid) -> total_energy:
    r = grid.r
    EH = 1/2 * 4 * np.pi * scipy.integrate.simpson(VH * np.sum(density_s, axis=1) * r ** 3) * grid.u_step
    Vxc_n = 4 * np.pi * scipy.integrate.simpson(np.sum(Vxc * density_s, axis=1) * r ** 3) * grid.u_step
    Exc = 4 * np.pi * scipy.integrate.simpson(np.sum(exc * density_s, axis=1) * r ** 3) * grid.u_step
    eig_sum = 0.0
    for orb in levels:
        eig_sum += orb.occ * orb.eig
    Eenuc = -Z * 4 * np.pi * scipy.integrate.simpson(np.sum(density_s, axis=1) * r ** 2) * grid.u_step
    Ekin = eig_sum - 2 * EH - Eenuc - Vxc_n
    Etot = eig_sum - EH + Exc - Vxc_n
    return total_energy(Etot, Ekin, EH, Eenuc, Exc, eig_sum)
