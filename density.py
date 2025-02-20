import numpy as np
import scipy
import configuration
import radial_wave


def get_density_r2s(levels: list[configuration.Level], N: int, nspin: int):
    density_r2s = np.zeros([N, nspin])
    for level in levels:
        density_r2s[:, level.s] += level.occ * level.wfn[:, 0] ** 2
    return density_r2s


def density_mixing_linear(density_old: np.ndarray, density_new: np.ndarray, mixing: float):
    return (1 - mixing) * density_old + mixing * density_new


def least_square(left_vecs, right_vec: np.ndarray, inner_prod) -> np.ndarray:
    if isinstance(left_vecs, list):
        N = len(left_vecs)
    elif isinstance(left_vecs, np.ndarray):
        N = left_vecs.shape[0]
    else:
        raise Exception(f'unknown type of left_vecs! ({type(left_vecs)})')
    left_mat = np.zeros([N, N])
    right_vec_ = np.zeros(N)
    for i in range(0, N):
        for j in range(i+1, N):
            left_mat[i, j] = inner_prod(left_vecs[i], left_vecs[j])
    left_mat = left_mat + left_mat.T
    for i in range(0, N):
        left_mat[i, i] = inner_prod(left_vecs[i], left_vecs[i])
        right_vec_[i] = inner_prod(left_vecs[i], right_vec)
    return np.linalg.solve(left_mat, right_vec_)


def density_mixing_anderson(history_density: list[np.ndarray], history_residue: list[np.ndarray], beta: float) -> np.ndarray:
    # ordering in history: oldest first, newest last
    # len(history_density) = N
    # len(history_residue) = N
    # len(beta) = N
    N = len(history_density)
    if N <= 0:
        raise Exception(f'anderson mixing requires at least one points, but I found {N} points!')
    if len(history_residue) != N:
        raise Exception(f'number of residue vectors ({len(history_residue)}) must equal to number of density vectors ({N})!')

    def inner_prod(a, b):
        return np.tensordot(a, b)

    # get alpha
    alpha = []
    if N == 1:
        alpha.append(1.0)
    else:
        df = []
        for i in range(N - 1):
            df.append(history_residue[i + 1] - history_residue[i])
        gamma = least_square(df, history_residue[N - 1], inner_prod)
        for i in range(N-1):
            if i == 0:
                alpha.append(gamma[i])
            else:
                alpha.append(gamma[i] - gamma[i-1])
        alpha.append(1.0 - gamma[-1])
    x = np.zeros_like(history_density[0])
    for i in range(N):
        x += alpha[i] * (history_density[i] + beta * history_residue[i])
    return x


def get_approx_density_gaussian(grid: radial_wave.radial_grid, Z: float, ra: float, nspin: int):
    density_r2s = np.zeros([grid.N, nspin])
    r = grid.r
    for ispin in range(nspin):
        density_r2s[:, ispin] = r ** 2 / np.power(2 * np.pi, 3 / 2) / ra ** 3 * np.exp(-(r / ra) ** 2 / 2) / nspin
    density_r2s *= Z
    # nelect = 4.0 * np.pi * scipy.integrate.simpson(np.sum(density_r2s, axis=1) * r) * grid.u_step
    # print(nelect)
    return density_r2s


def get_approx_density_thomas_fermi(grid: radial_wave.radial_grid, Z: float, nspin: int) -> np.ndarray:
    r = grid.r
    x = r * np.power(128.0 * Z / (9 * np.pi**2), 1.0/3.0)
    a = 0.7280642371
    b = -0.5430794693
    c = 0.3612163121
    Zeff = Z * (1.0 + a * np.sqrt(x) + b * x * np.exp(-c * np.sqrt(x))) ** 2 * np.exp(-2.0 * a * np.sqrt(x))
    V = -Zeff / r
    rho = -1.0 / (3.0 * np.pi ** 2) * np.power(-2.0 * V, 3.0/2.0)
    # normalize
    nelect = 4.0 * np.pi * scipy.integrate.simpson(rho * r ** 3) * grid.u_step
    rho = rho / nelect * Z
    density_s = np.zeros([grid.N, nspin])
    for s in range(nspin):
        density_s[:, s] = rho / nspin
    return density_s
