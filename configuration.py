import dataclasses
import numpy as np
import orbital_type


@dataclasses.dataclass
class Level:
    nr: int  # radial quantum number
    l: int  # angular quantum number
    s: int  # spin quantum number
    eig: float  # eigen value
    occ: int  # occupation
    wfn: np.ndarray  # wave function (times radius)
    spin_polarized: bool

    @property
    def n(self):
        # principle quantum number
        return self.nr + self.l + 1

    @property
    def max_occ(self):
        if self.spin_polarized:
            return 2 * self.l + 1
        else:
            return 2 * (2 * self.l + 1)


def set_occupation(orbitals: list[Level], Z: int):
    # assume that the atoms try to maximize spin (Hund's rule)
    elec_remain = Z
    for orb in orbitals:
        if orb.max_occ >= elec_remain:
            orb.occ = elec_remain
            break
        else:
            orb.occ = orb.max_occ
            elec_remain -= orb.max_occ


def extract_configuration_from_str(config_str: str) -> dict[(int, int), int]:
    config_map = {}
    for config_str in config_str.split():
        n = int(config_str[0])
        l = orbital_type.orb_type_map_inv[config_str[1]]
        occ = int(config_str[2:])
        # sanity check
        if occ > 2 * (2 * l + 1):
            raise ValueError(f'wrong configuration: {config_str}!')
        config_map[(n, l)] = occ
    return config_map


def get_configuration_str_from_config_map(config_map: dict[(int, int), int]):
    config_str = ''
    for (n, l) in sorted(config_map, key=lambda x: x[0]):
        config_str += f'{n}{orbital_type.orb_type_map[l]}{config_map[(n, l)]} '
    if config_str != '':
        config_str = config_str[:-1]
    return config_str


def get_configuration(Z: int) -> dict[(int, int), int]:
    # for neutral atoms only
    # Aufbau principle from https://en.wikipedia.org/wiki/Aufbau_principle
    # also refer to https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Inorganic_Chemistry_(LibreTexts)/02%3A_Atomic_Structure/2.02%3A_The_Schrodinger_equation_particle_in_a_box_and_atomic_wavefunctions/2.2.03%3A_Aufbau_Principle
    config_map = {}
    nl = 1
    elec_remain = Z
    # 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p, 8s, 5g
    while elec_remain > 0:
        for nr in range(nl):
            l = nl - nr - 1
            if l <= nr:
                max_occ = 2 * (2 * l + 1)
                if elec_remain > max_occ:
                    occ = max_occ
                else:
                    occ = elec_remain
                if occ == 0:
                    break
                config_map[(nr+1, l)] = occ
                elec_remain -= occ
        nl += 1
    # exceptions in d-block
    if Z in [24, 29]:  # Cr, Cu
        config_map[(3, 2)] += 1
        config_map[(4, 0)] -= 1
    elif Z in [41, 42, 44, 45, 47]:  # Nb, Mo, Ru, Rh, Ag
        config_map[(4, 2)] += 1
        config_map[(5, 0)] -= 1
    elif Z in [46]:  # Pd
        config_map[(4, 2)] += 2
        del config_map[(5, 0)]
    elif Z in [78, 79]:  # Pt, Au
        config_map[(5, 2)] += 1
        config_map[(6, 0)] -= 1
    elif Z in [103]:  # Lr
        del config_map[(6, 2)]
        config_map[(7, 1)] = 1
    # exceptions in f-block
    elif Z in [57]:  # La
        del config_map[(4, 3)]
        config_map[(5, 2)] = 1
    elif Z in [58, 64]:  # Ce, Gd
        config_map[(4, 3)] -= 1
        config_map[(5, 2)] = 1
    elif Z in [89, 90]:  # Ac, Th
        config_map[(6, 2)] = config_map[(5, 3)]
        del config_map[(5, 3)]
    elif Z in [91, 92, 93, 96]:  # Pa, U, Np, Cm
        config_map[(5, 3)] -= 1
        config_map[(6, 2)] = 1
    return config_map


def set_fixed_occupation(orbitals: list[Level], config_map: dict[(int, int), int], spin_polarized: bool):
    # assume that the atoms try to maximize spin (Hund's rule)
    for orb in orbitals:
        n = orb.n
        l = orb.l
        if (n, l) in config_map:
            if spin_polarized:
                if config_map[(n, l)] <= orb.max_occ:
                    if orb.s == 0:
                        orb.occ = config_map[(n, l)]
                    else:
                        orb.occ = 0
                else:
                    if orb.s == 0:
                        orb.occ = orb.max_occ
                    else:
                        orb.occ = config_map[(n, l)] - orb.max_occ
            else:
                orb.occ = config_map[(n, l)]


def get_nr_max_for_l(config_map: dict[(int, int), int]) -> dict[int, int]:
    nr_max_l_map = {}
    for (n, l) in config_map:
        nr = n - l - 1
        if l in nr_max_l_map:
            if nr_max_l_map[l] < nr:
                nr_max_l_map[l] = nr
        else:
            nr_max_l_map[l] = nr
    return nr_max_l_map


def get_orbital_name(orb: Level, spin_polarized: bool, latex: bool = False):
    if latex:
        if spin_polarized:
            spin_str = '\\uparrow$' if orb.s == 0 else '\\downarrow$'
        else:
            spin_str = '$'
        return f'${orb.n}{orbital_type.orb_type_map[orb.l]}^' + '{' + f'{orb.occ}' + '}' + spin_str
    else:
        if spin_polarized:
            spin_str = ' u' if orb.s == 0 else ' d'
        else:
            spin_str = ''
        return f'{orb.n}{orbital_type.orb_type_map[orb.l]}{orb.occ}' + spin_str


def print_orbitals(orbitals: list[Level], spin_polarized: bool, unoccupied: bool = False):
    for orb in orbitals:
        if orb.occ > 0 or unoccupied:
            print(f'{get_orbital_name(orb, spin_polarized)} {orb.eig}')
