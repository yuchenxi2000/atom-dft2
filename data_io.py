import configuration
import numpy as np
import re


def write_orbitals_to_dir(orbitals: list[configuration.Orbital], out_dir) -> None:
    with open(f'{out_dir}/orbital.txt', 'w') as file_orb:
        for orb in orbitals:
            file_orb.write(f'{orb.get_str(include_occ=True, latex=False)} {orb.eig}\n')
    for orb in orbitals:
        orb.wfn.tofile(f'{out_dir}/wfn_{orb.get_str(include_occ=False, latex=False)}.bin')


def read_orbitals_from_dir(calc_dir) -> list[configuration.Orbital]:
    orbitals = []
    with open(f'{calc_dir}/orbital.txt', 'r') as file_orb:
        for line in file_orb:
            data = line.split()
            match_res = re.match(r'([0-9]+)([a-z])([0-9]*)(.*)', data[0])
            eig = float(data[1])
            orb_str_no_occ = match_res.group(1) + match_res.group(2) + match_res.group(4)
            wfn_file = f'{calc_dir}/wfn_{orb_str_no_occ}.bin'
            wfn = np.fromfile(wfn_file, dtype=float, count=-1)
            wfn = wfn.reshape([-1, 2])
            orbitals.append(configuration.Orbital.from_str(data[0], eig, wfn))
    return orbitals
