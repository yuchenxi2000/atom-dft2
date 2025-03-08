"""
Density functional theory of atoms (v2.0)
Copyright (c) 2024 YCX. Licensed under GPL (v3) License.

Local density approximation, non-relativistic calculation, spin unpolarized & polarized
Can quantitatively match the results on website:
https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7

Acknowledgements:
Thanks to the project of aromanro https://compphys.go.ro/dft-for-an-atom/
and the article https://arxiv.org/pdf/1209.1752v2
"""
import matplotlib.pyplot as plt
import numpy as np
import pylibxc
import configuration
import density
import energy
import poisson
import radial_wave
import xc_func
import tomllib
import argparse
import periodic_table
import data_io

# parse arguments
parser = argparse.ArgumentParser(description='Density functional theory of isolated atoms (v2.0), author: YCX')
parser.add_argument('-c', '--config', default='atom_config.toml', help='config file in toml format')
args = parser.parse_args()

# read config file in toml format
config = tomllib.load(open(args.config, 'rb'))

N = config['params']['grid']['num_grids']
r0 = config['params']['grid']['r0']
rc = config['params']['rcut']
grid = radial_wave.radial_grid.exp_grid(N, r0, rc)
# core charge
if 'element' in config['params']:
    element = config['params']['element']
    Z = periodic_table.periodic_table.index(element)
elif 'core_charge' in config['params']:
    Z = config['params']['core_charge']
    element = periodic_table.periodic_table[Z]
elif 'atomic_number' in config['params']:
    Z = config['params']['atomic_number']
    element = periodic_table.periodic_table[Z]
else:
    raise ValueError('please provide atomic number, core charge or element name!')
r = grid.r
nspin = 2 if config['params']['spin_polarized'] else 1

# use C dylib to accelerate calculation
use_c = config['params']['use_c']

# exchange-correlation functional, LDA with Slater exchange and VWN correlation
spin_description = 'polarized' if nspin == 2 else 'unpolarized'
x_func = pylibxc.LibXCFunctional('LDA_X', spin_description)
c_func = pylibxc.LibXCFunctional('LDA_C_VWN', spin_description)

# initial density
density_s = density.get_approx_density_thomas_fermi(grid, Z, nspin)

# electron configuration
if config['occupation']['type'] == 'aufbau':
    config_map = configuration.get_configuration(Z)
elif config['occupation']['type'] == 'fixed':
    config_map = configuration.extract_configuration_from_str(config['occupation']['configuration'])
else:
    raise Exception(f'unknown occupation type {config['occupation']['type']}')
nr_max_l_map = configuration.get_nr_max_for_l(config_map)

# initial Vext
VH, dVH = poisson.solve_poisson(np.sum(density_s, axis=1), grid, use_c=use_c)
Vxc, exc = xc_func.get_Vxc_exc_lda(density_s, x_func, c_func)
Vext = np.expand_dims(VH, axis=1) + Vxc

orbitals = []
step = 0
max_step = config['calculation']['max_steps']
Eprev = 0.0
history_Vext = []
history_residue = []
history_max_num = config['density']['mixing']['history_max_num']
energy_eps = config['calculation']['energy_eps']
Vext_eps = config['calculation']['vext_eps']
mixing_beta = config['density']['mixing']['mixing']
while step < max_step:
    # get orbitals
    orbitals = configuration.get_orbitals(Z, nspin, Vext, grid, nr_max_l_map, use_c=use_c)
    # set occupation
    configuration.set_fixed_occupation(orbitals, config_map, nspin == 2)
    # get new density
    density_s = density.get_density_r2s(orbitals, N, nspin) / np.expand_dims(r, axis=1) ** 2
    # get VH
    VH, dVH = poisson.solve_poisson(np.sum(density_s, axis=1), grid, use_c=use_c)
    # get Vxc, exc
    Vxc, exc = xc_func.get_Vxc_exc_lda(density_s, x_func, c_func)
    # Vext
    Vext_new = np.expand_dims(VH, axis=1) + Vxc
    # get energy
    Etot = energy.get_total_energy(orbitals, exc, Vxc, VH, density_s, Z, grid).Etot
    # converged?
    norm = np.linalg.norm(Vext_new - Vext)
    if abs(Etot - Eprev) < energy_eps and norm < Vext_eps:
        break
    else:
        Eprev = Etot
    # record history
    history_Vext.append(Vext_new)
    history_residue.append(Vext_new - Vext)
    if len(history_Vext) > history_max_num:
        history_Vext.pop(0)
        history_residue.pop(0)
    # Vext mixing
    Vext = density.density_mixing_anderson(history_Vext, history_residue, mixing_beta)
    # print info
    print(f'step = {step}, energy = {Eprev}, Vext residue = {norm}')
    # next step
    step += 1

# print energy
if config['output']['print']['energy']:
    energy_parts = energy.get_total_energy(orbitals, exc, Vxc, VH, density_s, Z, grid)
    print(f'Etot = {energy_parts.Etot}')
    print(f'Ekin = {energy_parts.Ekin}')
    print(f'Ecoul = {energy_parts.EH}')
    print(f'Eenuc = {energy_parts.Eenuc}')
    print(f'Exc = {energy_parts.Exc}')
# print levels
if config['output']['print']['configuration']:
    orbitals.sort(key=lambda lvl: lvl.eig)
    for orb in orbitals:
        if orb.occ > 0:
            print(f'{orb.get_str(include_occ=True, latex=False)} {orb.eig}')

# output files
if 'dir' in config['output']:
    out_dir = config['output']['dir']
    density_s.tofile(f'{out_dir}/density_s.bin')
    exc.tofile(f'{out_dir}/exc.bin')
    grid.r.tofile(f'{out_dir}/grid.bin')
    VH.tofile(f'{out_dir}/VH.bin')
    dVH.tofile(f'{out_dir}/dVH.bin')
    Vxc.tofile(f'{out_dir}/Vxc.bin')
    data_io.write_orbitals_to_dir(orbitals, out_dir)
    with open(f'{out_dir}/energy.txt', 'w') as fenergy:
        fenergy.write(f'spin = {nspin}\n')
        fenergy.write(f'Z = {Z}\n')
        fenergy.write(f'element = {element}\n')
        fenergy.write(f'configuration = {configuration.get_configuration_str_from_config_map(config_map)}\n')
        fenergy.write(f'steps = {step}\n')
        energy_parts = energy.get_total_energy(orbitals, exc, Vxc, VH, density_s, Z, grid)
        fenergy.write(f'Etot = {energy_parts.Etot}\n')
        fenergy.write(f'Ekin = {energy_parts.Ekin}\n')
        fenergy.write(f'Ecoul = {energy_parts.EH}\n')
        fenergy.write(f'Eenuc = {energy_parts.Eenuc}\n')
        fenergy.write(f'Exc = {energy_parts.Exc}\n')
        fenergy.write(f'eig_sum = {energy_parts.eig_sum}\n')

# plot
if config['output']['plot']['show']:
    use_latex = config['output']['plot']['use_latex']
    plt.suptitle(f'Atom {element}')
    # density (r2 n)
    plt.subplot(211)
    plt.title('Density ($r^2n$)' if use_latex else 'Density (r2 n)')
    if nspin == 2:
        plt.plot(r, density_s[:, 0] * r ** 2, label='Spin up')
        plt.plot(r, density_s[:, 1] * r ** 2, label='Spin down')
        plt.legend()
    else:
        plt.plot(r, density_s[:, 0] * r ** 2)
    plt.xlim([-1, 9])
    # wave (r phi)
    plt.subplot(212)
    plt.title('Wave function ($r\\psi$)' if use_latex else 'Wave function (r psi)')
    for orb in orbitals:
        if orb.occ > 0:
            if config['output']['plot']['show_orbitals'] == 'all' or orb.eig > -1.837:
                # the criterion for outer shells (eig larger than -1.837 Hatree or -50.0 eV) is not precise.
                # may be improved in the future.
                plt.plot(r, orb.wfn[:, 0], label=orb.get_str(include_occ=True, latex=use_latex))
    plt.ylim([-0.4, 0.4])
    plt.xlim([-1, 9])
    plt.legend()
    # tight layout
    plt.tight_layout()
    plt.show()
