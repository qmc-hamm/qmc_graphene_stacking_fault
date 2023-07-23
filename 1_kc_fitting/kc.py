'''
Perform KC fitting to the raw data from either DFT or QMC

Author: Kittithat Krongchon, Lucas K. Wagner
'''
import h5py
import numpy as np
import numpy.linalg as la
import os
import pandas as pd
import scipy.optimize
import shutil
import argparse

from ase.calculators.lammpsrun import LAMMPS
import gen_geom
import read

def get_params_hybrid():
    '''
    Defines REBO and KC potentials
    '''
    parameters = {
        'pair_style': 'hybrid/overlay rebo kolmogorov/crespi/full 16.0 1',
        'pair_coeff': [
            '* * rebo CH.rebo C C',
            '* * kolmogorov/crespi/full CH_taper.KC C C'
            ],
        'atom_style': 'full',
        'specorder': ['C', 'C'],
        }
    files = ['CH.rebo', 'CH_taper.KC']
    return parameters, files

def get_params_kc_only():
    '''
    Defines only the KC potential
    '''
    parameters = {
        'pair_style': 'hybrid/overlay kolmogorov/crespi/full 16.0 1',
        'pair_coeff': [
            '* * kolmogorov/crespi/full CH_taper.KC C C'
            ],
        'atom_style': 'full',
        'specorder': ['C', 'C'],
        }
    files = ['CH_taper.KC']
    return parameters, files

def format_params(params, sep=' ', prec='.15f'):
    '''
    Returns a string joining KC parameters
    '''
    l = [f'{param:{prec}}' for param in params]
    s = sep.join(l)
    return s

def write_kc_potential(z0, C0, C2, C4, C, delta, lamda, A, kc_filename='CH_taper.KC', cite=False):
    params = [z0, C0, C2, C4, C, delta, lamda, A]
    headers = '               '.join(['', 'z0', 'C0', 'C2', 'C4', 'C', 'delta', 'lambda', 'A'])
    with open(kc_filename, 'w') as f:
        lines = ['# Refined parameters for Kolmogorov-Crespi Potential with taper function', '#']
        if cite:
            lines += ['# Cite as Krongchon, K., Rakib, T., Pathak, S., Ertekin, E., Johnson, H. T., & Wagner, L. K. (2023). arXiv preprint arXiv:2307.07210.', '#']
        lines += [f'# {headers}         S     rcut', f'C C {format_params(params)} 1.0    2.0']
        f.write('\n'.join(lines))

def _eval_energy(distance, disregistry, z0, C0, C2, C4, C, delta, lamda, A, E0, keep_tmp_files=True, tmp_dir='lmp_tmp'):
    '''
    Finds the energy for a given geometry (defined by `distance` and `disregistry`) and KC paramters
    '''
    atoms = gen_geom.create_graphene_geom(distance, disregistry)
    atoms.set_array('mol-id', np.array([0, 0, 1, 1]))

    write_kc_potential(z0, C0, C2, C4, C, delta, lamda, A)
    parameters, files = get_params_kc_only()
    # lammps_options = '' if keep_tmp_files else '-echo log -screen none -log /dev/stdout'
    atoms.calc = LAMMPS(files=files, keep_tmp_files=True, tmp_dir=tmp_dir, lammps_options='', **parameters)
    e = atoms.get_potential_energy()/len(atoms) + E0
    if not keep_tmp_files:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        for fn in files:
            if os.path.isfile(fn):
                os.remove(fn)
        lmp_output = 'log.lammps'
        if os.path.isfile(lmp_output):
            os.remove(lmp_output)
    return e

def eval_energy(df, z0, C0, C2, C4, C, delta, lamda, A, E0, keep_tmp_files=True):
    '''
    Finds the energy for a given geometry (in `df`) and KC paramters
    '''
    tmp_dir = 'lmp_tmp'
    shutil.rmtree(tmp_dir, ignore_errors=True)

    energy = []
    for _, row in df.iterrows():
        e = _eval_energy(row['d'], row['disregistry'], z0, C0, C2, C4, C, delta, lamda, A, E0, keep_tmp_files=keep_tmp_files, tmp_dir=tmp_dir)
        energy.append(e)

    en =  np.array(energy)
    return en

en_prev = 0
def eval_energy_track(df, z0, C0, C2, C4, C, delta, lamda, A, E0):
    '''
    Finds the energy for a given geometry (in `df`) and KC paramters
    Also keep track of the energy in the previous iteration
    '''
    global en_prev
    en = eval_energy(df, z0, C0, C2, C4, C, delta, lamda, A, E0, keep_tmp_files=True)
    en_diff = la.norm(en - en_prev)
    en_prev = en.copy()
    params = [z0, C0, C2, C4, C, delta, lamda, A, E0]
    print(f'{en_diff: 25.15f} ' + format_params(params, prec=' 18.15f'))
    return en

def fit(method='QMC', kT='inf', model_id=0, starting_guess='old_fit'):
    '''
    `kT`: temperature in the boltzmann factor, used to assign weights
    bound all params = [0, np.inf]
    '''
    now = datetime.datetime.now()
    print(now)
    df = read.read_data(method)
    ydata = np.random.normal(loc=df['energy'], scale=df['energy_err'])

    dir_prefix = f'fit_bootstrap/{method}_{int(model_id):02}_kT'
    if kT == 'inf':
        kT_str = 'inf'
        weights = np.nan
        sigma = None
    else:
        # set weights such that points near the minimum (lower energy) have more weights, e.g. kT=0.002, 0.004, 0.010
        kT = float(kT)
        kT_str = f'{kT:.3f}'
        weights = np.exp(-(ydata-ydata.min())/kT)
        sigma = 1/np.sqrt(weights)
    dirname = dir_prefix + kT_str
    workdir = os.getcwd()
    os.makedirs(dirname, exist_ok=True)
    os.chdir(dirname)

    if starting_guess == 'ouyang':
        print('starting guess from Ouyang, Mandelli, Urbakh, and Hod, Nanoserpents: Graphene Nanoribbon Motion on Two-Dimensional Hexagonal Materials, Nano letters, 2018')
        energy_inf = read.get_energy_inf(df)
        p0 = [3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588, energy_inf]
    elif starting_guess == 'old_fit':
        print('use a starting guess from the old fit')
        hdf_filename = os.path.join(workdir, f'fit/{method}_kT{kT_str}/result.hdf5')
        with h5py.File(fn, 'r') as f:
            p0 = f['popt'][()]
        print('starting guess parameters: ', p0)
        energy_inf = p0[-1]

    energy_range = 0.01
    popt, pcov = scipy.optimize.curve_fit(eval_energy_track, df, ydata, p0=p0, method='trf', sigma=sigma,
        bounds = (
            [0, 0, 0, 0, 0, 0, 0, 0, energy_inf - energy_range],
            [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, energy_inf + energy_range]
            )
        )

    with h5py.File('result.hdf5', 'w') as f:
        f['popt'] = popt
        f['pcov'] = pcov
        f['weights'] = weights

    os.chdir(workdir)
    now = datetime.datetime.now()
    print(now)

if __name__ == '__main__':
    import datetime
    parser = argparse.ArgumentParser()
    parser.add_argument('-kT', default='inf')
    parser.add_argument('--method', default='QMC')
    parser.add_argument('--model_id', default=0)
    args = parser.parse_args()
    fit(method=args.method, kT=args.kT, model_id=args.model_id)
