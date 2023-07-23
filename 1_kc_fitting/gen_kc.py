'''
Generate KC data points from the fitted parameters
This is done after the fitting process (`kc.py` and `submit.py`) has been completed
'''
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import seaborn as sns

import kc

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def read_hdf(fn):
    with h5py.File(fn, 'r') as f:
        popt = f['popt'][()]
        pcov = f['pcov'][()]
    return popt, pcov

def gen_kc(dirname, freq=0.01):
    '''
    Evaluates the KC potential at four disregistries, i.e. 0.0, 0.16667, 0.5, 0.66667, and interlayer distance at every 0.01 ang.
    Creates the output data in `kc.csv` in the given `dirname` directory.
    '''
    workdir = os.getcwd()
    os.makedirs(dirname, exist_ok=True)
    os.chdir(dirname)

    if dirname == 'fit/Ouyang':
        popt = [3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588, 0]
    else:
        popt, pcov = read_hdf('result.hdf5')

    l = []
    for disregistry in [0.0, 0.16667, 0.5, 0.66667]:
        for distance in np.arange(2.8, 7.2 + freq, freq):
            distance = np.round(distance, 2)
            print(disregistry, distance)
            e = kc._eval_energy(distance, disregistry, *popt, keep_tmp_files=True)
            l.append({
                'label': dirname.replace('_', '-'),
                'disregistry': disregistry,
                'd': distance,
                'energy': e,
                'energy_inf': popt[-1],
                })
    k = pd.DataFrame(l)
    k['stacking'] = k['disregistry'].map({0.0: 'AB', 0.16667: 'SP', 0.5: 'Mid', 0.66667: 'AA'})
    k = k.sort_values(by=['label', 'disregistry', 'd'])
    k.to_csv('kc.csv', index=False)
    os.chdir(workdir)
    return k

def gen_kc_wrap(method, kT, freq=0.01):
    dirname = f'fit/{method}_kT{kT}' if method != 'Ouyang' else 'fit/Ouyang'
    csv_name = f'{dirname}/kc.csv'
    if os.path.isfile(csv_name):
        k = pd.read_csv(csv_name)
    else:
        k = gen_kc(dirname, freq=freq)
    k['method'] = f'KC-{method}'
    k['kT'] = int(float(kT)*1000) if kT not in ['inf', 'NA'] else kT
    return k

def gen_kc_all():
    '''
    Loops over different methods and kT to generate `kc.csv` within those directories.
    '''
    l = []
    for method in [
        'QMC',
        'Ouyang',
        'DFT-D2',
        'DFT-D3',
        'DFT-MBD'
        ]:
        if method == 'Ouyang':
            kTs = ['NA']
        else:
            kTs = [f'{kT:.3f}' for kT in np.arange(2, 21, 1)*0.001] + ['inf']
        for kT in kTs:
            l.append(gen_kc_wrap(method, kT))
    k = pd.concat(l, ignore_index=True)
    k = k.loc[:, ['method', 'kT', 'stacking', 'disregistry', 'd', 'energy', 'energy_inf']]
    print(k)
    os.makedirs('processed', exist_ok=True)
    k.to_csv('processed/kc.csv', index=False)

def print_params():
    '''
    Prints latex table
    '''
    l = []
    for method in ['QMC', 'Ouyang']:
        if method == 'Ouyang':
            popt = ['KC-Ouyang', 3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588, 0]
            l.append(popt)
        else:
            kTs = [f'{kT:.3f}' for kT in np.arange(2, 21, 1)*0.001] + ['inf']
            for kT in kTs:
                print(kT)
                fn = f'fit/{method}_kT{kT}/result.hdf5'
                popt, pcov = read_hdf(fn)
                kT_str = f'{int(float(kT)*1000)} meV' if kT != 'inf' else 'unweighted'
                l.append([f'KC-{method} ({kT_str})'] + popt.tolist())
    d = pd.DataFrame(l)
    d.columns = ['Potential', 'z0', 'C0', 'C2', 'C4', 'C', 'delta','lambda', 'A', 'Einf']
    d = d.drop(['Einf'], axis=1)
    print(d.round(5).to_latex(index=False))

if __name__ == '__main__':
    # gen_kc_all()
    print_params()
