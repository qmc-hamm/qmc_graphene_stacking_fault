'''
Evaluated KC potential energy from multiple models for evaluation of quantities by `bootstrap_kc_find_min.py`
'''
import h5py
import pandas as pd
import numpy as np
import os

import kc

def gen_kc(dirname, freq=0.01):
    '''
    Evaluates the KC potential at four disregistries, i.e. 0.0, 0.16667, 0.5, 0.66667, and interlayer distance at every 0.01 ang.
    Creates the output data in `kc.csv` in the given `dirname` directory.
    '''
    workdir = os.getcwd()
    os.makedirs(dirname, exist_ok=True)
    os.chdir(dirname)

    with h5py.File('result.hdf5', 'r') as f:
        popt = f['popt'][()]

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

if __name__ == '__main__':
    for model_id in range(0, 20):
        dirname = f'fit_bootstrap/QMC_{model_id:02}_kTinf'
        gen_kc(dirname)
