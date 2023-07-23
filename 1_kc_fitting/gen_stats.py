'''
Find RMS of fitting in the specified range `lim` = `min` or `far`
'''
import h5py
import json
import numpy as np
import pandas as pd
import os

import kc, read

def read_hdf(fn):
    with h5py.File(fn, 'r') as f:
        popt = f['popt'][()]
        pcov = f['pcov'][()]
    return popt, pcov

def get_dirname(method, kT):
    kT_str = f'{kT*0.001:.3f}' if type(kT) != str else kT
    return f'fit/{method}_kT{kT_str}'

def get_d_range(lim):
    d_range_map = {
        'min': (3.1, 3.9),
        'far': (2.9, 7.1)
        }
    return d_range_map[lim]

def get_r2_rms(df, popt, lim):
    d_range = get_d_range(lim)
    df = df.loc[(df.d > d_range[0]) & (df.d < d_range[1]), :].reset_index(0, drop=True)
    ydata = df['energy']

    residuals = ydata - kc.eval_energy(df, *popt, keep_tmp_files=False)
    ss_res = np.sum(residuals**2)
    rms = (ss_res / len(ydata))**0.5
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return r2, rms

def write_json(fn, d):
    with open(fn, 'w') as f:
        f.write(json.dumps(d, indent=4))

def get_stats(method, kT, overwrite=False):
    dirname = get_dirname(method, kT)
    json_file = f'{dirname}/stats.json'
    if os.path.isfile(json_file) and not overwrite:
        with open(json_file) as f:
            stats = json.load(f)
    else:
        df_raw = read.read_data(method)
        popt, pcov = read_hdf(f'{dirname}/result.hdf5')
        stats = []
        for lim in ['min', 'far']:
            r2, rms = get_r2_rms(df_raw, popt, lim)
            stat = {
                'method': method,
                'kT': kT,
                'lim': lim,
                'r2': r2,
                'rms': rms
                }
            print(stat)
            stats.append(stat)
        write_json(json_file, stats)
    return stats

def loop_dirs(methods, kTs):
    for method in methods:
        for kT in kTs:
            yield method, kT

def collect_stats():
    methods = ['QMC', 'DFT-D2', 'DFT-D3', 'DFT-MBD']
    kTs = np.arange(2, 21, 1).tolist() + ['inf']
    l = []
    for method, kT in loop_dirs(methods, kTs):
        stats = get_stats(method, kT, overwrite=False)
        l += stats
    d = pd.DataFrame(l)
    os.makedirs('processed', exist_ok=True)
    d.to_csv('processed/stats.csv', index=False)

if __name__ == '__main__':
    collect_stats()
