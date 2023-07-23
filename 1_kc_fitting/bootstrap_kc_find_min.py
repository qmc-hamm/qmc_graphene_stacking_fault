'''
Finds the statistical errors of quantities such as minimum interlayer spacing or binding energy
by fitting multiple models to data points, each of which is a gaussian variable with standard deviation equal to the QMC error bar
'''
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def read_hdf(fn):
    with h5py.File(fn, 'r') as f:
        popt = f['popt'][()]
        pcov = f['pcov'][()]
    return popt, pcov


def quad(x, c2, c1, c0):
    return c2*x**2 + c1*x + c0

def get_quad_min(coeffs):
    return -coeffs[1]/(2*coeffs[0])

def find_min_xy(x, y):
    coeffs = np.polyfit(x, y, 2)
    d_eq_quad = get_quad_min(coeffs)
    be_quad = quad(d_eq_quad, *coeffs)
    return d_eq_quad, be_quad, coeffs

def find_min(g):
    d_range = 0.05
    d_eq = g.loc[g.energy == g.energy.min(), 'd'].values[0]
    mask = (d_eq - d_range <= g.d) & (g.d <= d_eq + d_range)
    gg = g.loc[mask, :]
    d_eq_quad, be_quad, coeffs = find_min_xy(gg['d'], gg['energy'])
    return d_eq_quad, be_quad, coeffs

def print_aggregated_dmin():
    q = pd.read_csv('data/qmc.csv')
    fig, axs = plt.subplots(nrows=1, ncols=4, sharey=True, sharex=True)
    dmin = []
    for model_id in range(20):
        csvname = f'fit_bootstrap/QMC_{model_id:02}_kT0.004/kc.csv'
        d = pd.read_csv(csvname)
        d = d.drop(['disregistry', 'label'], axis=1)
        for i, stacking in enumerate(d['stacking'].unique()):
            ax = axs[i]
            dd = d.loc[d['stacking'] == stacking, :]
            qq = q.loc[q['stacking'] == stacking, :]
            d_eq_quad, be_quad, coeffs = find_min(dd)
            dmin.append({
                'model_id': model_id,
                'stacking': stacking,
                'dmin': d_eq_quad
                })

    a = pd.DataFrame(dmin)
    print(a.groupby('stacking')['dmin'].apply(np.std))
    print(a.groupby('stacking')['dmin'].apply(np.average))

def print_aggregated_be():
    be_list = []
    for model_id in range(20):
        csvname = f'fit_bootstrap/QMC_{model_id:02}_kTinf/kc.csv'
        d = pd.read_csv(csvname)
        d = d.drop(['disregistry', 'label'], axis=1)
        energy_inf = d['energy_inf'].values[0]
        for i, stacking in enumerate(d['stacking'].unique()):
            dd = d.loc[d['stacking'] == stacking, :]
            d_eq_quad, be_quad, coeffs = find_min(dd)
            be = energy_inf - be_quad
            be_list.append({
                'model_id': model_id,
                'stacking': stacking,
                'be': be
                })
    a = pd.DataFrame(be_list)
    print(a.groupby('stacking')['be'].apply(np.std)*1000)
    print(a.groupby('stacking')['be'].apply(np.average)*1000)

if __name__ == '__main__':
    print_aggregated_dmin()
    print_aggregated_be()
