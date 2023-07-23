'''
Reads in raw data from calculations (`dft.csv` or `qmc.csv`) or the KC-fitted data
'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.optimize
import seaborn as sns

def read_dft(label):
    d = pd.read_csv('../0_interlayer_energy/data/dft.csv')

    d = d.loc[d['vdw_corr'] != 'none', :]
    d['label'] = d['vdw_corr'].map({'none': 'PBE', 'dft-d': 'DFT-D2', 'dft-d3': 'DFT-D3', 'ts': 'DFT-TS'})
    d['energy_err'] = 0
    d = d[['stacking', 'disregistry', 'label', 'd', 'energy', 'energy_err']]
    dd = d.loc[d.label == label, :].reset_index(0, drop=True)
    return dd

def read_qmc():
    qmc = pd.read_csv('../0_interlayer_energy/data/qmc.csv')
    return qmc

def read_data(method):
    '''
    Wrapper function of data reading functions
    '''
    if 'DFT' in method or 'PBE' in method:
        d = read_dft(method)
    elif method == 'QMC':
        d = read_qmc()
    d['method'] = method
    return d

def morse(r, D, r0, a, E0):
    '''
    Morse potential.
    '''
    return D*(1 - np.exp(-a*(r - r0)))**2 + E0

def morse_fit(x, y):
    '''
    Performs a morse fit from given `x` and `y`.
    Returns a list of the fitting parameters: D, r0, a, E0
    '''

    # find a starting guess value `D` for the morse fit by doing a quadratic fit near the minimum
    d_eq_idx = y.idxmin()
    d_eq = x[d_eq_idx]
    d_range = 1
    mask = (d_eq - d_range < x) & (x < d_eq + d_range)
    xx = x[mask]
    yy = y[mask]
    fit = np.polyfit(xx, yy, 2)

    # the starting guess for D is D0 = fit[0]
    # this is because at near minimum r = r0 + xi, and xi << 1, y(xi) = D*a**2*xi**2 + E0
    p0 = [fit[0], 3.4, 1.0, np.min(y)]
    popt, pcov = scipy.optimize.curve_fit(morse, x, y, p0)
    return popt

def get_energy_inf(d, gen_plot=False):
    '''
    Find a rough estimate of the energy at infinite interlayer spacing by fitting to a morse function
    To be used as a starting guess for KC fitting
    '''
    if gen_plot:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        palette = sns.color_palette()
        for i, (stacking, g) in enumerate(d.groupby('stacking')):
            x = g.d
            y = g.energy
            popt = morse_fit(x, y)
            ax.plot(x, y, 'o', label=stacking, color=palette[i])
            lin = np.linspace(np.min(x), np.max(x), 500)
            ax.plot(lin, morse(lin, *popt), '-', label=None, color=palette[i])
        ax.legend()
        plt.savefig('fit_morse.pdf', bbox_inches='tight')
        plt.close()

    energy_inf_l = []
    for i, (stacking, g) in enumerate(d.groupby('stacking')):
        x = g.d
        y = g.energy
        popt = morse_fit(x, y)
        energy_inf_l.append(popt[0] + popt[3])
    energy_inf = np.mean(energy_inf_l)
    return energy_inf

def read_kc(method, kTs=[]):
    '''
    Retrieves the KC data points from the given method and kTs (fitting weights).
    '''
    k = pd.read_csv('processed/kc.csv')
    k = k.replace([np.inf], 'inf')
    if method == 'KC-Ouyang':
        mask = (k.method == method)
    else:
        mask = (k.method == method) & (k.kT.isin(kTs))
    kk = k.loc[mask, :]
    return kk

def get_best_kTs(method, lim='min'):
    '''
    Provides a mapping for which value of kT is the best choice for the method and regions of interest (`min` or `far`)
    `min`: interlayer distance from 3.2 to 3.8
    `far`: all interlayer distances
    '''
    if lim == 'min':
        best_kTs = {
            'KC-QMC': 4,
            'KC-DFT-D2': 3,
            'KC-DFT-D3': 20
        }
    elif lim == 'far':
        best_kTs = {
            'KC-QMC': 'inf',
            'KC-DFT-D2': 'inf',
            'KC-DFT-D3': 'inf',
        }
    return best_kTs[method]

def read_kc_best(align=None, lim='min'):
    '''
    Reads the KC data points from the best fitting weights for the regions of interest (`min` or `far).
    '''
    d = pd.concat([
        read_kc('KC-QMC', kTs=[get_best_kTs('KC-QMC', lim=lim)]),
        read_kc('KC-Ouyang', kTs=[]),
        read_kc('KC-DFT-D2', kTs=[get_best_kTs('KC-DFT-D2', lim=lim)]),
        read_kc('KC-DFT-D3', kTs=[get_best_kTs('KC-DFT-D3', lim=lim)]),
        ], ignore_index=True)
    if align == 'inf':
        d.energy -= d.energy_inf
    elif align == 'min':
        d.energy -= d.groupby(['method'])['energy'].transform(min)

    d.energy *= 1000
    d.energy_inf *= 1000
    d = d[['method', 'stacking', 'disregistry', 'd', 'energy', 'energy_inf']]
    return d

def test_read_kc():
    k = read_kc('KC-QMC', kTs=[2, 4, 'inf'])
    print(k)
    k = read_kc('KC-Ouyang', kTs=[])
    print(k)

def test_read_data():
    methods = ['DFT-D2', 'DFT-D3', 'QMC']
    for method in methods:
        d = read_data(method)
        energy_inf = get_energy_inf(d)
        print(energy_inf)

if __name__ == '__main__':
    print(read_kc_best(align='inf', lim='min'))
    print(read_kc_best(align='inf', lim='far'))
