'''
Evalutates the minimum of the energy curve by quadratic fitting around the minimum
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

import read

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def quad(x, c2, c1, c0):
    return c2*x**2 + c1*x + c0

def get_quad_min(coeffs):
    return -coeffs[1]/(2*coeffs[0])

def find_r2(ydata, ypredicted):
    residuals = ydata - ypredicted
    ss_res = np.sum(residuals**2)
    rms = (ss_res / len(ydata))**0.5
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return r2, rms

def find_min_xy(x, y):
    coeffs = np.polyfit(x, y, 2)
    r2, rms = find_r2(y, quad(x, *coeffs))
    npoints = len(y)
    d_eq_quad = get_quad_min(coeffs)
    be_quad = quad(d_eq_quad, *coeffs)
    return d_eq_quad, be_quad, coeffs, r2, npoints

def fit_within(df, d_min, d_max):
    mask = (d_min <= df.d) & (df.d <= d_max)
    ddf = df.loc[mask, :]
    d_eq_quad, be_quad, coeffs, r2, npoints = find_min_xy(ddf.d, ddf.energy)
    return d_eq_quad, be_quad, coeffs, r2, npoints

def plot_ranges(df):
    method = df.method.values[0]
    stacking = df.stacking.values[0]
    l = []
    for d_range_tmp in np.arange(0.02, 0.10+0.001, 0.01):
        d_eq = df.loc[df.energy == df.energy.min(), 'd'].values[0]
        d_eq_quad, be_quad, coeffs, r2, npoints = fit_within(df, d_eq-d_range_tmp, d_eq+d_range_tmp)
        l.append({'d_range': d_range_tmp, 'd_eq': d_eq_quad, 'be': be_quad, 'r2': r2, 'npoints': npoints})
    eq = pd.DataFrame(l)
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3, 6), sharex=True)
    axs[0].plot(eq.d_range, eq.d_eq, 'o')
    axs[0].set(ylabel='deq (ang)')
    axs[1].plot(eq.d_range, eq.r2, 'o')
    axs[1].set(ylabel='R2')
    axs[2].plot(eq.d_range, eq.npoints, 'o')
    axs[2].set(ylabel='npoints')
    axs[2].set(xlabel='range (ang)')
    os.makedirs('min_plots', exist_ok=True)
    plt.savefig(f'min_plots/{method}_{stacking}_range.pdf', bbox_inches='tight')
    plt.close()

def plot_quad(df, d_range):
    method = df.method.values[0]
    stacking = df.stacking.values[0]
    d_eq = df.loc[df.energy == df.energy.min(), 'd'].values[0]
    d_eq_quad, be_quad, coeffs, r2, npoints = fit_within(df, d_eq-d_range, d_eq+d_range)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    ax.plot(df.d, df.energy, 'o', ms=5, label=f'{method}: deq={d_eq:.4f}')

    # quad
    lin = np.linspace(d_eq-d_range, d_eq+d_range, 500)
    ax.plot(lin, quad(lin, *coeffs), label=f'quadratic: deq={d_eq_quad:.4f}')
    ax.axvline(x=d_eq_quad)
    ax.set_xlim(3.0, 3.9)
    ax.set_ylabel('Binding energy (eV/atom)')
    ax.set_xlabel('d (ang)')
    ax.legend(frameon=False)
    plt.savefig(f'min_plots/{method}_{stacking}_quad.pdf', bbox_inches='tight')
    plt.close()

def find_min(df, d_range=0.05, gen_plot_quad=False, gen_plot_ranges=False):
    d_eq = df.loc[df.energy == df.energy.min(), 'd'].values[0]
    d_eq_quad, be_quad, coeffs, r2, npoints = fit_within(df, d_eq-d_range, d_eq+d_range)
    if gen_plot_ranges:
        plot_ranges(df)
    if gen_plot_quad:
        plot_quad(df, d_range=d_range)
    return d_eq_quad, be_quad

def test_find_min():
    d = read.read_kc_best()
    dd = d.loc[(d.method == 'KC-QMC') & (d.stacking == 'AB'), :]
    print(find_min(dd, gen_plot_quad=True, gen_plot_ranges=True))

def gen_min(lim='min'):
    d = read.read_kc_best(align='inf', lim=lim)
    l = []
    for method in d.method.unique():
        dd_ab = d.loc[(d.method == method) & (d.stacking == 'AB'), :]
        d_eq_ab, e_min_ab = find_min(dd_ab)
        for stacking in d.stacking.unique():
            dd = d.loc[(d.method == method) & (d.stacking == stacking), :]
            d_eq, e_min = find_min(dd, gen_plot_quad=True, gen_plot_ranges=True)
            l.append({
                'method': method,
                'stacking': stacking,
                'd_eq': d_eq,
                'gsfe': e_min - e_min_ab,
                'be': -e_min
                })
    m = pd.DataFrame(l)
    print(m)
    m.to_csv(f'processed/{lim}.csv', index=False)

if __name__ == '__main__':
    gen_min(lim='far')
