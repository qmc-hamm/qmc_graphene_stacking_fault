'''
Plot potential energy surface of different methods
'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import seaborn as sns

import read

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_clean_ticks(vmin, vmax, freq):
    '''
    customize axis ticks so that the ticks are at the given frequency
    e.g. vmin = 3.15, vmax = 3.85, freq = 0.2 will result in ticks at 3.2, 3.4, 3.6, 3.8
    '''
    tick_min = np.ceil(vmin/freq)*freq
    tick_max = np.floor(vmax/freq)*freq
    ticks = np.arange(tick_min, tick_max+freq, freq)
    return ticks

def custom_ticks(axs, xmin, xmax, ymin, ymax, xfreq, yfreq):
    '''
    Sets the tick min max and frequency on the given axes
    '''
    for ax in axs:
        xticks = get_clean_ticks(xmin, xmax, xfreq)
        yticks = get_clean_ticks(ymin, ymax, yfreq)
        ax.set(xticks=xticks, yticks=yticks)

def plot_near_eq(xlim=(3.15, 3.85), ylim=(-1, 20), tick_freqs=(0.2, 5)):
    '''
    Plots stacking-fault energy (relative to energy at the minimum of AB)
    '''
    d = read.read_kc_best(align='min')

    stackings = d.stacking.unique()
    methods = d.method.unique()
    fig, axs = plt.subplots(nrows=1, ncols=len(stackings), sharey=True, figsize=(7, 3))
    for i, stacking in enumerate(stackings):
        ax = axs[i]
        for j, method in enumerate(methods):
            mask = (d.stacking == stacking) & (d.method == method)
            dd = d.loc[mask, :]
            ax.plot(dd.d, dd.energy, '-', alpha=1, lw=1.5, ms=1, label=method)

        alphabet = chr(ord('a') + i)
        ax.set(xlabel='$d~(\\mathrm{\\AA})$', xlim=xlim, title=f'({alphabet}) {stacking}')
        if i == 0:
            ax.set(ylabel='$E - E_{\\mathrm{min}}$ (meV/atom)', ylim=ylim)

    custom_ticks(axs, xlim[0], xlim[1], ylim[0], ylim[1], tick_freqs[0], tick_freqs[1])
    axs[0].legend(frameon=False)
    fig.tight_layout()
    plt.savefig('pes_methods_eq.pdf', bbox_inches='tight')

def plot_far(xlim=(2.9, 7.1), ylim=(-26, 0), tick_freqs=(1, 5)):
    '''
    Plots binding energy curve (relative to energy at infinite interlayer distance)
    '''
    d = read.read_kc_best(align='far')
    d['energy'] -= d['energy_inf']
    stackings = d.stacking.unique()
    methods = d.method.unique()
    fig, axs = plt.subplots(nrows=1, ncols=len(stackings), sharey=True, figsize=(7, 3))
    for i, stacking in enumerate(stackings):
        ax = axs[i]
        for j, method in enumerate(methods):
            mask = (d.stacking == stacking) & (d.method == method)
            dd = d.loc[mask, :]
            ax.plot(dd.d, dd.energy, '-', alpha=1, lw=1.5, ms=1, label=method)

        alphabet = chr(ord('a') + i)
        ax.set(xlabel='$d~(\\mathrm{\\AA})$', xlim=xlim, title=f'({alphabet}) {stacking}')
        if i == 0:
            ax.set(ylabel='$E - E_{\\infty}$ (meV/atom)', ylim=ylim)

    custom_ticks(axs, xlim[0], xlim[1], ylim[0], ylim[1], tick_freqs[0], tick_freqs[1])
    axs[-1].legend(frameon=False)
    fig.tight_layout()
    plt.savefig('pes_methods_far.pdf', bbox_inches='tight')

if __name__ == '__main__':
    plot_near_eq()
    plot_far()
