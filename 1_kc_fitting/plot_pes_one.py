'''
Plot potential energy surface within one subplot for the workflow figure.
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import read

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_color(stacking):
    palette = sns.color_palette()
    color_map = {
        'AB': palette[0],
        'SP': palette[2],
        'Mid': palette[1],
        'AA': palette[3]
    }
    return color_map[stacking]

def order_legend(ax, fontsize=14):
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [handles[get_legend_order(label)] for label in labels]
    new_labels = [labels[get_legend_order(label)] for label in labels]
    print(new_labels)
    ax.legend(new_handles, new_labels, frameon=False, fontsize=fontsize, bbox_to_anchor=(0.65, 0.5), handletextpad=-0.3, borderaxespad=0.0, borderpad=0)

def get_legend_order(label):
    order_map = {
        'AB': 3,
        'SP': 2,
        'Mid': 1,
        'AA': 0
    }
    return order_map[label]

def get_clean_ticks(vmin, vmax, freq):
    '''
    customize axis ticks so that the ticks are at the given frequency
    e.g. vmin = 3.15, vmax = 3.85, freq = 0.2 will result in ticks at 3.2, 3.4, 3.6, 3.8
    '''
    tick_min = np.ceil(vmin/freq)*freq
    tick_max = np.floor(vmax/freq)*freq
    ticks = np.arange(tick_min, tick_max+freq, freq)
    return ticks

def plot(method, kT, align, fmt='o-', ext='pdf'):
    d = read.read_data(method)
    k = read.read_kc(f'KC-{method}', [kT])
    d.energy *= 1000
    d.energy_err *= 1000
    k.energy *= 1000
    if align == 'inf':
        energy_inf = k.energy_inf.values[0]*1000
        d.energy -= energy_inf
        k.energy -= energy_inf
    elif align == 'min':
        energy_min = k.energy.min()
        d.energy -= energy_min
        k.energy -= energy_min

    d.to_csv(f'processed/data_{method}_{kT}_{align}.csv', index=False)
    k.to_csv(f'processed/kc_{method}_{kT}_{align}.csv', index=False)

    fig, ax = plt.subplots(figsize=(3, 3))
    for stacking in d.stacking.unique():
        dd = d.loc[d.stacking == stacking, :]
        kk = k.loc[(k.stacking == stacking) & (k.kT == kT), :]
        label = f'{stacking}'
        if '-' in fmt:
            ax.plot(kk.d, kk.energy, '-', lw=1, color=get_color(stacking), label=None)
        if 'o' in fmt:
            ax.errorbar(data=dd, x='d', y='energy', yerr='energy_err', fmt='o', alpha=1, ms=5, mew=1, mec='white', color=get_color(stacking), label=label)

    order_legend(ax, fontsize=15)
    ax.set_xlabel('$d~(\\mathrm{\\AA})$', fontsize=15)
    ax.set_ylabel('$E~(\\mathrm{eV/atom})$', fontsize=15)
    if kT == 4:
        xmin, xmax, xfreq = (2.85, 5.15, 1)
        ymin, ymax, yfreq = (-1, 20, 5)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        xticks = get_clean_ticks(xmin, xmax, xfreq)
        yticks = get_clean_ticks(ymin, ymax, yfreq)
        ax.set(xticks=xticks, yticks=yticks)
    elif kT == 'inf':
        ax.set_ylim(-25, 0)
    plt.savefig(f'pes_one_{method}_{kT}_{fmt}.{ext}', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    plot('QMC', 4, 'min', ext='png', fmt='-o')
