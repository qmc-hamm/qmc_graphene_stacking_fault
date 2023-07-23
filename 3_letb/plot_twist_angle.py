import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd

import load

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

twist_angle = '0-99'
pot = 'qmc'

def get_pot_label(pot):
    pot_labels = {
        'qmc': 'KC-QMC',
        'ouyang': 'KC-Ouyang',
        'dft_d2': 'KC-DFT-D2',
        'dft_d3': 'KC-DFT-D3'
    }
    return pot_labels[pot]

def get_pot_color(pot):
    colors = sns.color_palette()
    pot_colors = {
        'qmc': colors[0],
        'ouyang': colors[1],
        'dft_d2': colors[2],
        'dft_d3': colors[3]
    }
    return pot_colors[pot]

def get_zorder(pot):
    zorder_map = {
        'qmc': 20,
        'ouyang': 10,
        'dft_d2': 0,
        'dft_d3': 0
    }
    return zorder_map[pot]

def get_marker(pot):
    markers = {
        'qmc': 'o',
        'ouyang': 'D',
        'dft_d2': '^',
        'dft_d3': 's'
    }
    return markers[pot]

def get_ylim(twist_angle):
    ylims = {
        '0-84': (-20, 30),
        '0-93': (-10, 15),
        '0-99': (-1, 6),
        '1-05': (-5, 20),
        '1-08': (-10, 25),
        '1-16': (-30, 40)
    }
    return ylims[twist_angle]

def order_legend(ax, order):
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    for new_label in order:
        for handle, label in zip(handles, labels):
            if label == new_label:
                new_handles.append(handle)
    # ax.legend(new_handles, order, borderpad=0, fancybox=False, edgecolor='w', bbox_to_anchor=(0.42, 0.34), handletextpad=-0.3, borderaxespad=0.0)
    ax.legend(new_handles, order, frameon=True, fancybox=False, edgecolor='k', fontsize=10.5, bbox_to_anchor=(0.33, 0.9), framealpha=1.0)
    # ax.legend(new_handles, order, frameon=False, fancybox=False, edgecolor='k', fontsize=10, bbox_to_anchor=(0.5, 0.5), framealpha=1.0)



def get_gaps(twist_angle, pot):
    evals = load.load_evals(twist_angle, pot, max_energy=200)
    k_dist, k_node = load.get_k(twist_angle, pot)
    evals_elec = evals.copy()
    evals_flat = evals.copy()
    evals_hole = evals.copy()

    ylim = get_ylim(twist_angle)
    evals_elec[evals_elec <= ylim[1]] = np.nan
    evals_flat[evals_flat > ylim[1]] = np.nan
    evals_flat[evals_flat < ylim[0]] = np.nan
    evals_hole[evals_hole >= ylim[0]] = np.nan

    elec_min = np.nanmin(evals_elec)
    flat_max = np.nanmax(evals_flat)
    flat_min = np.nanmin(evals_flat)
    hole_max = np.nanmax(evals_hole)

    elec_gap = elec_min - flat_max
    flat_gap = flat_max - flat_min
    hole_gap = flat_min - hole_max
    return elec_gap, flat_gap, hole_gap

def gen_gaps():
    twist_angles = ['0-84', '0-93', '0-99', '1-05', '1-08', '1-16']
    pots = ['qmc', 'ouyang', 'dft_d2', 'dft_d3']
    l = []
    for twist_angle in twist_angles:
        for pot in pots:
            elec_gap, flat_gap, hole_gap = get_gaps(twist_angle, pot)
            l.append({
                'twist_angle': float(twist_angle.replace('-', '.')),
                'pot': pot,
                'elec_gap': elec_gap,
                'flat_gap': flat_gap,
                'hole_gap': hole_gap
                })
    d = pd.DataFrame(l)
    d.to_csv('processed/gaps.csv', index=False)

def plot_gaps():
    d = pd.read_csv('processed/gaps.csv')
    print(d)
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3, 6), sharex=True)
    for pot in d.pot.unique():
        dd = d.loc[d.pot == pot, :]
        zorder = get_zorder(pot)
        label = get_pot_label(pot)
        color = get_pot_color(pot)
        marker = get_marker(pot)
        style = {'zorder': zorder, 'label': label, 'color': color, 'ms': 5, 'mew': 1, 'mec': 'white'}

        ax = axs[0]
        ax.plot(dd['twist_angle'], dd['flat_gap'], f'{marker}-', **style)
        ax.set_ylabel('Flat bandwidth (meV)')
        ax.set_ylim(0, 60)
        ax.set_title(f'(a)', x=0.03, y=0.87, horizontalalignment='left', transform=ax.transAxes, fontsize=12)

        ax = axs[1]
        ax.plot(dd['twist_angle'], dd['elec_gap'], f'{marker}-', **style)
        ax.set_ylabel('Electron gap (meV)')
        ax.set_ylim(0, 60)
        ax.set_title(f'(b)', x=0.03, y=0.87, horizontalalignment='left', transform=ax.transAxes, fontsize=12)

        ax = axs[2]
        ax.plot(dd['twist_angle'], dd['hole_gap'], f'{marker}-', **style)
        ax.set_ylabel('Hole gap (meV)')
        ax.set_ylim(0, 60)
        ax.set_title(f'(c)', x=0.03, y=0.87, horizontalalignment='left', transform=ax.transAxes, fontsize=12)
    ax.set_xlabel('Twist angle (degrees)')
    ax.set_xlim(0.8, 1.2)
    # axs[0].legend(frameon=False, fancybox=False, edgecolor='k', fontsize=9, bbox_to_anchor=(0.4, 0.45), framealpha=1.0)

    order_legend(axs[2], ['KC-DFT-D2', 'KC-QMC', 'KC-Ouyang', 'KC-DFT-D3'])

    plt.savefig('gaps.pdf', bbox_inches='tight')

if __name__ == '__main__':
    # gen_gaps()
    plot_gaps()
