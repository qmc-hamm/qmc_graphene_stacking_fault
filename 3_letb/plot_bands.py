import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd

import load

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_pot_label(pot):
    pot_labels = {
        'qmc_rigid': 'KC-QMC Rigid',
        'qmc': 'KC-QMC',
        'ouyang': 'KC-Ouyang',
        'dft_d2': 'KC-DFT-D2',
        'dft_d3': 'KC-DFT-D3'
    }
    return pot_labels[pot]

def get_pot_color(pot):
    colors = sns.color_palette()
    pot_colors = {
        'qmc_rigid': colors[0],
        'qmc': colors[0],
        'ouyang': colors[1],
        'dft_d2': colors[2],
        'dft_d3': colors[3]
    }
    return pot_colors[pot]

def get_pot_line_style(pot):
    pot_line_styles = {
        'qmc_rigid': '-',
        'qmc': '-',
        'ouyang': '-',
        'dft_d2': '-',
        'dft_d3': '-'
    }
    return pot_line_styles[pot]

def get_zorder(pot):
    zorder_map = {
        'qmc_rigid': 20,
        'qmc': 20,
        'ouyang': 10,
        'dft_d2': 0,
        'dft_d3': 0
    }
    return zorder_map[pot]

def get_ylim(twist_angle):
    ylims = {
        '0-84': (-25, 30),
        '0-93': (-10, 15),
        '0-99': (-1, 6),
        '1-05': (-5, 20),
        '1-08': (-10, 25),
        '1-16': (-30, 40)
    }
    return ylims[twist_angle]

def plot_bands_zoom(twist_angle, pots, label=True, ext='png'):
    fig, axs = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True, figsize=(6.5, 1.5))
    for j, pot in enumerate(pots):
        ax = axs[j]
        k_dist, k_node = load.get_k(twist_angle, pot)
        evals = load.load_evals(twist_angle, pot, max_energy=200)

        palette = sns.color_palette()
        for i in range(evals.shape[0]):
            ax.plot(k_dist, evals[i, :], get_pot_line_style(pot), color=palette[0], label=None)

        if j == 0:
            ax.set_ylabel('$E - E_{{\\mathrm{F}}}$ (meV)')
            ax.set_ylim(get_ylim(twist_angle))
        letter = chr(ord('a') + j + 4)
        if label:
            ax.set_title(f'({letter})', x=0.09, y=0.84, horizontalalignment='left', transform=ax.transAxes, fontsize=10)

        ax.set_xticks(k_node)
        ax.set_xticklabels(["K", "$\\Gamma$", "M", "$\\mathrm{K}^{\\prime}$"])
        ax.set_xlim(k_node[0], k_node[-1])

    fig.savefig(f'{twist_angle}_bands_zoom.{ext}', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    pots = ['qmc', 'ouyang', 'dft_d2', 'dft_d3']
    # plot_bands_zoom('0-84', pots, ext='pdf', label=False)
    # plot_bands_zoom('0-93', pots, ext='pdf', label=False)
    plot_bands_zoom('0-99', pots, ext='pdf', label=True)
    # plot_bands_zoom('1-05', pots, ext='pdf', label=False)
    # plot_bands_zoom('1-08', pots, ext='pdf', label=False)
    # plot_bands_zoom('1-16', pots, ext='pdf', label=False)
