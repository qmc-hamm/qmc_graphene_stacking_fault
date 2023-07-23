'''
Plot potential energy surface at different values of kTs
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import read

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_label(kT):
    if kT == 'inf':
        kT_str = '$\\infty$'
    else:
        legend_format = '{kT:>6}'
        kT_str = legend_format.format(kT=int(kT))
    return kT_str

def draw_arrows(ax, xy1, xy2, color):
    arrowprops = dict(arrowstyle='<->, head_width=0.2, head_length=0.3', color=color)
    ax.annotate(text='', xy=xy1, xytext=xy2, arrowprops=arrowprops)

def get_relative_pos(ratio, vmin, vmax):
    return ratio*(vmax - vmin) + vmin

def plot(method, kTs=[4, 20], ylim=None, stripe=True):
    d = read.read_data(method)
    k = read.read_kc(f'KC-{method}', kTs)
    if ylim is None:
        ymin = k.energy.min()
        ymax = k.energy.max()
        ylim = (ymin, ymax)

    stackings = d.stacking.unique()
    ncols = len(stackings)
    fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(6, 2.5), sharey=True)
    for i, stacking in enumerate(stackings):
        ax = axs[i]
        dd = d.loc[d.stacking == stacking, :]
        print(dd)
        ax.errorbar(data=dd, x='d', y='energy', yerr='energy_err', fmt='o', alpha=1, ms=5, mew=1, mec='white', color='k', label=None, zorder=30)
        for kT in kTs:
            kk = k.loc[(k.stacking == stacking) & (k.kT == kT), :]
            ax.plot(kk.d, kk.energy, '-', label=get_label(kT), lw=1.2)
            ax.set_ylim(ylim)

        if i == 0:
            ax.set_ylabel('$E$ (eV/atom)')

        ax.set_xlabel('$d~(\\mathrm{\\AA})$')
        ax.set_xlim((2.85, 7.15))
        alphabet = chr(ord('a') + i)
        ax.set_title(f'({alphabet}) {stacking}')

        if stripe:
            ax.axvspan(3.2, 3.8, color='red', alpha=0.2)
            offset = 0.1
            y = get_relative_pos(0.95, ylim[0], ylim[1])

    axs[-1].legend(fancybox=False, edgecolor='k', title='$k_{\\mathrm{B}} T$ (meV)', bbox_to_anchor=(0.4, 0.45), framealpha=1)
    fig.tight_layout()
    plt.savefig(f'pes_{method}.pdf', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    plot('QMC', [2, 4, 'inf'], ylim=(-154.615, -154.590))
    plot('DFT-D2', [2, 3, 'inf'], ylim=(-154.920, -154.890))
    plot('DFT-D3', [2, 20, 'inf'], ylim=(-154.895, -154.870))
