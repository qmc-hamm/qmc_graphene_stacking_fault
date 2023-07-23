'''
Plot RMS as a function of kT (fitting weight) parameters to choose which kT is the best for a method
'''
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_color(lim):
    return {
        'min': 'red',
        'far': 'black'
        }[lim]

def get_label(lim):
    return {
        'min': '$d \\in [3.2, 3.8]~\\mathrm{\\AA}$',
        'far': '$d \\in [3.0, 7.0]~\\mathrm{\\AA}$'
        }[lim]

def plot_rms(method, ylim=None, x_eq=None):
    d = pd.read_csv('processed/stats.csv')
    d['d_range'] = d['lim'].map(get_label)
    d = d.loc[d.method == method, :]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2, 2))
    style = {'lw': 1, 'mec': 'white', 'mew': 1, 'ms': 5}

    for lim in d.lim.unique():
        if lim != 'min':
            continue
        dd = d.loc[d.lim == lim, :]
        min_idx = dd.rms.argmin()
        print(method, lim, dd.iloc[min_idx, :]['kT'])
        label = get_label(lim)
        color = get_color(lim)
        ax.plot(dd.kT, dd.rms*1000, marker='o', color=color, label=label, **style)
        if x_eq is not None:
            ax.axvline(x=x_eq, ymin=0, ymax=0.28, ls=':', color=color)
    ax.set(xlabel='$k_{\\mathrm{B}} T$ (meV)')
    ax.set(ylabel='RMS (meV)', ylim=ylim)
    ax.set(title='(e)')
    plt.savefig(f'rms_{method}.pdf', bbox_inches='tight')

if __name__ == '__main__':
    plot_rms('QMC', ylim=(0.26, 0.36), x_eq=4)
    plot_rms('DFT-D2', x_eq=3)
    plot_rms('DFT-D3', x_eq=20)
