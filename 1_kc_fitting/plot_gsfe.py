'''
Plot stacking-fault energy as a function of disregistry from different methods
'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize
import seaborn as sns

import read_min

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def get_params_energy():
    '''
    Gets fitting parameters for the stacking-fault energy curve from
    S Zhou, J Han, S Dai, J Sun, DJ Srolovit, Physical Review B (2015)
    https://doi.org/10.1103/PhysRevB.92.155438
    '''
    conv = read_min.zhou_to_ev_per_atom()
    c0 = 21.336*conv
    c1 = -6.127*conv
    c2 = -1.128*conv
    c3 = 0.143*conv
    return c0, c1, c2, c3

def get_params_spacing():
    '''
    Gets fitting parameters for the relaxed interlayer spacing curve from
    S Zhou, J Han, S Dai, J Sun, DJ Srolovit, Physical Review B (2015)
    https://doi.org/10.1103/PhysRevB.92.155438
    '''
    c0 = 3.47889
    c1 = -0.02648
    c2 = -0.00352
    c3 = 0.00037
    return c0, c1, c2, c3

def get_color(method):
    palette = sns.color_palette()
    color_map = {
        'KC-QMC': palette[0],
        'KC-Ouyang': palette[1],
        'KC-DFT-D2': palette[2],
        'KC-DFT-D3': palette[3],
        'RPA': palette[4]
        }
    return color_map[method]

def get_zorder(method):
    zorder_map = {
        'KC-QMC': 25,
        'KC-Ouyang': 20,
        'RPA': 15,
        'KC-DFT-D2': 10,
        'KC-DFT-D3': 5
        }
    return zorder_map[method]

def get_clean_ticks(vmin, vmax, freq):
    '''
    customize axis ticks so that the ticks are at the given frequency
    e.g. vmin = 3.15, vmax = 3.85, freq = 0.2 will result in ticks at 3.2, 3.4, 3.6, 3.8
    '''
    tick_min = np.ceil(vmin/freq)*freq
    tick_max = np.floor(vmax/freq)*freq
    ticks = np.arange(tick_min, tick_max+freq, freq)
    return ticks

def label_stackings(axs):
    '''
    Annotates the stacking types on the given axes
    '''
    style = {'color': '#2f2f2f', 'ls': 'dotted', 'alpha': 0.6, 'zorder': 0, 'lw': 1}
    for curr_ax in axs:
        curr_ax.axvline(x=0, **style)
        curr_ax.axvline(x=1/6, **style)
        curr_ax.axvline(x=0.5, **style)
        curr_ax.axvline(x=2/3, **style)
        curr_ax.axvline(x=1/3, **style)
        curr_ax.axvline(x=1, **style)

    dy = 0.4
    y = 8 + dy
    dx = -0.04
    fontsize = 9
    axs[0].text(0+dx, y, 'AB', fontsize=fontsize)
    axs[0].text(1/6+dx, y, 'SP', fontsize=fontsize)
    axs[0].text(0.5+dx, y, 'Mid', fontsize=fontsize)
    axs[0].text(2/3+dx, y, 'AA', fontsize=fontsize)
    axs[0].text(1/3+dx, y, 'BA', fontsize=fontsize)
    axs[0].text(1+dx, y, 'AB', fontsize=fontsize)
    axs[1].set(xticks=[0, 1/6, 1/2, 2/3, 1], xticklabels=['0', '1/6', '1/2', '2/3', '1'])

def order_legend(ax, order):
    '''
    Custom sort legend items
    '''
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    for new_label in order:
        for handle, label in zip(handles, labels):
            if label == new_label:
                new_handles.append(handle)
    ax.legend(new_handles, order, borderpad=0, fancybox=False, edgecolor='w', bbox_to_anchor=(0.42, 0.34), handletextpad=-0.3, borderaxespad=0.0)

def print_table(consts, order):
    order_map = {order[i]:i for i in range(len(order))}
    d = pd.DataFrame(consts)
    d['method_order'] = d.method.map(order_map)
    d = d.sort_values(by='method_order')
    d = d[['method', 'c0', 'c1', 'c2', 'c3']]
    s = d.round(5).to_latex(index=False)
    print(s)

def F(y, c0, c1, c2, c3):
    '''
    Functioal forms of the stacking-fault energy and the relaxed interlayer spacing curves from
    S Zhou, J Han, S Dai, J Sun, DJ Srolovit, Physical Review B (2015)
    https://doi.org/10.1103/PhysRevB.92.155438
    '''
    a0 = 2.46 # ang
    A = 2*np.pi/a0

    x = 0 # registry in the zigzag direction
    y = 3**0.5*a0 - y # registry in the armchair direction

    c4 = 3**0.5*c1
    c5 = -3**0.5*c3

    F1 = np.cos(A*(x + 1/3**0.5*y)) + np.cos(A*(x - 1/3**0.5*y)) + np.cos(A*2/3**0.5*y)
    F2 = np.cos(A*(x + 3**0.5*y)) + np.cos(A*(x - 3**0.5*y)) + np.cos(A*2*x)
    F3 = np.cos(A*(2*x + 2/3**0.5*y)) + np.cos(A*(2*x - 2/3**0.5*y)) + np.cos(A*4/3**0.5*y)
    F4 = np.sin(A*(x - 1/3**0.5*y)) - np.sin(A*(x + 1/3**0.5*y)) + np.sin(A*2/3**0.5*y)
    F5 = np.sin(A*(2*x - 2/3**0.5*y)) - np.sin(A*(2*x + 2/3**0.5*y)) + np.sin(A*4/3**0.5*y)
    return c0 + c1*F1 + c2*F2 + c3*F3 + c4*F4 + c5*F5

def plot(methods, output='gsfe.pdf'):
    e = read_min.get_data('energy')
    r = read_min.get_data('d_eq')
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.5, 5), sharex=True)
    consts_energy = []
    consts_rmin = []
    for method in e.method.unique():
        print(method)
        if method not in methods:
            continue

        ee = e.loc[e.method == method, :]
        rr = r.loc[r.method == method, :]

        a0 = 2.46
        lin = np.linspace(-0.5, 1.5, 500)
        color = get_color(method)
        zorder = get_zorder(method)

        subplot_title_x = 0.1
        subplot_title_y = 0.83

        popt, pcov = scipy.optimize.curve_fit(F, ee.registry_ang, ee.energy, p0=get_params_energy())
        consts_energy.append({'method': method, 'c0': popt[0], 'c1': popt[1], 'c2': popt[2], 'c3': popt[3]})
        axs[0].plot(lin, F(lin*3**0.5*a0, *popt), '-', lw=1, color=color, zorder=zorder)
        axs[0].plot(ee.registry, ee.energy, 'o', ms=5, mew=1, mec='white', color=color, zorder=zorder, label=method)
        ymin = -0.5
        ymax = 8
        yticks = get_clean_ticks(ymin, ymax, 2)
        axs[0].set(ylabel='SFE (meV/atom)', ylim=(ymin, ymax), yticks=yticks)
        axs[0].set_title('(a)', x=subplot_title_x, y=subplot_title_y)

        popt, pcov = scipy.optimize.curve_fit(F, rr.registry_ang, rr.d_eq, p0=get_params_spacing())
        consts_rmin.append({'method': method, 'c0': popt[0], 'c1': popt[1], 'c2': popt[2], 'c3': popt[3]})
        axs[1].plot(lin, F(lin*3**0.5*a0, *popt), '-', lw=1, color=color, zorder=zorder)
        axs[1].plot(rr.registry, rr.d_eq, 'o', ms=5, mew=1, mec='white', color=color, zorder=zorder)
        axs[1].set(ylabel='$d_{\\mathrm{min}}$ ($\\mathrm{\\AA}$)', ylim=(3.2, 3.8))
        axs[1].set_title('(b)', x=subplot_title_x, y=subplot_title_y)

        axs[1].set(xlabel='Registry ($\\sqrt{{3}} a$)', xlim=(-0.05, 1.05))

    order = methods

    print_table(consts_energy, order)
    print_table(consts_rmin, order)

    label_stackings(axs)
    order_legend(axs[0], order)

    fig.tight_layout()
    plt.savefig(output, dpi=400, bbox_inches='tight')


if __name__ == '__main__':
    plot(['KC-DFT-D2', 'KC-Ouyang', 'KC-QMC', 'RPA', 'KC-DFT-D3'], output='gsfe.pdf')
