'''
Plot the coorugation in 3D

Author: Kittithat Krongchon
'''
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource
from matplotlib import cbook
from matplotlib import colors
from matplotlib import cm
import pandas as pd
import numpy as np
import scipy.interpolate
import seaborn as sns
import re

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def rotate_coords(coords, theta):
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
        ])
    rotated = rotation.dot(coords.T).T
    return rotated

def dumpfile_to_df(dumpfile):
    '''
    Reads lammps dump file to data frame
    '''
    with open(dumpfile) as f:
        text = f.read()
    m = re.search('ITEM: BOX BOUNDS pp pp pp\n0.0000000000000000e\+00\s+(\S+)\n0.0000000000000000e\+00\s+(\S+)\s+0.0000000000000000e\+00\s+(\S+)', text)
    latvec = [
        [float(m.group(1)), 0.0, 0.0],
        [0.0, float(m.group(2)), 0.0],
        [0.0, 0.0, float(m.group(3))]
        ]

    lines = [l.split() for l in text.split('\n')[9:] if l]
    df = pd.DataFrame(lines).astype(np.float64)
    df.columns = ['atom_id', 'atom_type', 'x', 'y', 'z', 'c_csym', 'stress_1', 'stress_2', 'stress_3', 'stress_4', 'stress_5', 'stress_6', 'energy']
    df['atom_id'] = df['atom_id'].astype(int)
    df['atom_type'] = df['atom_type'].astype(int)
    df = df.sort_values(by='atom_id').reset_index(0, drop=True)

    return df, latvec

def get_data(twist_angle, potential, relaxed=True, atom_type=1, rotate=False):
    '''
    Creates a data frame from `twist_angle` and `potential`.
    `atom_type` (int):
        1 for the bottom layer
        2 for the top layer
    '''
    which_dump = 'final' if relaxed else 'initial'
    d_all, latvec = dumpfile_to_df(f"kc_{potential}/{twist_angle}/dump_{which_dump}.txt")

    # select only one layer
    d = d_all.loc[d_all.atom_type == atom_type, :].reset_index(0, drop=True)

    if rotate:
        # divide the angle by two because the geometry is setup such that
        # the top layer is twisted counter-clockwise by `twist_angle`/2,
        # and the bottom layer is twisted clockwise by `twist_angle`/2
        theta = (float(twist_angle.replace('-', '.'))/2) *np.pi/180
    else:
        theta = 0

    d_rotated = d.copy()
    d_rotated[['x', 'y', 'z']] = rotate_coords(d[['x', 'y', 'z']], theta)
    d_rotated['potential'] = potential

    # tile the atoms using periodic boundary condition
    tile1 = d_rotated.copy()
    tile1['x'] -= latvec[0][0]

    tile2 = d_rotated.copy()
    tile2['y'] += latvec[1][1]

    tile3 = d_rotated.copy()
    tile3['x'] -= latvec[0][0]
    tile3['y'] += latvec[1][1]
    d = pd.concat([tile1, tile2, tile3, d_rotated], ignore_index=True)
    return d

def get_label(pot):
    label_map = {
        'qmc': 'KC-QMC',
        'ouyang': 'KC-Ouyang',
        'dft_d2': 'KC-DFT-D2',
        'dft_d3': 'KC-DFT-D3'
    }
    return label_map[pot]

def set_axes_equal(ax: plt.Axes):
    '''Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale.
    Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
    '''
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        # ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.3 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    # ax.set_zlim3d([z - radius, z + radius])

def get_spline(d, n=200, fix_z=False):
    x = d.x
    y = d.y
    z = d.z

    x_grid = np.linspace(x.min(), x.max(), n)
    y_grid = np.linspace(y.min(), y.max(), n)
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')

    if fix_z:
        Z = np.ones((n, n))*z[0]
    else:
        spline = scipy.interpolate.Rbf(x, y, z, function='thin_plate', smooth=5, episilon=5)
        Z = spline(X, Y)

    return X, Y, Z

def plot_surface(ax, d, alpha=1, fix_z=False):
    X, Y, Z = get_spline(d, fix_z=fix_z)
    ls = LightSource(270, 45)
    rgb = ls.shade(Z, cmap=cm.coolwarm, vert_exag=0.1, blend_mode='soft')
    ax.plot_surface(X, Y, Z, color=cm.coolwarm(0))


def get_divnorm(d, vcenter=None):
    vmin = d.z.min()
    vmax = d.z.max()
    if vcenter is None:
        vcenter = d.z.mean()
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    return divnorm

def label_stackings(ax, d):

    def draw(x, y, z, s, color, offset=0, zmax=4.5, fontsize=11):
        Cx = d.x.mean()
        Cy = d.y.mean()
        X = [Cx + x + offset, Cx + x + offset]
        Y = [Cy + y + offset, Cy + y + offset]
        Z = [z+0.1, zmax]
        ax.plot(X, Y, Z, '-', lw=0.8, zorder=40, color=color)
        ax.text(Cx + x - 35, Cy + y - 10, zmax + 0.2, s, zorder=40, fontsize=fontsize, color=color)

    palette = sns.color_palette()
    Lx = d.x.max() - d.x.min()
    Ly = d.y.max() - d.y.min()
    # draw(0, 0, d.z.max(), 'AA', palette[3])
    draw(-Lx/4, -Ly/4, d.z.max(), 'AA', palette[3])
    # draw(Lx/4, -Ly/4, d.z.max(), 'AA', palette[3])
    draw(0, -Ly/6, (d.z.min() + d.z.max())/2 + 0.1, 'AB', palette[0])

def rescale(v, vmin, vmax):
    # return (v - v.min())/(v.max() - v.min())*(vmax - vmin) + vmin
    return  (v - v.min())/(vmax - vmin)*(v.max() - v.min()) + v.min()

def plot_tbg(ax, twist_angle, pot, label, relax=True, label_on=False, vmin=None, vmax=None):
    b = get_data(twist_angle, pot, relaxed=relax, atom_type=1)
    t = get_data(twist_angle, pot, relaxed=relax, atom_type=2)
    center = (b.z.mean() + t.z.mean())/2
    b.z -= center
    t.z -= center
    mid_t = (t.z.max() + t.z.min())/2
    b.z -= mid_t
    t.z -= mid_t
    mid_b = (b.z.max() + b.z.min())/2

    if relax:
        print(pot)
        print(t.z.max())
        print(t.z.min())
        ax.scatter(b.x, b.y, b.z, marker='o', s=0.5, c=b.z, cmap='coolwarm', zorder=10, vmin=vmin + mid_b, vmax=vmax + mid_b)
        p = ax.scatter(t.x, t.y, t.z, marker='o', s=0.5, c=t.z, cmap='coolwarm', zorder=20, vmin=vmin, vmax=vmax)

        if label_on:
            label_stackings(ax, t)
    else:
        plot_surface(ax, t, fix_z=True)
        plot_surface(ax, b, fix_z=True)

    ax.set_xlabel('$x~(\\mathrm{\\AA})$')
    ax.set_ylabel('$y~(\\mathrm{\\AA})$')
    ax.set_zlim(-4, 4)
    ax.set_zlabel('$z~(\\mathrm{\\AA})$')
    ax.set_box_aspect([1, 1, 1])
    set_axes_equal(ax)
    ax.axis('off')

    ax.set_title(label, y=0.98, transform=ax.transAxes, fontsize=12)
    ax.view_init(30, 300 + 5)
    return p

def create_colorbar(fig, ax, zlabel):
    '''
    Makes color bar for the scatterplotter and quiverplotter
    '''
    pos = [
        ax.get_position().x1 + 0.1, # x-coordinate corner
        ax.get_position().y0 + -0.05, # y-coordinate corner
        0.01, # width
        ax.get_position().height*1.2 # height
        ]
    cax = fig.add_axes(pos)
    cax.set_title(zlabel, x=1.15, y=1.015, transform=ax.transAxes)
    return cax

def get_clean_ticks(vmin, vmax, freq):
    tick_min = np.ceil(vmin/freq)*freq
    tick_max = np.floor(vmax/freq)*freq
    ticks = np.arange(tick_min, tick_max+freq, freq)
    return ticks

def plot_tbg_wrapper(twist_angle, pots, ext='png', vmin=-0.06, vmax=0.06, with_cbar=True):
    ncols = len(pots)
    fig, axs = plt.subplots(ncols=ncols, subplot_kw=dict(projection='3d'), figsize=(7, 3))

    for i, pot in enumerate(pots):
        if ncols == 1:
            ax = axs
        else:
            ax = axs[i]
        letter = chr(ord('a') + i)
        label = f'({letter}) {get_label(pot)}'
        label_on = True if pot == 'qmc' else False
        p = plot_tbg(ax, twist_angle, pot, label, relax=True, vmin=vmin, vmax=vmax, label_on=label_on)

    if with_cbar:
        cbar = ax.figure.colorbar(p, cax=create_colorbar(fig, ax, '$\\delta z~(\\mathrm{\\AA})$'))
        cbar.set_ticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])

    fig.tight_layout()
    cbar_label = '_cbar' if with_cbar else ''
    plt.savefig(f'{twist_angle}_3d{cbar_label}.{ext}', dpi=600, bbox_inches='tight', transparent=True)

if __name__ == '__main__':
    twist_angle = '0-99'
    plot_tbg_wrapper(twist_angle, ['qmc', 'ouyang', 'dft_d2', 'dft_d3'])
