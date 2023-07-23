import numpy as np
import os

def get_k(twist_angle, pot):
    dirname = f'bands/{twist_angle}_{pot}'
    k_dist = np.loadtxt(f'{dirname}/k_dist.txt')
    k_node = np.loadtxt(f'{dirname}/k_node.txt')
    return k_dist, k_node

def trim(l, max_energy):
    '''
    Given a sorted list `l`, return the first/last index of the element that is greater/less than -max_energy/max_energy
    For filtering energy with in the range [-max_energy, max_energy]
    Ex:
        l = np.array([-12, -7, -3, -1, 0, 2, 5, 9, 15])
        max_energy = 5
        first, last = trim(l, max_energy)
        print(l[first:last]) # [-3 -1  0  2  5]
    '''
    last_idx = sum(l <= max_energy)
    first_idx = sum(l <= -max_energy)
    return first_idx, last_idx

def test_trim():
    l = np.array([-12, -7, -3, -1, 0, 2, 5, 9, 15])
    max_energy = 5
    first, last = trim(l, max_energy)
    print(l[first:last])

def concat_evals(twist_angle, pot):
    '''
    Combines all the band energy at each k point to a single matrix, where energies are rows and k-points are columns
    '''
    l = []
    dirname = f'bands/{twist_angle}_{pot}'
    k_dist = np.loadtxt(f'{dirname}/k_dist.txt')
    nk = k_dist.shape[0]
    for k_idx in range(nk):
        evals_i = np.load(f'{dirname}/bands_{k_idx:02}.npy')
        l.append(evals_i)
    evals = np.concatenate(l, axis=1)
    np.savetxt(f'{dirname}/bands.txt', evals)

def load_evals(twist_angle, pot, max_energy=None):
    dirname = f'bands/{twist_angle}_{pot}'
    bands_txt = f'{dirname}/bands.txt'
    if not os.path.isfile(bands_txt):
        concat_evals(twist_angle, pot)
    evals = np.loadtxt(bands_txt)
    evals *= 1000 # convert eV to meV

    # center the band structure at the fermi level
    evals -= min(evals[int(evals.shape[0]/2), :])

    # filter the energy outside of the range [-max_energy, max_energy]
    first_idxs = []
    last_idxs = []
    nk = evals.shape[1]
    for k_idx in range(nk):
        if max_energy is not None:
            first_idx, last_idx = trim(evals[:, k_idx], max_energy)
        else:
            first_idx = 0
            last_idx = -1
        first_idxs.append(first_idx)
        last_idxs.append(last_idx)
    evals_within_max_energy = evals[min(first_idxs):max(last_idxs), :]
    return evals_within_max_energy
