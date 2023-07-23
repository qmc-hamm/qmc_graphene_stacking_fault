from concurrent.futures import ProcessPoolExecutor
import json
import os
import sys

import ase
from pythtb import *
from bilayer_letb.api import tb_model

def read_poscar(poscar_path):
    with open(poscar_path) as f:
        text = f.read()
    lines = text.split('\n')

    latvec = [[float(coord) for coord in l.split()] for l in lines[2:5]]
    atoms = [[float(coord) for coord in l.split()] for l in lines[8:-1]]
    return latvec, atoms

def create_letb(twist_angle, pot, relax_keyword):
    '''
    `relax_keyword` is either `relax` or `rigid`
    '''
    poscar_path = f'../2_optimized_geometry/kc_{pot}/raw/POSCAR_{twist_angle}{relax_keyword}_hex.txt'

    if os.path.isfile(poscar_path):
        # compute hoppings
        lattice_vectors, atomic_basis = read_poscar(poscar_path)
        natoms = len(atomic_basis)
        ase_atoms = ase.Atoms(['C']*natoms, positions=atomic_basis, cell=lattice_vectors, pbc=True)
        letb = tb_model(ase_atoms)
        return letb
    else:
        print(f'poscar does not exist!!!: {poscar_path}')

def solve_bands(letb, dirname, k_vec, k_idx):
    output = f'bands_{k_idx:02}.npy'

    workdir = os.getcwd()
    os.chdir(dirname)

    evals = letb.solve_one(k_vec[k_idx]).reshape(-1, 1)
    np.save(output, arr=evals)
    os.chdir(workdir)

if __name__ == '__main__':
    if len(sys.argv) == 4:
        twist_angle = sys.argv[1]
        pot = sys.argv[2]
        relax_keyword = sys.argv[3]

        with ProcessPoolExecutor(max_workers=100) as client:
            nk = 100
            letb = create_letb(twist_angle, pot, relax_keyword)
            k = [[1/3., 2/3.], [0.0, 0.0], [0.5, 0.0], [2/3., 1/3.]]
            k_vec, k_dist, k_node = letb.k_path(k, nk)
            dirname = f'bands/{twist_angle}_{pot}_{relax_keyword}'
            os.makedirs(dirname, exist_ok=True)
            np.savetxt(f'{dirname}/k_vec.txt', k_vec)
            np.savetxt(f'{dirname}/k_dist.txt', k_dist)
            np.savetxt(f'{dirname}/k_node.txt', k_node)
            jobs = [client.submit(solve_bands, letb, dirname, k_vec, k_idx) for k_idx in range(nk)]
