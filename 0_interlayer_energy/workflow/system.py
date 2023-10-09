'''
Define geometry of bilayer graphene.
'''
from nexus import generate_physical_system
import numpy as np

def get_basis(a, d, c, disregistry, zshift='CM'):

    '''
    `disregistry` is unitless, defined such that the distance to translate from AB to AB again is 1.0.
    In angstrom, this distance corresponds to 3*bond_length = 3/sqrt(3)*lattice_constant = sqrt(3)*lattice_constant
    See Fig. 1 in [paper]
    '''
    disregistry_ang = 3**0.5*a*disregistry
    orig_basis = np.array([
        [0, 0, 0],
        [0, a/3**0.5, 0],
        [0, a/3**0.5 + disregistry_ang, d],
        [a/2, a/(2*3**0.5) + disregistry_ang, d]
        ])

    # for open boundary condition in the z-direction
    # move the first layer to the middle of the cell
    if zshift == 'first_layer':
        z = c/2
    # or move the center of mass to the middle of the cell
    elif zshift == 'CM':
        z = c/2 - d/2
    shift_vector = np.array([0, 0, z])
    shifted_basis = orig_basis + shift_vector
    return shifted_basis.tolist()

def get_lattice_vectors(a, c):
    '''
    Defines lattice vectors for bilayer graphene
    '''
    return [
        [a, 0, 0],
        [1/2*a, 1/2*3**0.5*a, 0],
        [0, 0, c]
        ]

def create_graphene_geom(a, d, disregistry, qmc_kgrid, tiling, c=20, kshift=0.0, zshift='CM'):
    '''
    Creates a Nexus geometry object for the Nexus-monitored SCF and Monte Carlo workflow
    '''
    return generate_physical_system(
        units = 'A',
        axes = get_lattice_vectors(a, c),
        elem = ['C']*4,
        pos = get_basis(a, d, c, disregistry, zshift=zshift),
        kgrid = qmc_kgrid,
        kshift = [kshift, kshift, 0.0],
        net_charge = 0,
        net_spin = 0,
        tiling = tiling,
        C = 4
        )
