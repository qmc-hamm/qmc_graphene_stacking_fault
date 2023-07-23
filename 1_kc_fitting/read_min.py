'''
Prepare stacking-fault energy data for plotting
'''
import copy
import pandas as pd
import numpy as np

def zhou_to_ev_per_atom(a=2.46e-10):
    '''
    Converts the units of mJ/m^2 to eV/atom for bilayer graphene
    '''
    A = 3**0.5/2*a**2 # m^2
    # mJ/m^2 * A m^2 / unitcell * 1 unitcell / 4 atoms = 0.001*convert(1.0, 'J', 'eV')*A/4
    j_to_ev = 1/1.602176634e-19
    return j_to_ev*0.001*A/4

def read_zhou(which):
    '''
    Reads the stacking-fault energy (`energy`) or relaxed interlayer spacing (`d_eq`) as digitized from
    S Zhou, J Han, S Dai, J Sun, DJ Srolovit, Physical Review B (2015)
    https://doi.org/10.1103/PhysRevB.92.155438
    '''
    if which == 'energy':
        d = pd.read_csv('refs/zhou_rpa_energy.csv')
        d.energy *= zhou_to_ev_per_atom()*1000 # convert to meV
        yvar = 'energy'
    elif which == 'd_eq':
        d = pd.read_csv('refs/zhou_rpa_spacing.csv')
        d['d_eq'] = d['spacing']
        yvar = 'd_eq'
    d['method'] = 'RPA'
    d['registry'] = d['disregistry']
    d = d[['method', 'registry', yvar]]
    return d

def dupe_registries(d):
    '''
    Duplicates the energy data to the equivalent registry (due to symmetry),
    e.g. the energy at s=0 (AB) should be the same as s=1/3 (BA),
    and the energy at s=0.5 (Mid) should be the same as s=5/6.
    '''
    l = []
    for row in d.to_dict(orient='records'):
        stacking = row['stacking']
        if stacking == 'AB':

            # append original
            l.append(row)

            # append BA
            row2 = copy.deepcopy(row)
            row2['stacking'] = 'BA'
            row2['registry'] = 1/3
            l.append(row2)

            # append AB
            row3 = copy.deepcopy(row)
            row3['registry'] = 1
            l.append(row3)

        elif stacking == 'Mid':

            # append original
            l.append(row)

            # append 5/6
            row2 = copy.deepcopy(row)
            row2['registry'] = 5/6
            l.append(row2)
        else:
            l.append(row)

    d2 = pd.DataFrame(l)
    return d2

def get_data(which, lim='min', a0=2.46):
    '''
    Returns the stacking fault energy data including refs
    '''
    if lim == 'min':
        csv_name = 'processed/min.csv'
    elif lim == 'far':
        csv_name = 'processed/far.csv'
    d = pd.read_csv(csv_name)
    d['energy'] = d['gsfe']
    d['registry'] = d['stacking'].map({'AB': 0.0, 'SP': 1/6, 'Mid': 0.5, 'AA': 2/3})
    d = d.groupby('method').apply(dupe_registries).reset_index(0, drop=True)
    if which == 'energy':
        d = d[['method', 'registry', 'energy']]
    elif which == 'd_eq':
        d = d[['method', 'registry', 'd_eq']]
    elif which == 'be':
        d = d[['method', 'registry', 'be']]

    if which != 'be':
        rpa = read_zhou(which)
        d = pd.concat([d, rpa], ignore_index=True)
    d['registry_ang'] = d['registry']*3**0.5*a0
    d['stacking'] = d['registry'].map({0: 'AB', 1/6: 'SP', 1/2: 'Mid', 2/3: 'AA'})
    return d

def test():
    print(get_data('energy'))
    print(get_data('d_eq'))
    print(get_data('be'))

def print_table(d, value):
    d = d.drop(['registry', 'registry_ang'], axis=1)
    d = d.dropna()
    d = d.pivot(index='method', columns='stacking', values=value)
    d = d[['AB', 'SP', 'Mid', 'AA']].reset_index(0)
    d = d.round(3).to_latex(index=False)
    print(d)

if __name__ == '__main__':
    print_table(get_data('d_eq', lim='min'), 'd_eq')
    print_table(get_data('be', lim='far'), 'be')
