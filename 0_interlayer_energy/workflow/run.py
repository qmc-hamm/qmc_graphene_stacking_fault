'''
Main workflow for executing DMC calculations

Author: Kittithat Krongchon and Lucas K. Wagner
'''
import os
import itertools 
import numpy as np

import nexus
import qmcpack_input

import system
import config
import check
import opt


def get_dmc_block(timestep, dmc_blocks):
    return nexus.dmc(
        blocks = dmc_blocks,
        warmupsteps = int(2.0/timestep),
        steps = int(1.0/timestep),
        timestep = timestep,
        nonlocalmoves = True,
        gpu = False
        )

def get_vmc_block(samplesperthread):
    return nexus.vmc(
        usedrift = True,
        walkers = 1,
        warmupsteps = 50,
        blocks = 5,
        steps = 1,
        substeps = 10,
        timestep = 0.8,
        samplesperthread = samplesperthread,
        gpu = False
        )

def get_jastrow_resources(ntiling, opt_samples):
    resources = {
        2: {'qmc_nodes': 16, 'minutes': 10, 'max_loop': 12},
        3: {'qmc_nodes': 32, 'minutes': 30, 'max_loop': 12},
        4: {'qmc_nodes': 64, 'minutes': 50, 'max_loop': 12},
        5: {'qmc_nodes': 128, 'minutes': 120, 'max_loop': 15},
        6: {'qmc_nodes': 512, 'minutes': 120, 'max_loop': 15}
        }
    return resources[ntiling]

def get_dmc_resources(ntiling):
    resources = {
        2: {'dmc_nodes_per_twist': 3, 'samplesperthread': 50},
        3: {'dmc_nodes_per_twist': 6, 'samplesperthread': 22},
        4: {'dmc_nodes_per_twist': 12, 'samplesperthread': 13},
        5: {'dmc_nodes_per_twist': 12, 'samplesperthread': 13},
        6: {'dmc_nodes_per_twist': 12, 'samplesperthread': 13},
        }
    return resources[ntiling]

def generate_workflow(dirname, d, disregistry,
    ntiling = 3,
    qmc_nk = 4,
    a = 2.46,
    c = 20,
    timesteps = [0.01, 0.02],
    kshift = 0.076923,
    zshift = 'CM',
    lr_handler = 'opt_breakup',
    lr_dim_cutoff = 15,
    pw_nk = 20,
    pw_nkz = 1,
    dmc_blocks = None,
    vdw_corr = 'dft-d',
    dft_only = False,
    opt_method = 'basic',
    opt_samples = 50000,
    vmc_only = False,
    correction_name = 'none',
    ):

    qmc_kgrid = [qmc_nk, qmc_nk, 1]
    tiling = [ntiling, ntiling, 1]
    geom = system.create_graphene_geom(a, d, disregistry, qmc_kgrid, tiling, c=c, kshift=kshift, zshift=zshift)

    if vdw_corr in ['ts', 'mbd']:
        pseudo = 'C.ccECP.upf'
    else:
        pseudo = 'C.tn.upf'

    pw_options = {
        'pseudos': [pseudo],
        'system': geom,
        'ecutwfc': 200,
        'vdw_corr': vdw_corr,
        'input_dft': 'pbe',
        'conv_thr': 1e-6,
        'mixing_beta': 0.3,
        'use_folded': True
        }

    scf = nexus.generate_pwscf(
        identifier = 'scf',
        path = os.path.join(dirname, 'scf'),
        job = config.pw_job(pw_nodes=1, minutes=10),
        calculation = 'scf',
        kgrid = [pw_nk, pw_nk, pw_nkz],
        nosym = False,
        noinv = False,
        wf_collect = False,
        **pw_options
        )
    if dft_only:
        return [scf]

    nscf = nexus.generate_pwscf(
        identifier = 'nscf',
        path = os.path.join(dirname, 'nscf'),
        job = config.pw_job(pw_nodes=1, minutes=10),
        calculation = 'nscf',
        nosym = True,
        noinv = True,
        wf_collect = True,
        dependencies = [(scf, 'charge_density')],
        **pw_options
        )

    p2q = nexus.generate_pw2qmcpack(
        identifier = 'p2q',
        path = os.path.join(dirname, 'nscf'),
        job = config.p2q_job(p2q_nodes=1, minutes=10),
        write_psir = False,
        dependencies = [(nscf, 'orbitals')]
        )

    qmc_options = {
        'input_type': 'basic',
        'system': geom,
        'bconds': 'ppn',
        'randomsrc': 'ion0',
        'pseudos': ['C.tn.xml'],
        'lr_handler': lr_handler,
        'lr_dim_cutoff': lr_dim_cutoff
        }

    jastrow_resources = get_jastrow_resources(ntiling, opt_samples)
    jastrow = nexus.generate_qmcpack(
        identifier = 'jastrow',
        path = os.path.join(dirname, 'jastrow'),
        job = config.qmc_job(qmc_nodes=jastrow_resources['qmc_nodes'], minutes=jastrow_resources['minutes']),
        jastrows = [
            ('J1', 'bspline', 8),
            ('J2', 'bspline', 8, 'coeff', [8*[0], 8*[0]])
            ],
        calculations = opt.strats(opt_method, opt_samples, jastrow_resources['max_loop']),
        twistnum = 0,
        dependencies = [(p2q, 'orbitals')],
        **qmc_options
        )

    if dmc_blocks is None:
        dmc_blocks = max(int(16*200/np.prod(qmc_kgrid))+1, 20)

    if correction_name == 'none':
        corrections = []
    else:
        corrections = [correction_name]
    dmc_resources = get_dmc_resources(ntiling)
    vmc_block = get_vmc_block(dmc_resources['samplesperthread'])
    dmcs = [nexus.generate_qmcpack(
        identifier = 'dmc',
        path = os.path.join(dirname, f'dmc{timestep}'),
        job = config.qmc_job(qmc_nodes=dmc_resources['dmc_nodes_per_twist']*np.prod(qmc_kgrid), minutes=120),
        twistnum = None,
        jastrows = [],
        calculations = [vmc_block] if vmc_only else [vmc_block, get_dmc_block(timestep, dmc_blocks=dmc_blocks)],
        corrections = corrections,
        estimators = [
            qmcpack_input.estimator(
                name = 'skall',
                type = 'skall',
                source = 'ion0',
                target = 'e',
                hdf5 = 'yes'
                ),
            qmcpack_input.estimator(
                name = 'LocalEnergy',
                hdf5 = 'no'
                )
           ],
        dependencies = [(p2q, 'orbitals'), (jastrow, 'jastrow')],
        **qmc_options
        )
        for timestep in timesteps]
    return [scf, p2q, jastrow]

def run_workflows(params, timesteps=[0.01, 0.02], data_path='', generate_only=False, run_check=True, force_run=False):
    '''
    Submits all calculations defined in `params`.
    '''
    config.apply_settings(generate_only=generate_only)
    workflows = []
    for p in itertools.product(*params.values()):
        options = dict(zip(params.keys(), p))
        dirname = os.path.join(data_path, '_'.join(k+str(v) for k, v in options.items()))
        if run_check:
            check.check_outputs(dirname, timesteps, force_clean=force_run)
        workflows += generate_workflow(dirname, timesteps=timesteps, **options)
    nexus.run_project(workflows)

def example():
    '''
    Submits all the Cartesian product combinations of parameters below.
    '''
    params = {
        'd': [3.0, 3.35, 5.0],
        'disregistry': [0.0, 0.5],
        'ntiling': [2],
        'qmc_nk': [2]
        }
    run_workflows(params)

if __name__ == '__main__':
    example()
