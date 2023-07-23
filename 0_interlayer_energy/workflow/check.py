'''
This utility file recovers and cleans the calculation directories
to prepare for the subsequent submissions in case of crashing during the QMCPACK Nexus workflow.

Author: Kittithat Krongchon
'''
import os

def check_pw_output(out):
    if not os.path.isfile(out):
        print(f'{out} does not exist. cleaning this dir')
        return 'failed'

    with open(out) as f:
        text = f.read()
    if ('Maximum CPU time exceeded' in text) or ('JOB DONE' not in text):
        return 'failed'
    else:
        return 'finished'

def check_jastrow_output(out):
    if not os.path.isfile(out):
        print(f'{out} does not exist. cleaning this dir')
        return 'failed'

    with open(out) as f:
        text = f.read()
    if 'Total Execution time' not in text:
        return 'failed'
    else:
        return 'finished'

def check_dmc_output(out):
    if not os.path.isfile(out):
        print(f'{out} does not exist. cleaning this dir')
        return 'failed'

    if not os.path.getsize(out):
        return 'failed'
    else:
        return 'finished'

def get_outputs(dirname, prefix):
    '''
    Preares a set of output files according to the type of caluclation (`prefix`).
    '''
    dirprefix = 'nscf' if prefix == 'p2q' else prefix
    file_prefix = 'dmc' if 'dmc' in prefix else prefix
    outputs = {
        'sim': f'sim_{file_prefix}',
        'out': f'{file_prefix}.out',
        'err': f'{file_prefix}.err',
        'old': f'{file_prefix}.old'
        }
    return {file_type: os.path.join(dirname, dirprefix, fn) for file_type, fn in outputs.items()}

def clean(dirname, prefix):
    '''
    Removes old calculation files
    '''
    print(f'{dirname} {prefix} is being cleaned.')
    outputs = get_outputs(dirname, prefix)
    items = [outputs['sim'], outputs['out'], outputs['err']]
    old = outputs['old']
    for item in items:
        if os.path.exists(item):
            if item.endswith('.out'):
                os.system(f'mv {item} {old}')
            else:
                os.system(f'rm -r {item}')

def check_output(dirname, prefix, force_clean=False):
    '''
    Inspects the directory and cleans it if needed.
    '''
    dirprefix = 'nscf' if prefix == 'p2q' else prefix
    if not os.path.isdir(os.path.join(dirname, dirprefix)):
        print(f'{dirname} {prefix} dir not found. no need to clean')
        return

    outputs = get_outputs(dirname, prefix)
    out = outputs['out']
    scf_fail = (prefix == 'scf') and (check_pw_output(out) == 'failed')
    nscf_fail = (prefix == 'nscf') and (check_pw_output(out) == 'failed')
    p2q_fail = (prefix == 'p2q') and (check_pw_output(out) == 'failed')
    jastrow_fail = (prefix == 'jastrow') and (check_jastrow_output(out) == 'failed')
    dmc_fail = ('dmc' in prefix) and (check_dmc_output(out) == 'failed')
    if force_clean or (scf_fail or nscf_fail or p2q_fail or jastrow_fail or dmc_fail):
        clean(dirname, prefix)
    else:
        print(f'{dirname} {prefix} is finished! no need to clean')

def check_outputs(dirname, timesteps, force_clean=False):
    '''
    Wrapper of `check_output`.
    Inspects the directory, checks the calculation status, and cleans accordingly.
    `force_clean` if True will remove old results regardless of the status.
    '''
    things_to_check = ['scf', 'nscf', 'p2q', 'jastrow'] + [f'dmc{timestep}' for timestep in timesteps]
    for prefix in things_to_check:

        # use the given dirname for an absolute path
        if dirname.startswith('/'):
            dir_path = dirname

        # prepend `runs` to the dirname for a relative path becasuse Nexus does so
        else:
            dir_path = f'runs/{dirname}'

        check_output(dir_path, prefix, force_clean=force_clean)

def test():
    done = 'test_jastrow/disregistry0.0_ntiling2_d3.5_qmc_nk3_opt_methodbasic_opt_samples50000/'
    check_outputs(done, [0.02])

    new = 'test_jastrow/disregistry0.0_ntiling2_d3.5_qmc_nk3_opt_methodbasic_opt_samples2000/'
    check_outputs(new, [0.02])

    failed = 'test_jastrow/disregistry0.0_ntiling2_d3.5_qmc_nk3_opt_methodbasic_opt_samples20000/'
    check_outputs(failed, [0.02])


if __name__ == '__main__':
    test()
