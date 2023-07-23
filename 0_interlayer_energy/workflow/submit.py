'''
Submit actual production runs.
Relavant parameters should be carefully studied and converged before running this file
'''
import numpy as np

import run

def get_disregistry(stacking):
    return {
        'ab': 0.0,
        'sp': 0.16667,
        'mid': 0.5,
        'aa': 0.66667
        }[stacking]

def submit_production(stacking, ntiling, timesteps=[0.02]):
    params = {
        'disregistry': [get_disregistry(stacking)],
        'ntiling': [ntiling],
        'd': [3.0, 3.2, 3.35, 3.5, 3.65, 3.8, 4.0, 4.5, 5.0, 6.0, 7.0],
        'qmc_nk': [4]
        }
    run.run_workflows(params, timesteps=timesteps, data_path='production')

if __name__ == '__main__':
    submit_production('ab', 6, timesteps=[0.01, 0.02])
