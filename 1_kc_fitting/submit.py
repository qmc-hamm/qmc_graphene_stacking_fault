'''
Generates submission scripts and submit jobs on cluster to fit the KC potential
'''
import os

def get_hours(queue):
    hours_map = {'secondary': 4, 'qmchamm': 4, 'physics': 14, 'test': 1}
    return hours_map[queue]

def submit_slurm(dirname, cmd, queue='secondary'):
    submit_fn = f'{dirname}/kc.slurm'
    hours = get_hours(queue)

    with open(submit_fn, 'w') as f:
        f.write(
f'''#! /bin/bash
#SBATCH --job-name="{dirname}"
#SBATCH --time={hours}:00:00
#SBATCH --partition="{queue}"
#SBATCH --cpus-per-task=20
#SBATCH --output={dirname}/kc.out

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
export ASE_LAMMPSRUN_COMMAND=/home/krongch2/projects/lammps/lammps/src/lmp_mpi
{cmd}
''')
    os.system(f'sbatch {submit_fn}')

def submit_windows(dirname, cmd):
    os.system(f'{cmd} > {dirname}/kc.out')

def submit(method, kT, model_id=0, queue='secondary'):
    if kT == 'inf':
        kT_str = 'inf'

    else:
        kT = float(kT)
        kT_str = f'{kT:.3f}'

    dirname = f'fit_bootstrap/{method}_{int(model_id):02}_kT' + kT_str
    cmd = f'python -u kc.py --method {method} -kT {kT_str} --model_id {model_id}'
    workdir = os.getcwd()
    os.makedirs(dirname, exist_ok=True)

    submit_slurm(dirname, cmd, queue=queue)

def fit_multi_kT():
    for method in [
        'QMC',
        'DFT-D2',
        'DFT-D3',
        'DFT-MBD'
        ]:
        for kT in [
            'inf',
            0.002,
            0.003,
            0.004,
            0.005,
            0.006,
            0.007,
            0.008,
            0.009,
            0.010,
            0.011,
            0.012,
            0.013,
            0.014,
            0.015,
            0.016,
            0.017,
            0.018,
            0.019,
            0.020,
            ]:
            submit(method=method, kT=kT, queue='secondary')

def fit_bootstrap(method='QMC'):
    for model_id in range(20):
        for kT in [
            'inf',
            0.004
            ]:
            submit(method=method, kT=kT, model_id=model_id, queue='secondary')

def submit_gen_kc():
    cmd = f'python -u bootstrap_kc.py'
    submit_slurm('.', cmd, queue='qmchamm')

if __name__ == '__main__':
    submit_gen_kc()
