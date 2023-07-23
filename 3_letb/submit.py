import os

def get_resources_slurm(twist_angle):
    if twist_angle in ['4-4', '6-0']:
        time = '00:20:00'
        queue = 'secondary'
    else:
        time = '48:00:00'
        queue = 'qmchamm'
    return time, queue

def get_resources_lsf(twist_angle):
    if twist_angle in ['4-4', '6-0']:
        time = '00:10'
        queue = 'batch'
    elif twist_angle in ['0-93', '0-84', '0-5']:
        time = '24:00'
        queue = 'batch-hm'
    else:
        time = '6:00'
        queue = 'batch-hm'
    return time, queue

def submit_slurm(submit_fn, job_name, time, queue, output, cmd):
    with open(submit_fn, 'w') as f:
        f.write(
f'''#!/bin/bash
#SBATCH --job-name="{job_name}"
#SBATCH --time={time}
#SBATCH --partition="{queue}"
#SBATCH --cpus-per-task=20
#SBATCH --output={output}

cd $SLURM_SUBMIT_DIR
{cmd}
''')
    os.system(f'sbatch {submit_fn}')

def submit_lsf(submit_fn, job_name, time, queue, output, cmd):
    with open(submit_fn, 'w') as f:
        f.write(
f'''#!/bin/bash
#BSUB -P mat221
#BSUB -W {time}
#BSUB -nnodes 1
#BSUB -J {job_name}
#BSUB -o {output}
#BSUB -q {queue}

eval "$(conda shell.bash hook)"
conda activate my_ase_env
module load gcc/9.3.0 python/3.8-anaconda3 netlib-lapack netlib-scalapack cuda/11.0.3 fftw
export OMP_NUM_THREADS=1
jsrun -n 1 -c 42 -a 1 --bind="none" {cmd}
''')
    os.system(f'bsub {submit_fn}')

def submit_lsf_new(submit_fn, job_name, time, queue, output, cmd):
    with open(submit_fn, 'w') as f:
        f.write(
f'''#!/bin/bash
#BSUB -P mat221
#BSUB -W {time}
#BSUB -nnodes 1
#BSUB -J {job_name}
#BSUB -o {output}
#BSUB -q {queue}

export OMP_NUM_THREADS=1
module load gcc/9.3.0 python/3.8-anaconda3 netlib-lapack netlib-scalapack cuda/11.0.3 fftw
source activate my_ase_env
jsrun -n 1 -c 42 -a 1 --bind="none" {cmd}
''')
    os.system(f'bsub -L $SHELL {submit_fn}')

def submit(machine, twist_angle, pot, relax_keyword):
    dirname = f'bands/{twist_angle}_{pot}_{relax_keyword}'
    os.makedirs(dirname, exist_ok=True)

    job_name = f'{twist_angle}_{pot}'
    output = f'{dirname}/solve.out'
    cmd = f'python -u solve.py {twist_angle} {pot} {relax_keyword}'

    if machine == 'slurm':
        submit_fn = f'{dirname}/solve.slurm'
        time, queue = get_resources_slurm(twist_angle)
        submit_slurm(submit_fn, job_name, time, queue, output, cmd)
    elif machine == 'lsf':
        submit_fn = f'{dirname}/solve.lsf'
        time, queue = get_resources_lsf(twist_angle)
        submit_lsf(submit_fn, job_name, time, queue, output, cmd)

for twist_angle in [
    '0-84',
    ]:
    for pot in [
        'qmc',
        'ouyang',
        'dft_d2',
        'dft_d3'
        ]:
        submit('lsf', twist_angle, pot, 'rigid')
