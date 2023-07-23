#! /bin/bash
#SBATCH --job-name="fit_bootstrap/QMC_17_kT0.004"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20
#SBATCH --output=fit_bootstrap/QMC_17_kT0.004/kc.out

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
export ASE_LAMMPSRUN_COMMAND=/home/krongch2/projects/lammps/lammps/src/lmp_mpi
python -u kc.py --method QMC -kT 0.004 --model_id 17
