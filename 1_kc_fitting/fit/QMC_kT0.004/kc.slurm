#! /bin/bash
#SBATCH --job-name="fit/QMC_kT0.004"
#SBATCH --time=24:00:00
#SBATCH --partition="qmchamm"
#SBATCH --cpus-per-task=20
#SBATCH --output=fit/QMC_kT0.004/kc.out

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
export ASE_LAMMPSRUN_COMMAND=/home/krongch2/projects/lammps/lammps/src/lmp_mpi
python -u kc.py --method QMC -kT 0.004 --model_id 2
