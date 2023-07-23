#! /bin/bash
#SBATCH --job-name="fit_bootstrap/QMC_00_kTinf"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20
#SBATCH --output=fit_bootstrap/QMC_00_kTinf/kc.out

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
export ASE_LAMMPSRUN_COMMAND=/home/krongch2/projects/lammps/lammps/src/lmp_mpi
python -u kc.py --method QMC -kT inf --model_id 0
