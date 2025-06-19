: '
Project: Diversification Rates and ADE
Description:
Calculation of extinction, speciation, and net diversification rater using PyRate.
This script was written to be launched in a HPC Cluster, but can be run locally by simply using the last line withouth the srun command
and by providing a number from 1 to 10  in the -j argument (the maximum number depends on the number of replicates done when generating the PyRate input file)
'
#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.kocakova@pim.uzh.ch
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --output=fast_burnin.out
#SBATCH --array=1-10

module load anaconda3
source activate pyrate5

srun python3 ~/data/PyRate-master/PyRate.py ~/data/PyRate_Data_3/Spec_Ext/all_sp_15Myr_fast_burnin/all_species_15Myr_PyRate.py  -n 50000000 -s 10000 -j $SLURM_ARRAY_TASK_ID -qShift ~/data/PyRate_Data_3/ages.txt -min_dt 0.5 -pP 2 0 -singleton 1 -mG -fast_burnin 25000 -fQ 0.05 -fU 0.05 0.2 0

