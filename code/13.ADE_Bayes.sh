: '
Project: Diversification Rates and ADE
Author: Kristína Kocáková
Description:
Estimation of the Weibull distribution shape parameter using the ADE-Bayes method
-filter_lad argument needs to be changed to the appropriate value to represent the time bin of interest
This script was written to be launched in a HPC Cluster, but can be run locally by simply using the last line withouth the srun command
'

#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.kocakova@pim.uzh.ch
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=48:00:00
#SBATCH --output=ade_bayes_1.out

module load anaconda3
source activate pyrate5


srun python3 ~/data/PyRate-master/PyRate_2.py ~/data/PyRate_Data_3/ADE/2024/ADE1/all_species_PyRate.py -ADE 1 -qShift ~/data/PyRate_Data_3/ages.txt -filter_lad 145 95.61
