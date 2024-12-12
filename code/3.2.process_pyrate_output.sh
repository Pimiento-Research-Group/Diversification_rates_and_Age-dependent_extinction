: '
Project: Diversification Rates and ADE
Author: Kristína Kocáková
Description:
Processing of PyRate .log outputs
1. Combining 10 replicates into a single set of posterior values
2. Generating of the input file for plotting of the rates
3. Generating a file containing estimated speciation and extinction times per taxon
'
#1.
python3 ./PyRate-master/PyRate.py -combLogRJ /Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/all_species/
#2.
python3 ./PyRate-master/PyRate.py -plotRJ /Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/all_species/ -tag combined -grid_plot 0.1 -b 0.1
#3.
python3 ./PyRate-master/PyRate.py -ginput ./PyRate_Analysis/outputs/2024/fast_burnin_10_2024/all_species/combined_10_mcmc.log



