: '
Project: Diversification Rates and ADE
Description:
Processing of PyRate .log outputs
1. Combining 10 replicates into a single set of posterior values
2. Generating of the input file for plotting of the rates
3. Generating a file containing estimated speciation and extinction times per taxon
4. Generate file containing the preservation rate through time and plot it
'
#1.
python3 ./PyRate-master/PyRate.py -combLogRJ /path_to_output_folder/
#2.
python3 ./PyRate-master/PyRate.py -plotRJ /path_to_output_folder/ -tag combined -grid_plot 0.1 -b 0.1
#3.
python3 ./PyRate-master/PyRate.py -ginput /path_to_output_folder/combined_10_mcmc.log
#4
python3 ./PyRate-master/PyRate.py -plotQ /path_to_output_folder/species_1_Grj_mcmc.log -qShift /path_to_input_files/stages.txt -b 0.1
# to plot (might have adjust the path within the .r script for an absolute path)
Rscript /path_to_output_folder/species_1_Grj_mcmc_RTT_Qrates.r



