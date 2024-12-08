"""
Project: SHARK-XT Rates and ADE
Author: Kristína Kocáková
Description:
1. In case PyRate output mcmm.log files weren't formatted correctly (happens when using Cluster)
remove rows with incorrect formatting
2. In case numbers in speciation or extinction rates get formatted incorrectly (two decimal markers in one number), drop rows containing these
This can be done on the combined file as well.
3. Add taxa names to the output file produced by the -ginput function in PyRate
"""

from pandas import *
import re

# 1.
file_name = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
def rm_lines(file_name):

    # read the file and skip lines which have incorrect formatting (i.e. too many columns)
    log = read_csv(
        "/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/all_sp_fast_burn/all_species_15Myr_{i}_Grj_mcmc.log".format( i=file_name), sep="\t", header=None, on_bad_lines='skip')

    # remove lines which are too short (i.e. have nan in them)

    log = log.dropna()

    # save the modified file

    log.to_csv("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/all_sp_fast_burn/modded/all_species_15Myr_{i}_Grj_mcmc.log".format(i = file_name), sep="\t", header=None, index = None)

for i in file_name:
    rm_lines(i)

# 2.

log = open("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/all_sp_fast_burn/combined_10_sp_rates.log", "r")

log_lines = log.readlines()

log_lines_clean = []

for i in log_lines:
    if re.search(r'[0-9]*\.[0-9]*\.[0-9]', i):
        continue
    else:
        log_lines_clean.append(i)

with open("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/all_sp_fast_burn/mod_combined_10_sp_rates.log", "w") as f:
    for line in log_lines_clean:
        f.write(f"{line}")


# 3.

from pandas import *

te_ts = read_csv("/path/to/_se_est.txt", sep="\t")

species = read_csv("/path/to/_mcmc.log", sep="\t", )
columns = species.columns.to_list()  # put the column headers in a list

species_names = []  # empty list

# select only column names that contain "_TS", which will select only taxa names
for i in columns:
    if "TS" in i:
        species_names.append(i)

# remove "_TS" from the names
for i in range(len(species_names)):
    species_names[i] = species_names[i].replace("_TS", "")

te_ts["species"] = species_names  # add the species names list as a column to the se_ext.txt file

# write the new file
te_ts.to_csv("/path/to/output.txt", sep="\t")


















