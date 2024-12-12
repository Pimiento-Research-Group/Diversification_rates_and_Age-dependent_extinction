#
#Project: Diversification Rates and ADE
#Author: Kristína Kocáková
#Description:
#Produce a final version of PyRate input file using the initial .txt file created in step 1
#
f <- "path/to/PyRate-master/pyrate_utilities.r"

source(file = f)

extract.ages(file = "/Users/kristinakocakova/PycharmProjects/SHARKS-XT-Diversification_ADE/data/all_species/all_species_input.txt", replicates=10, random=TRUE)
