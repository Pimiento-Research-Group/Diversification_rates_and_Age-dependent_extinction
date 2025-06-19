#
#Project: Diversification Rates and ADE
#Description:
#Produce a final version of PyRate input file using the initial .txt file created in step 1
#
f <- "path/to/PyRate-master/pyrate_utilities.r"
path <- getwd()

source(file = f)

extract.ages(file = paste(path, "/data/all_species_input.txt", sep = ""), replicates=10, random=TRUE)
