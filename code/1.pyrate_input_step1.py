"""
Project: SHARK-XT Rates and ADE
Author: Kristína Kocáková
Description:
Create an input file for PyRate program from an .xlsx file, options for species and genus level
If adenn = True the file will simply not contain the collection number, which is required for ADE-NN input
"""

from pandas import *


def pyrate_input(path_to_database, taxonomic_rank, path_to_output, adenn = False):
    database = ExcelFile(path_to_database)
    occurrences = read_excel(database, "Occurrences")


    # SELECT ONLY THE VALID ENTRIES (considering both age and taxonomy)

    occurrences = occurrences.loc[occurrences["age_evaluation"] == "valid"]
    occurrences = occurrences.loc[occurrences["taxonomy_validation"] == "valid"]

    #CREATE AN INPUT FILE FOR SPECIES
    if taxonomic_rank == "species":
        occurrences = occurrences.loc[occurrences["rank"] == "species"]
        occurrences = occurrences.loc[occurrences["early_interval"] != "present"]

        #APPLY ANY ADDITIONAL FILTERS (Optional), e.g.:

        #Only occurrences with age resolution below 15 Myr
        occurrences = occurrences.loc[occurrences["age_range"] <= 15]

        #Select only required columns

        occurrences = occurrences[["accepted_name", "status", "max_ma", "min_ma", "collection_no"]]


        occurrences = occurrences.rename(columns={"accepted_name": "taxon_name"})

        #DROP DUPLICATES (Optional, occurrences of the same taxon from the same locality with identical ages could represent the same individual or population)
        occurrences = occurrences.drop_duplicates(keep="first")

        if adenn == True:
            occurrences = occurrences[["taxon_name", "status", "max_ma", "min_ma"]]

        #WRITE OUTPUT FILE
        occurrences.to_csv(path_to_output, sep="\t", index=False)

    #CREATE AN INPUT FILE FOR GENERA
    #INCLUDING OCCURRENCES IDENTIFIED DO A GENUS LEVEL AND THE GENUS OF OCCURRENCES IDENTIFIED TO A SPECIES LEVEL
    if taxonomic_rank == "genus":
        occurrences = occurrences.loc[(occurrences["rank"] == "species") | (occurrences["rank"] == "genus")]
        occurrences = occurrences.loc[occurrences["early_interval"] != "present"]

        #APPLY ANY ADDITIONAL FILTERS (Optional), e.g.:

        #Only occurrences with age resolution below 15 Myr
        occurrences = occurrences.loc[occurrences["age_range"] <= 15]

        #Select only required columns
        occurrences = occurrences[["genus", "genus_status", "max_ma", "min_ma", "collection_no"]]

        occurrences = occurrences.rename(columns={"genus": "taxon_name"})

        #DROP DUPLICATES (Optional, occurrences of the same taxon from the same locality with identical ages could represent the same individual or population)
        occurrences = occurrences.drop_duplicates(keep="first")

        if adenn == True:
            occurrences = occurrences.loc[["genus", "genus_status", "max_ma", "min_ma"]]

        #WRITE OUTPUT FILE
        occurrences.to_csv(path_to_output, sep="\t", index=False)



pyrate_input("/Users/kristinakocakova/Dropbox/Analyses/Data/Master files/fins.xlsx", "species",
             "//data/all_species/all_species_input.txt")
